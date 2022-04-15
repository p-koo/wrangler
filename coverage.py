import os, sys, h5py
import numpy as np
import pandas as pd
import tensorflow as tf
import utils


def whole_chrom_bed(bin_size, prefix_path, chrom_sizes_path, 
                    blacklist_path=None, stride=None, 
                    valid_chr='11,12', test_chr='7,8', 
                    ignore_chr_y=True, ignore_auxiliary_chr=True,
                    verbose=True):
  if not stride:   # if no stride, then make non-overlapping bins
    stride = bin_size

  # convert list of chromosome numbers to proper name
  valid_chr = ['chr'+chr_index for chr_index in valid_chr.split(',')]
  test_chr  = ['chr'+chr_index for chr_index in test_chr.split(',')]

  # get chrom sizes
  chrom_sizes = utils.get_chrom_sizes(chrom_sizes_path)

  # filter list of chromosomes 
  keep_chr = utils.filter_chr_list(chrom_sizes, ignore_chr_y, ignore_auxiliary_chr)

  if verbose:
    print('Generating binned bed files.')
 
  # create files to save
  bed_files = {}
  bed_files['train'] = open(prefix_path+'_train_all.bed', 'w')
  bed_files['valid'] = open(prefix_path+'_valid_all.bed', 'w')
  bed_files['test']  = open(prefix_path+'_test_all.bed', 'w')
  
  # create bed file with non-overlapping bins
  for chr in keep_chr:
    # select dataset to save bed coordinates
    if chr in valid_chr:
      fout = bed_files['valid'] 
    elif chr in test_chr:
      fout = bed_files['test']
    else:
      fout = bed_files['train']
    
    # write bed coordinates for chromosome
    MAX = chrom_sizes[chr]
    num_bins = (MAX - bin_size)//stride
    for i in range(num_bins):
      fout.write("%s\t%d\t%d\n"%(chr, i*stride, i*stride+bin_size))

  # close bed files
  for key in bed_files.keys():
    bed_files[key].close()

  # remove ENCODE blacklist sites if provided
  if blacklist_path:
    if verbose:
      print("Removing unmappable regions defined by: %s"%(blacklist_path))
    utils.filter_blacklist(prefix_path+'_train_all.bed', prefix_path+'_train.bed', blacklist_path)
    utils.filter_blacklist(prefix_path+'_valid_all.bed', prefix_path+'_valid.bed', blacklist_path)
    utils.filter_blacklist(prefix_path+'_test_all.bed' , prefix_path+'_test.bed' , blacklist_path)
  else:
    sys.command('mv %s %s'%(prefix_path+'_train_all.bed', prefix_path+'_train.bed'))
    sys.command('mv %s %s'%(prefix_path+'_valid_all.bed', prefix_path+'_valid.bed'))
    sys.command('mv %s %s'%(prefix_path+'_test_all.bed', prefix_path+'_test.bed'))
  
  # remove unnecessary files
  os.remove(prefix_path+'_train_all.bed') 
  os.remove(prefix_path+'_valid_all.bed') 
  os.remove(prefix_path+'_test_all.bed')  



def process_data_h5(prefix_path, bigwig_paths, genome_path, alphabet='ACGT', 
                    uncertain_N=True, batch_size=512, verbose=True):

  set_names = ['train', 'valid', 'test']

  # load bigwig files
  bigwig_files = []
  for bigwig_path in bigwig_paths:   
    bigwig_files.append(pyBigWig.open(bigwig_path))
  num_targets = len(bigwig_files)

  # process datasets
  with h5py.File(prefix_path+'.h5', 'w') as fout:
    for set_name in set_names:  # loop through data splits
      if verbose:
        print('Processing %s set'%(set_name))

      #--------------------------------------------------------------------------
      # process input sequences
      #--------------------------------------------------------------------------
      if verbose:
        print('  Processing inputs of %s set'%(set_name))

      # generate fasta file from bed coordinates
      fasta_path = prefix_path+".fa"    
      subprocess.call(
        "bedtools getfasta -fi {} -s -bed {} -fo {}".format(
          genome_path, prefix_path+'_'+set_name+'.bed', fasta_path
        ),
        shell=True,
      )

      # parse fasta files
      seqs, names = parse_fasta(fasta_path)
      #seqs = seqs[:10]
      num_seq = len(seqs)
      bin_size = len(seqs[0])
      os.remove(fasta_path)  # remove unnecessary files

      # create h5 input dataset
      dataset = fout.create_dataset('x_'+set_name, (num_seq, bin_size, len(alphabet)))

      # convert one-hot and store in h5 in batches (to maintain low memory footprint)
      num_batches = np.floor(num_seq/batch_size).astype(int)
      for i in range(num_batches):
        if verbose:
          if np.mod(i+1, 50) == 0:
            print("    batch %d out of %d"%(i+1, num_batches))
        one_hot = convert_one_hot_fast(seqs[i*batch_size:(i+1)*batch_size], alphabet=alphabet)#, uncertain_N=uncertain_N)
        dataset[i*batch_size:(i+1)*batch_size,:,:] = one_hot
      if num_seq > num_batches*batch_size:
        one_hot = convert_one_hot_fast(seqs[num_batches*batch_size:], alphabet=alphabet)#, uncertain_N=uncertain_N)
        dataset[num_batches*batch_size:,:,:] = one_hot    

      #--------------------------------------------------------------------------
      # process targets
      #--------------------------------------------------------------------------
      print('  Processing targets of %s set'%(set_name))

      # open bed file
      bed_df = pd.read_csv(prefix_path+'_'+set_name+'.bed', delimiter='\t', header=None)
      num_bed_entries = len(bed_df)

      # create h5 target dataset
      dataset = fout.create_dataset('y_'+set_name, (num_bed_entries, bin_size, num_targets))

      # go through each bed entry, get coverage value and store in h5
      counter = 0
      coverage = np.zeros((batch_size, bin_size, num_targets))
      coords = []
      for i, row in bed_df.iterrows():
        if np.mod(i+1, 25000) == 0:
          if verbose:
            print("    entry %d out of %d"%(i+1, num_bed_entries))

        # coordinates
        chr, start, end = row

        # save coordinates
        coords.append("%s:%d-%d"%(chr, start, end))

        # concatenate coverage and save
        for j, bw in enumerate(bigwig_files):
          coverage[counter,:,j] = bw.values(chr, start, end)
        counter += 1

        if counter == batch_size-1:
          dataset[i-batch_size:i,:,:] = coverage
          batch = np.minimum(num_bed_entries-i-1, batch_size)
          coverage = np.zeros((batch, bin_size, num_targets))
          counter = 0
      if counter:
        dataset[-batch:,:,:] = coverage

      coord_set = fout.create_dataset("coords_"+set_name, data=np.array(coords).astype("S"))



def process_data_tfr(prefix_path, bigwig_paths, genome_path, alphabet='ACGT', 
                     uncertain_N=True, batch_size=512, standard_coord=False, verbose=True):
 
  set_names = ['valid', 'test', 'train']

  # load bigwig files
  bigwig_files = []
  for bigwig_path in bigwig_paths:   
    bigwig_files.append(pyBigWig.open(bigwig_path))
  num_targets = len(bigwig_files)

  # process datasets
  for set_name in set_names:  # loop through data splits
    if verbose:
      print('Processing %s set'%(set_name))
    
    # open tf record_file
    tf_opts = tf.io.TFRecordOptions(compression_type='ZLIB')
    with tf.io.TFRecordWriter(prefix_path+'_'+set_name+'.tfr', tf_opts) as writer:
  
      # generate fasta file from bed coordinates
      fasta_path = prefix_path+".fa"    
      subprocess.call(
        "bedtools getfasta -fi {} -s -bed {} -fo {}".format(
          genome_path, prefix_path+'_'+set_name+'.bed', fasta_path
        ),
        shell=True,
      )

      # parse fasta files
      seqs, names = parse_fasta(fasta_path)
      bin_size = len(seqs[0])
      os.remove(fasta_path)  # remove unnecessary files

      # open bed file
      bed_df = pd.read_csv(prefix_path+'_'+set_name+'.bed', delimiter='\t', header=None)
      num_bed_entries = len(bed_df)

      for i, row in bed_df.iterrows():
        if np.mod(i+1, 25000) == 0:
          if verbose:
            print("    entry %d out of %d"%(i+1, num_bed_entries))

        # coordinates
        chr, start, end = row
        if standard_coord:
          coord = "%s:%d-%d"%(chr, start, end)
        else:            
          coord = "%s_%d_%d"%(chr, start, end)

        # one-hot input sequences
        one_hot = utils.convert_one_hot([seqs[i]], alphabet=alphabet, uncertain_N=uncertain_N)

        # concatenate coverage and save
        coverage = np.zeros((bin_size, num_targets))
        for j, bw in enumerate(bigwig_files):
          coverage[:,j] = bw.values(chr, start, end)

      example_proto = utils.serialize_example(one_hot, coverage, coord.encode())
      writer.write(example_proto)


