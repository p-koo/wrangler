import os, h5py, gzip
import subprocess
import numpy as np
import pandas as pd
import utils


def process_data(
  pos_path,
  neg_path,
  prefix_save_path='./sample',
  genome_path='./hg38.fa',
  seq_len=200,
  max_len_thresh=300,
  filter_N=False,
  uncertain_N=True,
  gc_match=True,
  neg_pos_ratio=1,
  valid_chr='6,7',
  test_chr='4,8',
  alphabet="ACGT",
  blacklist_path=None,
  standard_coord=False
  seed=None,
):
  """Preprocess data for a single-class task."""
  print(prefix_save_path)
  print("Processing positive label data")

  # set random number seed
  if seed is not None:
    np.random.seed(seed)

  # remove extremely large peaks
  tf_filtered_path = prefix_save_path + "_pos_filtered.bed"
  filter_max_length(pos_path, tf_filtered_path, max_len_thresh)

  # create new bed file with seq_len enforced
  pos_bed_path = prefix_save_path + "_pos_" + str(seq_len) + ".bed"
  enforce_constant_size(tf_filtered_path, pos_bed_path, seq_len)
  
  # remove ENCODE blacklist sites if provided
  if blacklist_path:
    print("removing unmappable regions defined by: %s"%(blacklist_path))

    utils.bedtools_intersect(
      a=pos_bed_path,
      b=blacklist_path,
      output_path=pos_bed_path, 
      write_a=False, nonoverlap=True
    )

  # get coordinates
  pos_coords = get_bed_coords(pos_bed_path, standard_coord=False)

  # extract sequences from bed file and save to fasta file
  pos_fasta_path = prefix_save_path + "_pos.fa"
  utils.bedtools_getfasta(
    pos_bed_path, genome_path, output_path=pos_fasta_path, strand=True
  )

  # remove ENCODE blacklist sites if provided
  if blacklist_path:
    print("removing unmappable regions defined by: %s"%(blacklist_path))

    utils.bedtools_intersect(
      a=pos_bed_path,
      b=blacklist_path,
      output_path=pos_bed_path, 
      write_a=False, nonoverlap=True
    )

  # parse sequence from fasta file
  pos_seq,_ = utils.parse_fasta(pos_fasta_path)

  if filter_N:
    N = len(pos_seq)
    # filter sequences with absent nucleotides
    pos_seq, good_index = filter_nonsense_sequences(pos_seq)
    pos_coords = pos_coords[good_index]
    print('    Filtered %d positive sequences with N'%(N-len(pos_seq)))

  # convert filtered sequences to one-hot representation
  pos_one_hot = utils.convert_one_hot(pos_seq, alphabet, uncertain_N=uncertain_N)


  print("Processing negative label data")

  # get non-overlap between pos peaks and neg peaks
  neg_bed_path = prefix_save_path + "_nonoverlap.bed"
  utils.bedtools_intersect(
    neg_path, pos_path, neg_bed_path, write_a=True, nonoverlap=True
  )

  # create new bed file with seq_len enforced
  neg_bed_path2 = prefix_save_path + "_neg_" + str(seq_len) + ".bed"
  enforce_constant_size(neg_bed_path, neg_bed_path2, seq_len)
  
  # remove ENCODE blacklist sites if provided
  if blacklist_path:
    print("removing unmappable regions defined by: %s"%(blacklist_path))

    utils.bedtools_intersect(
      a=neg_bed_path2,
      b=blacklist_path,
      output_path=neg_bed_path2, 
      write_a=False, nonoverlap=True
    )

  # get coordinates
  neg_coords = get_bed_coords(neg_bed_path2, standard_coord=False)

  # extract sequences from bed file and save to fasta file
  neg_fasta_path = prefix_save_path + "_neg.fa"
  utils.bedtools_getfasta(
    neg_bed_path2, genome_path, output_path=neg_fasta_path, strand=True
  )

  # parse sequence and chromosome from fasta file
  neg_seq,_ = utils.parse_fasta(neg_fasta_path)

  if filter_N:
    # filter sequences with absent nucleotides
    N = len(neg_seq)
    neg_seq, good_index = filter_nonsense_sequences(neg_seq)
    neg_coords = neg_coords[good_index]
    print('    Filtered %d negative sequences with N'%(N-len(neg_seq)))

  # convert filtered sequences to one-hot representation
  neg_one_hot = utils.convert_one_hot(neg_seq, alphabet, uncertain_N=uncertain_N)

  # match GC content
  if gc_match:
    print("    Matching GC content between negative and positive data.")
    neg_one_hot, match_index = match_gc_content(pos_one_hot, neg_one_hot, neg_pos_ratio=neg_pos_ratio)
    neg_coords = neg_coords[match_index]

  else:
    num_pos = len(pos_one_hot)
    num_neg = len(neg_one_hot)
    if num_pos < num_neg:
      num_keep = int(num_pos * neg_pos_ratio)
      index = np.random.permutation(len(neg_one_hot))[:num_keep]
      neg_one_hot = neg_one_hot[index]
      neg_coords = neg_coords[index]
      print("    Downsampling negative set to %d"%(num_keep))

  print("%d positive sequences and %d negative sequences"%(len(pos_one_hot), len(neg_one_hot)))

  # merge positive and negative labels
  one_hot = np.vstack([pos_one_hot, neg_one_hot])
  labels = np.vstack(
    [np.ones((len(pos_one_hot), 1)), np.zeros((len(neg_one_hot), 1))]
  )
  coords = np.concatenate([pos_coords, neg_coords])

  # split data into training, validation, and test sets
  print("Splittig dataset into:")
  print("    validation: %s"%(valid_chr))
  print("    test: %s"%(test_chr))
  train, valid, test = utils.split_dataset_chrom(
    one_hot, labels, coords, valid_chr, test_chr
  )

  # save to hdf5 file
  file_path = prefix_save_path + ".h5"
  print("Saving to: %s"%(file_path))
  with h5py.File(file_path, "w") as fout:
    fout.create_dataset("x_train", data=train[0], compression="gzip")
    fout.create_dataset("y_train", data=train[1], compression="gzip")
    fout.create_dataset("coords_train", data=train[2].astype("S"), compression="gzip")
    fout.create_dataset("x_valid", data=valid[0], compression="gzip")
    fout.create_dataset("y_valid", data=valid[1], compression="gzip")
    fout.create_dataset("coords_valid", data=valid[2].astype("S"), compression="gzip")
    fout.create_dataset("x_test", data=test[0], compression="gzip")
    fout.create_dataset("y_test", data=test[1], compression="gzip")
    fout.create_dataset("coords_test", data=test[2].astype("S"), compression="gzip")
  print('%d train sequences'%(len(train[0])))
  print('%d valid sequences'%(len(valid[0])))
  print('%d test sequences '%(len(test[0])))

  # clean up directory
  os.remove(tf_filtered_path)
  os.remove(pos_bed_path)
  os.remove(pos_fasta_path)

  os.remove(neg_bed_path)
  os.remove(neg_bed_path2)
  os.remove(neg_fasta_path)





#------------------------------------------------------------
# Useful functions for single-task
#------------------------------------------------------------



def filter_max_length(bed_path, output_path, max_len=1000):
  """Function to plot histogram of bed file peak sizes; automatically infers
  compression  from extension and allows for user-input in removing outlier sequence
  Parameters
  -----------
  bed_path : <str>
      Path to bed file.
  output_path : <int>
      Path to filtered bed file.
  max_len: <int>
      Cutoff for maximum length of peak -- anything above will be filtered out.
  Returns
  ----------
  None
  """

  # check if bedfile is compressed
  compression = "gzip" if _is_gzipped(bed_path) else None

  # load bed file
  df = pd.read_table(bed_path, header=None, sep="\t", compression=compression)
  start = df[1].to_numpy()
  end = df[2].to_numpy()

  # get peak sizes
  peak_sizes = end - start
  good_index = []
  for i, size in enumerate(peak_sizes):
    if size < max_len:
      good_index.append(i)

  # create dictionary for dataframe with filtered peaks
  data = {}
  for i in range(len(df.columns)):
    data[i] = df[i].to_numpy()[good_index]

  # create new dataframe
  df_new = pd.DataFrame(data)

  # save dataframe with fixed width window size to a bed file
  df_new.to_csv(output_path, sep="\t", header=None, index=False)



def enforce_constant_size(bed_path, output_path, window):
  """Generate a bed file where all peaks have same size centered on original peak
  Parameters
  ----------
  bed_path : <str>
      The path to the unfiltered bed file downloaded from
      ENCODE.
  output_path : <str>
      Filtered bed file name with peak size clipped to specified
      window size.
  window : <int>
      Size in bps which to clip the peak sequences.
  Returns
  ----------
  None
  Example
  ----------
  >>> window = 200
  >>> bed_path = './ENCFF252PLM.bed.gz'
  >>> output_path = './pos_'+str(window)+'.bed'
  >>> enforce_constant_size(pos_path, pos_bed_path, window, compression='gzip')
  """
  assert isinstance(window, int) and window > 0, "Enter positive integer window size."
  assert os.path.exists(bed_path), "No such bed file."

  # set up the compression argument
  if bed_path.split(".")[-1] == "gz" or bed_path.split(".")[-1] == "gzip":
    compression = "gzip"
  else:
    compression = None

  # load bed file
  df = pd.read_table(bed_path, header=None, sep="\t", compression=compression)
  # chrom = df[0].to_numpy().astype(str)  # TODO: unused variable
  start = df[1].to_numpy()
  end = df[2].to_numpy()

  # calculate center point and create dataframe
  middle = np.round((start + end) / 2).astype(int)
  left_window = np.round(window / 2).astype(int)
  right_window = window - left_window

  # calculate new start and end points
  start = middle - left_window
  end = middle + right_window

  # create dictionary for dataframe
  data = {}
  for i in range(len(df.columns)):
    data[i] = df[i].to_numpy()
  data[1] = start
  data[2] = end
  for i in range(3, df.shape[1]):
    data[i] = df[i].to_numpy()

  # create new dataframe
  df_new = pd.DataFrame(data)

  # save dataframe with fixed width window size to a bed file
  df_new.to_csv(output_path, sep="\t", header=None, index=False)



def get_bed_coords(bed_path, standard_coord=True):
  assert os.path.exists(bed_path), "No such bed file."

  # set up the compression argument
  if bed_path.split(".")[-1] == "gz" or bed_path.split(".")[-1] == "gzip":
    compression = "gzip"
  else:
    compression = None

  # load bed file
  df = pd.read_table(bed_path, header=None, sep="\t", compression=compression)
  chrom = df[0].to_numpy().astype(str)  # TODO: unused variable
  start = df[1].to_numpy().astype(str)
  end = df[2].to_numpy().astype(str)
  strand = df[5].to_numpy().astype(str)

  coords = []
  for a in zip(chrom, start, end, strand):
    if standard_coord:
      coords.append("%s:%s-%s(%s)" % (a[0], a[1], a[2], a[3]))
    else:
      coords.append("%s_%s_%s" % (a[0], a[1], a[2]))
  return np.array(coords)


def filter_nonsense_sequences(sequences):
  """Filter sequences with N.
  Parameters
  ----------
  sequences : <numpy.ndarray>
      A numpy vector of sequence strings.
  Returns
  -------
  filter_sequences : <numpy.ndarray>
      The parsed sequences from the input fasta file as a numpy array of sequences.
  good_index : <numpy.ndarray>
      A numpy array of indices corresponding to sequences without nonsense 'N'
      entries.
  Example
  -------
  >>> print(sequences)
  ['GGCTGAAATGGCCACTGGAA' 'ACGCTCTCTCATCAAGTGGT' 'GCAGAANANCGAACACCAAC'
  'NNCNNCANCNACNGGGGAAC' 'GCCTAGTCCAGACATAATTC']
  >>> print(filter_nonsense_sequences(sequences))
  (array(['GGCTGAAATGGCCACTGGAA', 'ACGCTCTCTCATCAAGTGGT',
          'GCCTAGTCCAGACATAATTC'], dtype='<U20'), array([0, 1, 4]))
  """

  # filter sequences if contains at least one 'N' character
  good_index = []
  filter_sequences = []
  for i, seq in enumerate(sequences):
    if "N" not in seq.upper():
      good_index.append(i)
      filter_sequences.append(seq)
  return np.array(filter_sequences), np.array(good_index)



def match_gc_content(pos_one_hot, neg_one_hot, neg_pos_ratio=1):
  """match the GC content of negative set (which is larger) as positive set"""
  N, L, A = pos_one_hot.shape
  gc_pos = np.sum(np.sum(pos_one_hot[:,:,[1,2]], axis=2), axis=1)/L
  gc_neg = np.sum(np.sum(neg_one_hot[:,:,[1,2]], axis=2), axis=1)/L
  print('    Average GC content for positive sequences: %.3f'%(np.mean(gc_pos)))
  print('    Average GC content for negative sequences: %.3f'%(np.mean(gc_neg)))

  pos_index = np.argsort(gc_pos)
  neg_index = np.argsort(gc_neg)
  num_neg = len(neg_index)
  num_pos = len(pos_index)

  match_index = []
  if num_neg > num_pos:
    k = 0
    status = True
    for i in pos_index:
      for j in range(k, num_neg):
        if gc_pos[i] < gc_neg[neg_index[j]]:
          if k > num_neg:
            status = False
            break
          else:
            # print("%.2f vs %.2f"%(gc_pos[i], gc_neg[neg_index[j]]))
            match_index.append(neg_index[j])
            k = j+1
            break
      if not status:
        break

  remainder = int(num_pos*neg_pos_ratio) - len(match_index)
  print('    Found %d GC-matched sequences.'%(len(match_index)))
  if remainder > 0:
    print('    Adding %d more random negative sequences.'%(remainder))
    remain_index = np.array(list(set(range(num_neg)) - set(match_index)))
    index = np.random.permutation(len(remain_index))[:remainder]     
    # index = np.argsort(gc_neg[remain_index])[::-1]
    for n in remain_index[index[:remainder]]:
      match_index.append(n)
              
  match_index = np.array(match_index)
  print('    Average GC content for sub-sampled negative sequences: %.3f'%(np.mean(gc_neg[match_index])))

  return neg_one_hot[match_index], match_index

