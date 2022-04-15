import os
import coverage
import encode_utils


# generate an experiment list file w/ paths of bed files to be merged
data_path = './data'
metadata_path = './data/metadata.tsv'
exp_list_bed_path = os.path.join(data_path, "exp_list_bed.tsv")  # 2-column list of experiments, eg. (name    dirpath)
criteria={'File assembly': "GRCh38",                             # criteria for filtering files to process
          'File format': 'bigWig',
          'Output type': 'fold change over control',
          'Biosample term name': 'GM12878',     
          'Biological replicate(s)': '1, 2',        
          }
label_list=['Experiment target', 'Biosample term name', 'Experiment accession']  # used to name experiment
summary = encode_utils.generate_exp_list_from_metadata(data_path, metadata_path, exp_list_bed_path, 
                                                criteria, label_list, download=True, ext='bigWig')


# generate whole chromosome bed file for train, valid, test -- tiled based on stride
bin_size = 2048
prefix_path = './data/name' 
chrom_sizes_path = '../genomes/hg38.chrom.sizes'                 # path to chrom sizes
genome_path = './hg38.fa'
blacklist_path = './hg38-blacklist.v2.bed.gz'
coverage.whole_chrom_bed(bin_size, prefix_path, chrom_sizes_path, blacklist_path,
                         stride=None, valid_chr='8,9', test_chr='11,12', 
                         ignore_chr_y=True, ignore_auxiliary_chr=True)


# generate h5 with one-hot inputs and coverage targets 
bigwig_paths = [s[1] for s in summary]      # get list of bigwig_paths
coverage.process_data_h5(prefix_path, bigwig_paths, genome_path, alphabet='ACGT', 
                         uncertain_N=True, standard_coords=False)

#coverage.process_data_tfr(prefix_path, bigwig_paths, genome_path, alphabet='ACGT', 
#                          uncertain_N=True, standard_coord=False)