import os
import coverage
import utils

bin_size = 2048
stride = None
valid_chr = '11,12'
test_chr = '8,9' 
ignore_chr_y = True
ignore_auxiliary_chr = True

data_path = utils.make_directory(base_path, cell_type)           # path to data
exp_list_bed_path = os.path.join(base_path, "exp_list_bed.tsv")  # filename to save list of experiments-2columns (name    datapath)
criteria={'File assembly': "GRCh38",                                             # criteria for narrowing down files to process
          'File format': 'bigWig',
          'Output type': 'fold change over control',
          'Biosample term name': 'GM12878',     
          'Biological replicate(s)': '1, 2',        
          }
label_list=['Experiment target', 'Biosample term name', 'Experiment accession']  # used to name experiment

prefix_path = './cell_type' 
chrom_sizes_path = '../genomes/hg38.chrom.sizes'                 # path to chrom sizes
genome_path = './hg38.fa'
blacklist_path = './hg38-blacklist.v2.bed.gz'


# generate an experiment list file w/ paths of bed files to be merged
summary = utils.generate_exp_list_from_metadata(data_path, metadata_path, exp_list_bed_path, 
                                                criteria, label_list, download=True, ext='bigWig')

# get list of bigwig_paths
bigwig_paths = [s[1] for s in summary]

# generate whole chromosome bed file for train, valid, test -- tiled based on stride
coverage.whole_chrom_bed(bin_size, prefix_path, chrom_sizes_path, blacklist_path,
                         stride=stride, valid_chr=valid_chr, test_chr=test_chr, 
                         ignore_chr_y=ignore_chr_y, ignore_auxiliary_chr=ignore_auxiliary_chr)

# generate h5 with one-hot inputs and coverage targets 
coverage.process_data_h5(prefix_path, bigwig_paths, genome_path, alphabet='ACGT', 
                         uncertain_N=True, standard_coords=False)

#coverage.process_data_tfr(prefix_path, bigwig_paths, genome_path, alphabet='ACGT', 
#                          uncertain_N=True, standard_coord=False)