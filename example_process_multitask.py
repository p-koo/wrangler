import os
import utils
import encode_utils
import multitask

####################################################################################
# Generate an experiment list file w/ paths of bed files to be merged
####################################################################################

data_path = './data'                                             # path to data (will be downloaded here if not present) 
metadata_path = './data/metadata.tsv'                            # path to metadata file from ENCODE
exp_list_bed_path = os.path.join(data_path, "exp_list_bed.tsv")  # filename to save list of experiments-2columns (name    datapath)
criteria={'File assembly': "GRCh38",                             # criteria for narrowing down files to process
          'File format': 'bed narrowPeak',
          'Output type': 'IDR thresholded peaks',
          'Biosample term name': cell_type,             
          }
label_list=['Experiment target',        # used to name experiment
            'Biosample term name', 
            'Experiment accession']  
summary = encode_utils.generate_exp_list_from_metadata(data_path, metadata_path, exp_list_bed_path, 
											           criteria, label_list, download=True, ext='bed')

####################################################################################
# Generate peak-based, merged binary classification dataset -- save as h5 or tfr
####################################################################################

# generate the dataset and save to h5 file
prefix_path = os.path.join(data_path, 'dataset_name')
genome_path = '../genomes/hg38.fa'                      # path to reference genome
chrom_sizes_path = '../genomes/hg38.chrom.sizes'        # path to chrom sizes

multitask.process_data(
    exp_bed_list_path=exp_list_bed_path,    
    prefix_save_path=data_path,
    genome_path=genome_path,
    chrom_sizes_path =chrom_sizes_path,
    seq_len=1000,                        # length of input sequences
    merge_overlap_len=200,               # peak merge parameter
    valid_chr='12,13',                   # held-out validation chromosomes (separated by comma)
    test_chr='7,8',                      # held-out test chromosomes (separated by comma)
    ignore_chr_y=True,                   # if true, don't include y chromosome in dataset
    ignore_auxiliary_chr=True,           # if true, don't include non-standard chromosomes (i.e. contigs)
    uncertain_N=True,                    # if true, set N to 0.25 for each nucleotide
    alphabet='ACGT',                     # alphabet order
    blacklist_path=blacklist_path,       # path to blacklist regions (i.e. unmappable)
    standard_coord=False,                # saves coordinates as chr1:start-end, otherwise chr_start_end
) 
