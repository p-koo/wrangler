import os
import pandas as pd
import numpy as np
import utils
import generate_singletask

#--------------------------------------------------------------------------------------


cell_type = 'A549'
base_path = '../data/tf/binary/'
data_path = '../data/tf/binary/all'
save_path = utils.make_directory(base_path, cell_type)


#--------------------------------------------------------------------------------------
# Parse TFs for a given cell type
#--------------------------------------------------------------------------------------

metadata_path = os.path.join(data_path, 'metadata.tsv')
exp_bed_list_path = os.path.join(save_path, "exp_bed_list.tsv")
criteria={'File assembly': "GRCh38",
          'File format': 'bed narrowPeak',
          'Biological replicate(s)': '1, 2',
          'Output type': 'IDR thresholded peaks',
          'Biosample term name': cell_type,
          }
label_list=['Experiment target', 'Biosample term name', 'Experiment accession']


# process transcription factor bed files for a given cell type
summary = utils.generate_exp_bed_list(data_path, metadata_path, exp_bed_list_path, 
									  criteria, label_list)

#--------------------------------------------------------------------------------------
# Parse accessibility data for the same cell type
#--------------------------------------------------------------------------------------

print('Searching for matching DNase-seq files')

# find a matching DNase/ATAC dataset
accessible_data_path = '../data/accessibility/mixed_cell_line'
neg_criteria={'Biosample term name': cell_type,
              'Assay': 'DNase-seq',
              'File assembly': "GRCh38",
              'File format': 'bed narrowPeak',
              'File analysis title': 'ENCODE4 v3.0.0-alpha.2 GRCh38', 
              }
filenames, filepaths = utils.get_path_from_metadata(accessible_data_path, neg_criteria)

if filenames is not None:
    neg_name = filenames[0] # take first file
    neg_path = filepaths[0]
    print("Using DNase-seq negative dataset: "+neg_path)
else:
    print('Searching for matching ATAC-seq files')
    neg_criteria['Assay'] = 'ATAC-seq'
    filenames, filepaths = utils.get_path_from_metadata(accessible_data_path, neg_criteria)
    if filenames is not None:
        neg_name = filenames[0] # take first files
        neg_path = filepaths[0]
        print("Using ATAC-seq negative dataset: "+neg_path)
    
#--------------------------------------------------------------------------------------
# Generate a binary classification dataset for each TF
#--------------------------------------------------------------------------------------

genome_path = '../data/genomes/hg38.fa'
seq_len = 200               # sequence length for the dataset
max_len_thresh = 350        # remove bed entries with size greater than this threshold
filter_N = False            # remove sequences with N characters
uncertain_N = True          # use 0.25 for N characters
gc_match = True             # sub-sample negative set to match GC content of positive set
neg_pos_ratio = 1.50       # ratio of the amount of negative samples to positive samples (Eg. ratio=1 means balanced dataset)
valid_chr = '12,13'         # validation set from held out chromosome
test_chr = '7,8'            # test set from held out chromosome 
alphabet = 'ACGT'           # alphabet for one-hot encoding -- note N is treated as null vector unless uncertain_N flag is true
seed = 9919                 # random number seed for reproducibility -- random sampling occurs when downsampling negative set

if filenames is not None:
    for name, pos_path in summary:
        generate_singletask.process_data(
            pos_path=pos_path,
            neg_path=neg_path,
            prefix_save_path=os.path.join(save_path, name+'_'+neg_name),
            genome_path=genome_path,
            seq_len=seq_len,
            max_len_thresh=max_len_thresh,
            filter_N=filter_N,
            uncertain_N=uncertain_N,
            gc_match=gc_match,
            neg_pos_ratio=neg_pos_ratio,
            valid_chr=valid_chr,
            test_chr=test_chr,
            valid_frac=None, 
            test_frac=None,    # if set to a number, randomly split data into train, valid, and test
            alphabet=alphabet,
            seed=seed,
        )

else:
    print("Error: Could not find a suitable negative dataset.")