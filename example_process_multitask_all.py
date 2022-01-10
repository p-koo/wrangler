import os
import pandas as pd
import numpy as np
import utils
import generate_multitask


base_path = '../tf'
data_path = '../tf/bed'
metadata_path = os.path.join(data_path, 'metadata.tsv')
exp_bed_list_path = os.path.join(base_path, "exp_bed_list.tsv")
criteria=None
label_list=['Experiment target', 'Biosample term name', 'Experiment accession']


# generate an experiment list file w/ paths of bed files to be merged
summary = utils.generate_exp_bed_list(data_path, metadata_path, exp_bed_list_path, 
											 criteria, label_list)


# generate the dataset and save to h5 file
generate_multitask.process_data(
    exp_bed_list_path=exp_bed_list_path,
    prefix_save_path=os.path.join(base_path, 'all'),
    genome_path='../genomes/hg38.fa',
    chrom_sizes_path ='../genomes/hg38.chrom.sizes',
    seq_len=1000,
    merge_overlap_len=200,
    valid_chr='12,13',
    test_chr='7,8',
    ignore_chr_y=False,
    ignore_auxiliary_chr=False,
    uncertain_N=True,
    alphabet='ACGT',
) 
