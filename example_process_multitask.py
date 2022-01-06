import os
import pandas as pd
import numpy as np
import encode_utils
import generate_multitask_binary


cell_type = 'HCT116'

base_path = './tf/binary/'
data_path = './tf/binary/all'
metadata_path = os.path.join(data_path, 'metadata.tsv')
exp_bed_list_path = os.path.join(base_path, "exp_bed_list.tsv")
criteria={'File assembly': "GRCh38",
          'File format': 'bed narrowPeak',
          'Output type': 'IDR thresholded peaks',
          'Biosample term name': cell_type,
          }
label_list=['Experiment target', 'Biosample term name', 'Experiment accession']

summary = encode_utils.generate_exp_bed_list(data_path, metadata_path, exp_bed_list_path, 
											 criteria, label_list)


generate_multitask_binary.process_data(
    exp_bed_list_path,
    prefix_path=os.path.join(base_path, cell_type),
    genome_path='./genomes/hg38.fa',
    chrom_sizes_path ='./genomes/hg38.chrom.sizes',
    feature_size=1000,
    merge_overlap=200,
    ignore_y=False,
    ignore_auxiliary=False,
    valid_chrom='8,9',
    test_chrom='6,7',
    alphabet='ACGT',
    uncertain_N=True,
) 
