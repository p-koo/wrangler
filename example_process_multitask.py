import utils
import multitask


cell_type = 'HepG2'
base_path = '../tf'                                                 
data_path = utils.make_directory(base_path, cell_type)                           # path to data
genome_path = '../genomes/hg38.fa'                                               # path to reference genome
chrom_sizes_path = '../genomes/hg38.chrom.sizes'                                 # path to chrom sizes
metadata_path = os.path.join(base_path, 'metadata.tsv')                          # path to ENCODE metadata table
exp_list_bed_path = os.path.join(base_path, "exp_list_bed.tsv")                  # filename to save list of experiments-2columns (name    datapath)
criteria={'File assembly': "GRCh38",                                             # criteria for narrowing down files to process
          'File format': 'bed narrowPeak',
          'Output type': 'IDR thresholded peaks',
          'Biosample term name': cell_type,             
          }
label_list=['Experiment target', 'Biosample term name', 'Experiment accession']  # used to name experiment


# generate an experiment list file w/ paths of bed files to be merged
summary = utils.generate_exp_list_from_metadata(data_path, metadata_path, exp_list_bed_path, 
											    criteria, label_list, download=False, ext='bed')


# generate the dataset and save to h5 file
multitask.process_data(
    exp_bed_list_path=exp_list_bed_path,
    prefix_save_path=data_path,
    genome_path=genome_path,
    chrom_sizes_path =chrom_sizes_path,
    seq_len=1000,
    merge_overlap_len=200,
    valid_chr='12,13',
    test_chr='7,8',
    ignore_chr_y=False,
    ignore_auxiliary_chr=False,
    uncertain_N=True,
    alphabet='ACGT',
    blacklist_path=None,
) 
