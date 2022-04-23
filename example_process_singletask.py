import singletask


tf_name = 'REST'  
prefix_save_path = os.path.join('../data/test_singletask', tf_name)     # path to save file (eg. prefix_save_path + '.h5')
pos_path = '../data/test_singletask/ENCFF262MRD.bed.gz'                 # path to positive bd file
neg_path = '../data/test_singletask/ENCFF759OLD.bed.gz'                 # path to negative bed file
genome_path = '../data/genomes/hg38.fa'                                 # path to reference genome
blacklist_path = '../data/hg38_unmap.fa'                                # path to blacklist (i.e. unmappable) regions
seq_len = 200               # sequence length for the dataset
max_len_thresh = 350        # remove bed entries with size greater than this threshold
filter_N = False            # remove sequences with N characters
uncertain_N = True          # use 0.25 for N characters
gc_match = True             # sub-sample negative set to match GC content of positive set
neg_pos_ratio = 1.350       # ratio of the amount of negative samples to positive samples (Eg. ratio=1 means balanced dataset)
valid_chr = '12,13'         # validation set from held out chromosome
test_chr = '7,8'            # test set from held out chromosome 
alphabet = 'ACGT'           # alphabet for one-hot encoding -- note N is treated as null vector unless uncertain_N flag is true
seed = 9919                 # random number seed for reproducibility -- random sampling occurs when downsampling negative set


singletask.process_data(
    pos_path=pos_path,
    neg_path=neg_path,
    prefix_save_path=prefix_save_path,
    genome_path=genome_path,
    seq_len=seq_len,
    max_len_thresh=max_len_thresh,
    filter_N=filter_N,
    uncertain_N=uncertain_N,
    gc_match=gc_match,
    neg_pos_ratio=neg_pos_ratio,
    valid_chr=valid_chr,
    test_chr=test_chr,
    alphabet=alphabet,
    blacklist_path=blacklist_path,
    standard_coord=False,
    seed=seed,
)
