"""
This script will generate a multitask binary classification dataset by processing
bed files from a tsv file (2 columns: name and path to bed files). The dataset is
saved to an h5 file.

Parameters:
-----------
exp_bed_list: str
    path to tsv file that has 2 columns of names and bed paths to be merged
prefix_save_path: str
    path + prefix for all output files
seq_len: str
    desired length of sequences
merge_overlap: str
    if two peak regions overlaps more than this amount, they will be re-centered and merged into a single sample.
genome_path: str
    Path to reference genome 
chrom_size:str
    Location for .chome.size file
valid_chr: str
    Held out chromosomes for validation set. Separate chromosomes with commas. Eg. 7,8,9
test_chr: str
    Held out chromosomes for validation set. Separate chromosomes with commas. Eg. 7,8,9

Example:
$ python generate_singletask.py -p ./pos.bed -n ./neg.bed -g ./genomes/hg38.fa -o ./test -s 200 -m 350 -v 12,13 -t 7,8

Test:
$ python generate_singletask.py -p ../data/test_singletask/ENCFF262MRD.bed.gz \
                                -n ../data/test_singletask/ENCFF759OLD.bed.gz \
                                -g ../data/genomes/hg38.fa \
                                -o ../data/test_singletask/REST_GM12878 \
                                -s 200 -m 350 -t 7,8 -v 12,13 \
                                --gc_match --uncertain_N --neg_pos_ratio 1.5 \
                                --seed 2019

"""

from optparse import OptionParser
import singletask




def main():
    parser = OptionParser()
    parser.add_option(
        '-p', "--pos_path",
        dest="pos_path",
        help="Path to bed file with positive peaks.",
    )
    parser.add_option(
        '-n', "--neg_path",
        dest="neg_path",
        help="Path to bed file with negative peaks.",
    )
    parser.add_option(
        '-g', "--genome_path",
        dest="genome_path",
        default="./hg38.fa",
        help="Path to reference genome -- needed to extract sequences from bed files.",
    )
    parser.add_option(
        '-o', "--prefix_save_path",
        dest="prefix_save_path",
        default="./experiment",
        help="Path and predix name for intermediate files (eg. merged peak bed file) and output h5 file.",
    )
    parser.add_option(
        '-s', "--seq_len",
        dest="seq_len",
        default=200,
        type='int',
        help="Length of selected sequence regions.",
    )
    parser.add_option(
        '-m', "--max_len_thresh",
        dest="max_len_thresh",
        default=350,
        type='int',
        help="Length of selected sequence regions.",
    )
    parser.add_option(
        '-f', '--filter_N', 
        dest = 'filter_N', 
        default=False, 
        action="store_true",
        help='Remove sequences with N charcters if set to True.'
    )
    parser.add_option(
        '-u', '--uncertain_N', 
        dest = 'uncertain_N', 
        default=False, 
        action="store_true",
        help='If True, sets positions in sequences with N charcters to 0.25. Otherwise, it remains a 0.'
    )
    parser.add_option(
        '-c', '--gc_match', 
        dest = 'gc_match', 
        default=False, 
        action="store_true",
        help='Match GC content of negative sequences with positive sequences.'
    )
    parser.add_option(
        '-r', "--neg_pos_ratio",
        dest="neg_pos_ratio",
        default=1.,
        type='float',
        help="Ratio of negative data to postivie data. If more negative data, will downsample to this ratio.",
    )
    parser.add_option(
        '-v', "--valid_chr", 
        dest="valid_chr", 
        default='6,7', 
        type="str", 
        help="Held out chromosomes for validation set. Separate chromosomes with commas."
    )
    parser.add_option(
        '-t', '--test_chr', 
        dest = 'test_chr', 
        default='5,8', 
        type='str', 
        help='Held out chromosomes for test set. Separate chromosomes with commas.'
    )
    parser.add_option(
        "--valid_frac",
        dest="valid_frac",
        default=0.,
        type='float',
        help="Valid fraction for random split -- only executed if test_frac not 0. Supercedes valid chrom.",
    )
    parser.add_option(
        "--test_frac",
        dest="test_frac",
        default=0.,
        type='float',
        help="Test fraction for random split -- only executed if not 0. Supercedes test chrom.",
    )
    parser.add_option(
        '-a', '--alphabet', 
        dest = 'alphabet', 
        default='ACGT', 
        type='str', 
        help='Alphabet for one-hot encoding.'
    )
    parser.add_option(
        '-z', "--seed",
        dest="seed",
        default=None,
        type='int',
        help="Random number seed for reproducibility.",
    )
    (options, args) = parser.parse_args()
    prefix_save_path = options.prefix_save_path

    singletask.process_data(
        pos_path=options.pos_path,
        neg_path=options.neg_path,
        prefix_save_path=options.prefix_save_path,
        genome_path=options.genome_path,
        seq_len=options.seq_len,
        max_len_thresh=options.max_len_thresh,
        filter_N=options.filter_N,
        uncertain_N=options.uncertain_N,
        neg_pos_ratio=options.neg_pos_ratio,
        gc_match=options.gc_match,
        valid_chr=options.valid_chr,
        test_chr=options.test_chr,
        valid_frac=options.valid_frac,
        test_frac=options.test_frac,
        alphabet=options.alphabet,
        seed=options.seed,
    )



#------------------------------------------------------------
# main function
#------------------------------------------------------------


    


#------------------------------------------------------------
# bedtools functions
#------------------------------------------------------------





if __name__ == "__main__":

    main()

