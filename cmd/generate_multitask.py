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
merge_overlap_len: str
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
$ python generate_multitask.py -i ./exp_bed_list.tsv -o ./test  -g ./genomes/hg38.fa -c ./genomes/hg38.chrom.sizes -s 1000 -m 200 -v 12,13 -t 7,8


test:
$ python generate_multitask.py -i ../data/test_multitask/test.txt -o ../data/test_multitask/test  -g ../data/genomes/hg38.fa -c ../data/genomes/hg38.chrom.sizes -s 1000 -m 200 -v 12,13 -t 7,8

"""

from optparse import OptionParser
import multitask

def main():
    parser = OptionParser()
    parser.add_option(
        '-i', "--exp_bed_list",
        dest="exp_bed_list",
        default="./exp_bed_list.tsv",
        help="Path to file that contains 2-column list: (1) experiment names and (2) bed paths to be merge.",
    )
    parser.add_option(
        '-o', "--prefix_save_path",
        dest="prefix_save_path",
        default="./merged_data",
        help="Path and predix name for intermediate files (eg. merged peak bed file) and output h5 file.",
    )
    parser.add_option(
        '-s', "--seq_len",
        dest="seq_len",
        default=1000,
	    type='int',
        help="Length of selected sequence regions.",
    )
    parser.add_option(
        '-m', "--merge_overlap_len",
        dest="merge_overlap_len",
        default=200,
        type='int',
        help="If two peak regions overlaps more than this amount, they will be"
             " re-centered and merged into a single sample."
    )
    parser.add_option(
        '-g', "--genome_path",
        dest="genome_path",
        default="./hg38.fa",
        help="Path to reference genome -- needed to extract sequences from bed files.",
    )
    parser.add_option(
        '-c', "--chrom_size", 
        dest="chrom_size", 
        default="./hg38.chrom.sizes", 
        help="Path to chromosome size file."
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
        '-y', '--ignore_chr_y', 
        dest = 'ignore_chr_y', 
        default=False, 
        action="store_true",
        help='Ignore Y chromosome.'
    )
    parser.add_option(
        '-x', '--ignore_auxiliary_chr', 
        dest = 'ignore_auxiliary_chr', 
        default=False, 
        action="store_true",
        help='Ignore auxiliary chromosome.'
    )
    parser.add_option(
        '-a', '--alphabet', 
        dest = 'alphabet', 
        default='ACGT', 
        type='str', 
        help='Alphabet for one-hot encoding. Default: ACGT'
    )
    parser.add_option(
        '-u', '--uncertain_N', 
        dest = 'uncertain_N', 
        default=False, 
        action="store_true",
        help='If True, sets positions in sequences with N charcters to 0.25. Otherwise, it remains a 0.'
    )
    parser.add_option(
        '-b', "--blacklist", 
        dest="blacklist_path", 
        default="./hg38_umap.bed", 
        help="Path to blacklisted bed file."
    )
    parser.add_option
    (options, args) = parser.parse_args()
    prefix_save_path = options.prefix_save_path


    multitask.process_data(
        exp_list_path=options.exp_bed_list,
        prefix_save_path=options.prefix_save_path,
        genome_path=options.genome_path,
        chrom_sizes_path=options.chrom_size,
        seq_len=options.seq_len,
        merge_overlap_len=options.merge_overlap_len,
        valid_chr=options.valid_chr,
        test_chr=options.test_chr,
        ignore_chr_y=options.ignore_chr_y,
        ignore_auxiliary_chr=options.ignore_auxiliary_chr,
        uncertain_N=options.uncertain_N,
        alphabet=options.alphabet,
        blacklist_path=options.blacklist_path,
    )


#------------------------------------------------------------
# main function
#------------------------------------------------------------









if __name__ == "__main__":

    main()

