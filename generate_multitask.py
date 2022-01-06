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
from collections import OrderedDict
import os, sys, gzip, re, h5py
import numpy as np
import subprocess


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
    (options, args) = parser.parse_args()
    prefix_save_path = options.prefix_save_path


    process_data(
        exp_bed_list_path=options.exp_bed_list,
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
    )


#------------------------------------------------------------
# main function
#------------------------------------------------------------



def process_data(
    exp_bed_list_path,
    prefix_save_path,
    genome_path,
    chrom_sizes_path,
    seq_len=1000,
    merge_overlap_len=200,
    valid_chr='5,7',
    test_chr='3,6',
    ignore_chr_y=False,
    ignore_auxiliary_chr=False,
    uncertain_N=True,
    alphabet='ACGT',
):

    # generate a merged multitask bed file and activity table (outputs: prefix_merged.bed and prefix_activity.tsv)
    print("Merging bed files from: %s"%(exp_bed_list_path))
    generate_multitask_bed(
        exp_bed_list_path=exp_bed_list_path,
        prefix_save_path=prefix_save_path,
        seq_len=seq_len,
        merge_overlap_len=merge_overlap_len,
        chrom_sizes_path=chrom_sizes_path,
        bed_path=prefix_save_path+'_merged.bed',
        activity_path=prefix_save_path+'_activity.tsv',
        ignore_chr_y=ignore_chr_y,
        ignore_auxiliary_chr=ignore_auxiliary_chr,
    )
    print("Saved to: %s_merged_bed and %s_activity.tsv"%(prefix_save_path, prefix_save_path))


    # convert merged bed file into fasta file (eg. prefix.fa)
    print("Extracting sequences from reference genome: %s"%(genome_path))
    subprocess.call(
        "bedtools getfasta -fi {} -s -bed {} -fo {}".format(
            genome_path, prefix_save_path+"_merged.bed", prefix_save_path+".fa"
        ),
        shell=True,
    )

    # load sequences from fasta file
    print("Processing files -- converting to one-hot and generating labels")
    seq_vecs = convert_fasta_to_onehot(prefix_save_path+".fa", alphabet, uncertain_N)
    os.remove(prefix_save_path+".fa")

    # load scores
    seq_targets, target_labels = parse_targets(prefix_save_path+"_activity.tsv")

    # align and construct input matrix
    seqs, targets, coords = align_seqs_targets(seq_vecs, seq_targets, sort=False)

    # split dataset by chromosome
    print("Splittig dataset into:")
    print("    validation: %s"%(valid_chr))
    print("    test: %s"%(test_chr))
    train, valid, test = split_dataset_chrom(seqs, targets, coords, valid_chr, test_chr)

    # data splits
    prefix_save_path = prefix_save_path+'.h5'
    print("Saving to: %s"%(prefix_save_path))
    with h5py.File(prefix_save_path, "w") as fout:
        fout.create_dataset("x_train", data=train[0], compression="gzip")
        fout.create_dataset("y_train", data=train[1], compression="gzip")
        fout.create_dataset("coords_train", data=train[2].astype("S"), compression="gzip")
        fout.create_dataset("x_valid", data=valid[0], compression="gzip")
        fout.create_dataset("y_valid", data=valid[1], compression="gzip")
        fout.create_dataset("coords_valid", data=valid[2].astype("S"), compression="gzip")
        fout.create_dataset("x_test", data=test[0], compression="gzip")
        fout.create_dataset("y_test", data=test[1], compression="gzip")
        fout.create_dataset("coords_test", data=test[2].astype("S"), compression="gzip")
    print('%d train sequences'%(len(train[0])))
    print('%d valid sequences'%(len(valid[0])))
    print('%d test sequences '%(len(test[0])))




#------------------------------------------------------------
# Useful functions
#------------------------------------------------------------



def generate_multitask_bed(
    exp_bed_list_path,
    prefix_save_path=None,
    seq_len=1000,
    merge_overlap_len=200,
    chrom_sizes_path=None,
    bed_path="./merged.bed",
    activity_path="./activity.tsv",
    ignore_chr_y=False,
    ignore_auxiliary_chr=False,
):

    """Merge multiple bed files to select sample sequence regions with at least one
    peak.
    This function outputs a .bed file in the specified directory containing seven
    columns: chromosome, sequence start, sequence end, name, score, strand, and indexes
    of experiments that have a peak in this region.
    Parameters
    ----------
    exp_bed_list_path: str
        Location of the sample file containing experiment label and their corresponding
        file locations. Should be a two column text file, first row contains label,
        second row contain directory for the .bed/.bed.gz file.
    seq_len: int, optional
        Length of the sequence region per sample in output. Default to 1000.
    merge_overlap_len: int, optional
        After adjusting peaks into seq_len, if two peak regions overlaps more than
        this amount, they will be re-centered and merged into a single sample. Defaults
        to 200.
    prefix_save_path: str, optional
        Location and naming of the output bed file. Default to 'features.bed'
    chr_lenghts_file: str, optional
        Location of the chrom.sizes file. Default to None.
    db_act_file: str, optional
        Location of the existing database activity table. Defaults to None.
    db_bed: str, optional
        Location of the existing database .bed file. Defaults to None.
    ignore_auxiliary_chr: bool, optional
        Ignore auxiliary chromosomes. Defaults to False.
    no_db_acticity: bool, optional
        Whether to pass along the activities of the database sequences. Defaults to
        False.
    ignor_y: bool, optional
        Ignore Y chromsosome features. Defaults to False.
    Returns
    -------
    None
    Examples
    --------
    >>> multitask_bed_generation(
        example_file,chr_lens_file='/data/hg38.chrom.size',
        seq_len=1000,merge_overlap_len=200,prefix_save_path='/data/bed')
    """

    if not exp_bed_list_path:
        raise Exception(
            "Must provide file labeling the targets and providing BED file paths."
        )

    # read in targets and assign them indexes into the db
    target_names = []
    target_bedpaths = []
    target_index = []
    for line in open(exp_bed_list_path):
        a = line.split()
        if len(a) != 2:
            print(a)
            print(
                "Each row of the target BEDS file must contain a label and BED file"
                " separated by whitespace",
                file=sys.stderr,
            )
            sys.exit(1)
        target_index.append(len(target_names))
        target_names.append(a[0])
        target_bedpaths.append(a[1])

    # read in chromosome lengths
    chr_lens = {}
    if chrom_sizes_path is not None:
        chr_lens = {}
        for line in open(chrom_sizes_path):
            a = line.split()
            chr_lens[a[0]] = int(a[1])
    else:
        print(
            "Warning: chromosome lengths not provided, so regions near ends may be"
            " incorrect.",
            file=sys.stderr,
        )

    #################################################################
    # save peaks from all bed files to unified chromosome-specific files
    #################################################################
    chr_file_dict = {}
    chr_fout_dict = {}
    for i, bedpath in enumerate(target_bedpaths):
        if bedpath[-3:] == ".gz":
            bed_in = gzip.open(bedpath, "rt")
        else:
            bed_in = open(bedpath)

        for line in bed_in:
            if not line.startswith("#"):
                a = line.split("\t")
                a[-1] = a[-1].rstrip()

                # hash by chrom/strand
                chrom = a[0]

                # reset start and end according to midpoint
                start = int(a[1])
                end = int(a[2])
                mid = int((start + end)/2)
                a[1] = str(mid)
                a[2] = str(mid + 1)

                # get strand info (if any)
                strand = "+"
                if len(a) > 5 and a[5] in "+-":
                    strand = a[5]

                # specify the target index
                while len(a) < 7:
                    a.append("")
                a[5] = strand
                a[6] = str(target_index[i])
    
                # open chromosome file if doesn't already exist
                chr_key = (chrom, strand)
                if chr_key not in chr_fout_dict:
                    chr_file_dict[chr_key] = "%s_%s_%s.bed" % (
                        prefix_save_path,
                        chrom,
                        strand,
                    )
                    chr_fout_dict[chr_key] = open(chr_file_dict[chr_key], "w")

                # append entry into chromosome specific file, merging all bed files    
                print("\t".join(a[:7]), file=chr_fout_dict[chr_key])

        bed_in.close()

    # close chromosome-specific files
    for chr_key in chr_fout_dict:
        chr_fout_dict[chr_key].close()

    # if ignore Y, filter data
    if ignore_chr_y:
        for orient in "+-":
            chr_key = ("chrY", orient)
            if chr_key in chr_file_dict:
                print("Ignoring chrY %s" % orient, file=sys.stderr)
                os.remove(chr_file_dict[chr_key])
                del chr_file_dict[chr_key]

    # if ignore auxiliary, filter data
    if ignore_auxiliary_chr:
        primary_re = re.compile("chr\\d+$")
        for chr_key in chr_file_dict.keys():
            chrom, strand = chr_key
            primary_m = primary_re.match(chrom)
            if not primary_m and chrom != "chrX":
                print("Ignoring %s %s" % (chrom, strand), file=sys.stderr)
                os.remove(chr_file_dict[chr_key])
                del chr_file_dict[chr_key]


    #################################################################
    # sort chromosome-specific files
    #################################################################
    for chr_key in chr_file_dict:
        chrom, strand = chr_key
        chrom_sbed = "%s_%s_%s_sort.bed" % (prefix_save_path, chrom, strand)
        sort_cmd = "sortBed -i %s > %s" % (chr_file_dict[chr_key], chrom_sbed)
        subprocess.call(sort_cmd, shell=True)
        os.remove(chr_file_dict[chr_key])
        chr_file_dict[chr_key] = chrom_sbed


    #################################################################
    # parse chromosome-specific files
    #################################################################
    final_bed_out = open("%s_merged.bed" % prefix_save_path, "w")

    for chr_key in chr_file_dict:
        chrom, strand = chr_key

        open_peaks = []
        for line in open(chr_file_dict[chr_key], "rt"):
            a = line.split("\t")
            a[-1] = a[-1].rstrip()

            # construct Peak
            peak_start = int(a[1])
            peak_end = int(a[2])
            peak_act = _activity_set(a[6])
            peak = Peak(peak_start, peak_end, peak_act)
            peak.extend(seq_len, chr_lens[chrom])

            # check if peak is in new region or already active region
            if len(open_peaks) == 0:
                # initialize open peak
                open_end = peak.end
                open_peaks = [peak]

            else:
                # if beyond existing open peak
                if open_end - merge_overlap_len <= peak.start:
                    # close open peak
                    mpeaks = _merge_peaks(
                        open_peaks,
                        seq_len,
                        merge_overlap_len,
                        chr_lens[chrom],
                    )

                    # print to file
                    for mpeak in mpeaks:
                        print(mpeak.bed_str(chrom, strand), file=final_bed_out)

                    # initialize open peak
                    open_end = peak.end
                    open_peaks = [peak]

                else:
                    # extend open peak
                    open_peaks.append(peak)
                    open_end = max(open_end, peak.end)

        if len(open_peaks) > 0:
            # close open peak
            mpeaks = _merge_peaks(
                open_peaks, seq_len, merge_overlap_len, chr_lens[chrom]
            )

            # print to file
            for mpeak in mpeaks:
                print(mpeak.bed_str(chrom, strand), file=final_bed_out)
    final_bed_out.close()

    # clean
    for chr_key in chr_file_dict:
        os.remove(chr_file_dict[chr_key])


    #################################################################
    # generate activity table (i.e. names and binary labels for each class)
    #################################################################
    final_act_out = open("%s_activity.tsv" % prefix_save_path, "w")

    # print header
    cols = [""] + target_names
    print("\t".join(cols), file=final_act_out)

    # get coordinates and labels from bed file and save as activity table
    for line in open("%s_merged.bed" % prefix_save_path):
        a = line.rstrip().split("\t")
        peak_id = "%s:%s-%s(%s)" % (a[0], a[1], a[2], a[5])

        # construct full activity vector
        peak_act = [0] * len(target_names)
        for ai in a[6].split(","):
            if ai != ".":
                peak_act[int(ai)] = 1

        # print line
        cols = [peak_id] + peak_act
        print("\t".join([str(c) for c in cols]), file=final_act_out)
    final_act_out.close()


class Peak:
    """Peak representation
    Attributes:
        start (int) : peak start
        end   (int) : peak end
        act   (set[int]) : set of target indexes where this peak is active.
    """

    def __init__(self, start, end, act):
        self.start = start
        self.end = end
        self.act = act

    def extend(self, ext_len, chr_len):
        """Extend the peak to the given length
        Args:
            ext_len (int) : length to extend the peak to
            chr_len (int) : chromosome length to cap the peak at
        """
        mid = _find_midpoint(self.start, self.end)
        self.start = max(0, mid - ext_len / 2)
        self.end = self.start + ext_len
        if chr_len and self.end > chr_len:
            self.end = chr_len
            self.start = self.end - ext_len

    def bed_str(self, chrom, strand):
        """Return a BED-style line
        Args:
            chrom (str)
            strand (str)
        """
        if len(self.act) == 0:
            act_str = "."
        else:
            act_str = ",".join([str(ai) for ai in sorted(list(self.act))])
        cols = (
            chrom,
            str(int(self.start)),
            str(int(self.end)),
            ".",
            "1",
            strand,
            act_str,
        )
        return "\t".join(cols)

    def merge(self, peak2, ext_len, chr_len):
        """Merge the given peak2 into this peak
        Args:
            peak2 (Peak)
            ext_len (int) : length to extend the merged peak to
            chr_len (int) : chromosome length to cap the peak at
        """
        # find peak midpoints
        peak_mids = [_find_midpoint(self.start, self.end)]
        peak_mids.append(_find_midpoint(peak2.start, peak2.end))

        # weight peaks
        peak_weights = [1 + len(self.act)]
        peak_weights.append(1 + len(peak2.act))

        # compute a weighted average
        merge_mid = int(0.5 + np.average(peak_mids, weights=peak_weights))

        # extend to the full size
        merge_start = max(0, merge_mid - ext_len / 2)
        merge_end = merge_start + ext_len
        if chr_len and merge_end > chr_len:
            merge_end = chr_len
            merge_start = merge_end - ext_len

        # merge activities
        merge_act = self.act | peak2.act

        # set merge to this peak
        self.start = merge_start
        self.end = merge_end
        self.act = merge_act


def _find_midpoint(start, end):
    """ Find the midpoint coordinate between start and end """
    mid = (start + end) / 2
    return int(mid)


def _merge_peaks(peaks, peak_size, merge_overlap_len, chr_len):
    """Merge and the list of Peaks.
    Repeatedly find the closest adjacent peaks and consider
    merging them together, until there are no more peaks
    we want to merge.
    Attributes:
        peaks (list[Peak]) : list of Peaks
        peak_size (int) : desired peak extension size
        chr_len (int) : chromsome length
    Returns:
        Peak representing the merger
    """
    max_overlap = merge_overlap_len
    while len(peaks) > 1 and max_overlap >= merge_overlap_len:
        # find largest overlap
        max_i = 0
        max_overlap = peaks[0].end - peaks[1].start
        for i in range(1, len(peaks) - 1):
            peaks_overlap = peaks[i].end - peaks[i + 1].start
            if peaks_overlap > max_overlap:
                max_i = i
                max_overlap = peaks_overlap

        if max_overlap >= merge_overlap_len:
            # merge peaks
            peaks[max_i].merge(peaks[max_i + 1], peak_size, chr_len)

            # remove merged peak
            peaks = peaks[: max_i + 1] + peaks[max_i + 2 :]

    return peaks


def _activity_set(act_cs):
    """Return a set of ints from a comma-separated list of int strings.
    Attributes:
        act_cs (str) : comma-separated list of int strings
    Returns:
        set (int) : int's in the original string
    """
    ai_strs = [ai for ai in act_cs.split(",")]

    if ai_strs[-1] == "":
        ai_strs = ai_strs[:-1]

    if ai_strs[0] == ".":
        aset = set()
    else:
        aset = set([int(ai) for ai in ai_strs])

    return aset


def convert_fasta_to_onehot(fasta_file, alphabet='ACGT', uncertain_N=True):
    seq_vecs = OrderedDict()
    for line in open(fasta_file):
        if line[0] == ">":
            name = line[1:].rstrip()
        else:
            seq_vecs[name] = convert_one_hot([line.rstrip().upper()], alphabet, uncertain_N)

    return seq_vecs




def convert_one_hot(sequences, alphabet="ACGT", uncertain_N=True):
    """Convert flat array of sequences to one-hot representation.
    **Important**: all letters in `sequences` *must* be contained in `alphabet`, and
    all sequences must have the same length.
    Parameters
    ----------
    sequences : numpy.ndarray of strings
        The array of strings. Should be one-dimensional.
    alphabet : str
        The alphabet of the sequences.
    Returns
    -------
    Numpy array of sequences in one-hot representation. The shape of this array is
    `(len(sequences), len(sequences[0]), len(alphabet))`.
    Examples
    --------
    >>> one_hot(["TGCA"], alphabet="ACGT")
    array([[[0., 0., 0., 1.],
            [0., 0., 1., 0.],
            [0., 1., 0., 0.],
            [1., 0., 0., 0.]]])
    """
    A = len(alphabet)
    alphabet += 'N' # to capture non-sense characters

    sequences = np.asanyarray(sequences)
    if sequences.ndim != 1:
        raise ValueError("array of sequences must be one-dimensional.")
    n_sequences = sequences.shape[0]
    seq_len = len(sequences[0])

    # Unpack strings into 2D array, where each point has one character.
    s = np.zeros((n_sequences, seq_len), dtype="U1")
    for i in range(n_sequences):
        s[i] = list(sequences[i])

    # Make an integer array from the string array.
    pre_onehot = np.zeros(s.shape, dtype=np.uint8)
    for i, letter in enumerate(alphabet):
        # do nothing on 0 because array is initialized with zeros.
        if i:
            pre_onehot[s == letter] = i

    # create one-hot representation
    n_classes = len(alphabet)
    one_hot = np.eye(n_classes)[pre_onehot]

    # remove nonsense character
    one_hot = one_hot[:,:,:A]

    # replace positions with N with 0.25 if true
    if uncertain_N:
        for n,x in enumerate(one_hot):
            index = np.where(np.sum(x, axis=-1) == 0)[0]
            one_hot[n,index,:] = 0.25
    return one_hot



def parse_targets(targets_file):
    seq_targets = OrderedDict()
    for line in open(targets_file):
        a = line.split()
        try:
            seq_targets[a[0]] = np.array([int(a[i]) for i in range(1, len(a))])
        except Exception:
            target_labels = a
    return seq_targets, target_labels


def align_seqs_targets(seq_vecs, seq_targets, sort=True):
    if sort:
        seq_headers = sorted(seq_vecs.keys())
    else:
        seq_headers = seq_vecs.keys()

    # construct lists of vectors
    targets = []
    seqs = []
    coords = []
    for header in seq_headers:
        coords.append(header)
        seqs.append(seq_vecs[header])
        targets.append(seq_targets[header])
    return np.vstack(seqs), np.vstack(targets), np.array(coords)
    



def split_dataset_chrom(x, y, coords, valid_chr, test_chr):

    def parse_held_out_chrom(x, y, coords, split_chr):
        """extract data according to chromosome"""

        # get list of chromosomes for each data entry
        chrom_list = np.array([n.split(':')[0] for n in coords])

        # parse held out chrom data
        index = []
        for chr_index in split_chr:
            index.append(np.where(chr_index == chrom_list)[0])
        index = np.concatenate(index)
        return x[index], y[index], coords[index]

    def remove_held_out_chrom(x, y, coords, remove_chr):
        """remove data according to chromosome"""

        # get list of chromosomes for each data entry
        chrom_list = np.array([n.split(':')[0] for n in coords])

        # filter out held out chroms in data
        chrom_filt = np.copy(chrom_list)
        x_filt = np.copy(x)
        y_filt = np.copy(y)
        coords_filt = np.copy(coords)
        for chr_index in remove_chr:
            index = np.where(chr_index == chrom_filt)[0]
            chrom_filt = np.delete(chrom_filt, index, axis=0)
            x_filt = np.delete(x_filt, index, axis=0)
            y_filt = np.delete(y_filt, index, axis=0)
            coords_filt = np.delete(coords_filt, index, axis=0)
        return x_filt, y_filt, coords_filt

    # get list of held out chromosomes 
    valid_chr = ['chr'+chr_index for chr_index in valid_chr.split(',')]
    test_chr  = ['chr'+chr_index for chr_index in test_chr.split(',')]

    # generate dataset
    test  = parse_held_out_chrom(x, y, coords, test_chr)
    valid = parse_held_out_chrom(x, y, coords, valid_chr)
    train = remove_held_out_chrom(x, y, coords, valid_chr+test_chr)

    return train, valid, test


if __name__ == "__main__":

    main()

