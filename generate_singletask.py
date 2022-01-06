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
import os, h5py, gzip
import subprocess
import numpy as np
import pandas as pd




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

    process_data(
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


def process_data(
    pos_path,
    neg_path,
    prefix_save_path='./sample',
    genome_path='./hg38.fa',
    seq_len=200,
    max_len_thresh=300,
    filter_N=False,
    uncertain_N=True,
    gc_match=True,
    neg_pos_ratio=1,
    valid_chr='6,7',
    test_chr='4,8',
    valid_frac=None,
    test_frac=None,
    alphabet="ACGT",
    seed=None,
):
    """Preprocess data for a single-class task."""
    print(prefix_save_path)
    print("Processing positive label data")

    # set random number seed
    if seed is not None:
        np.random.seed(seed)

    # remove extremely large peaks
    tf_filtered_path = prefix_save_path + "_pos_filtered.bed"
    filter_max_length(pos_path, tf_filtered_path, max_len_thresh)

    # create new bed file with seq_len enforced
    pos_bed_path = prefix_save_path + "_pos_" + str(seq_len) + ".bed"
    enforce_constant_size(tf_filtered_path, pos_bed_path, seq_len)
    pos_coords = get_bed_coords(pos_bed_path)

    # extract sequences from bed file and save to fasta file
    pos_fasta_path = prefix_save_path + "_pos.fa"
    bedtools_getfasta(
        pos_bed_path, genome_path, output_path=pos_fasta_path, strand=True
    )

    # parse sequence from fasta file
    pos_seq,_ = parse_fasta(pos_fasta_path)

    if filter_N:
        N = len(pos_seq)
        # filter sequences with absent nucleotides
        pos_seq, good_index = filter_nonsense_sequences(pos_seq)
        pos_coords = pos_coords[good_index]
        print('    Filtered %d positive sequences with N'%(N-len(pos_seq)))

    # convert filtered sequences to one-hot representation
    pos_one_hot = convert_one_hot(pos_seq, alphabet, uncertain_N=uncertain_N)


    print("Processing negative label data")

    # get non-overlap between pos peaks and neg peaks
    neg_bed_path = prefix_save_path + "_nonoverlap.bed"
    bedtools_intersect(
        neg_path, pos_path, neg_bed_path, write_a=True, nonoverlap=True
    )

    # create new bed file with seq_len enforced
    neg_bed_path2 = prefix_save_path + "_neg_" + str(seq_len) + ".bed"
    enforce_constant_size(neg_bed_path, neg_bed_path2, seq_len)
    neg_coords = get_bed_coords(neg_bed_path2)

    # extract sequences from bed file and save to fasta file
    neg_fasta_path = prefix_save_path + "_neg.fa"
    bedtools_getfasta(
        neg_bed_path2, genome_path, output_path=neg_fasta_path, strand=True
    )

    # parse sequence and chromosome from fasta file
    neg_seq,_ = parse_fasta(neg_fasta_path)

    if filter_N:
        # filter sequences with absent nucleotides
        N = len(neg_seq)
        neg_seq, good_index = filter_nonsense_sequences(neg_seq)
        neg_coords = neg_coords[good_index]
        print('    Filtered %d negative sequences with N'%(N-len(neg_seq)))

    # convert filtered sequences to one-hot representation
    neg_one_hot = convert_one_hot(neg_seq, alphabet, uncertain_N=uncertain_N)

    # match GC content
    if gc_match:
        print("    Matching GC content between negative and positive data.")
        neg_one_hot, match_index = match_gc_content(pos_one_hot, neg_one_hot, neg_pos_ratio=neg_pos_ratio)
        neg_coords = neg_coords[match_index]

    else:
        num_pos = len(pos_one_hot)
        num_neg = len(neg_one_hot)
        if num_pos < num_neg:
            num_keep = int(num_pos * neg_pos_ratio)
            index = np.random.permutation(len(neg_one_hot))[:num_keep]
            neg_one_hot = neg_one_hot[index]
            neg_coords = neg_coords[index]
            print("    Downsampling negative set to %d"%(num_keep))

    print("%d positive sequences and %d negative sequences"%(len(pos_one_hot), len(neg_one_hot)))

    # merge positive and negative labels
    one_hot = np.vstack([pos_one_hot, neg_one_hot])
    labels = np.vstack(
        [np.ones((len(pos_one_hot), 1)), np.zeros((len(neg_one_hot), 1))]
    )
    coords = np.concatenate([pos_coords, neg_coords])

    # split data into training, validation, and test sets
    if test_frac:
        print("Splittig dataset into:")
        print("    random validation: %.2f"%(valid_frac))
        print("    random test: %.2f"%(test_frac))
        # shuffle indices for train, validation, and test sets
        train, valid, test = split_dataset_random(
            one_hot, labels, coords, valid_frac=valid_frac, test_frac=test_frac
        )
    else:
        print("Splittig dataset into:")
        print("    validation: %s"%(valid_chr))
        print("    test: %s"%(test_chr))
        train, valid, test = split_dataset_chrom(
            one_hot, labels, coords, valid_chr, test_chr
        )

    # save to hdf5 file
    file_path = prefix_save_path + ".h5"
    print("Saving to: %s"%(file_path))
    with h5py.File(file_path, "w") as fout:
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

    # clean up directory
    os.remove(tf_filtered_path)
    os.remove(pos_bed_path)
    os.remove(pos_fasta_path)

    os.remove(neg_bed_path)
    os.remove(neg_bed_path2)
    os.remove(neg_fasta_path)
    


#------------------------------------------------------------
# Useful functions
#------------------------------------------------------------



def filter_max_length(bed_path, output_path, max_len=1000):
    """Function to plot histogram of bed file peak sizes; automatically infers
    compression  from extension and allows for user-input in removing outlier sequence
    Parameters
    -----------
    bed_path : <str>
        Path to bed file.
    output_path : <int>
        Path to filtered bed file.
    max_len: <int>
        Cutoff for maximum length of peak -- anything above will be filtered out.
    Returns
    ----------
    None
    """

    # check if bedfile is compressed
    compression = "gzip" if _is_gzipped(bed_path) else None

    # load bed file
    df = pd.read_table(bed_path, header=None, sep="\t", compression=compression)
    start = df[1].to_numpy()
    end = df[2].to_numpy()

    # get peak sizes
    peak_sizes = end - start
    good_index = []
    for i, size in enumerate(peak_sizes):
        if size < max_len:
            good_index.append(i)

    # create dictionary for dataframe with filtered peaks
    data = {}
    for i in range(len(df.columns)):
        data[i] = df[i].to_numpy()[good_index]

    # create new dataframe
    df_new = pd.DataFrame(data)

    # save dataframe with fixed width window size to a bed file
    df_new.to_csv(output_path, sep="\t", header=None, index=False)



def enforce_constant_size(bed_path, output_path, window):
    """Generate a bed file where all peaks have same size centered on original peak
    Parameters
    ----------
    bed_path : <str>
        The path to the unfiltered bed file downloaded from
        ENCODE.
    output_path : <str>
        Filtered bed file name with peak size clipped to specified
        window size.
    window : <int>
        Size in bps which to clip the peak sequences.
    Returns
    ----------
    None
    Example
    ----------
    >>> window = 200
    >>> bed_path = './ENCFF252PLM.bed.gz'
    >>> output_path = './pos_'+str(window)+'.bed'
    >>> enforce_constant_size(pos_path, pos_bed_path, window, compression='gzip')
    """
    assert isinstance(window, int) and window > 0, "Enter positive integer window size."
    assert os.path.exists(bed_path), "No such bed file."

    # set up the compression argument
    if bed_path.split(".")[-1] == "gz" or bed_path.split(".")[-1] == "gzip":
        compression = "gzip"
    else:
        compression = None

    # load bed file
    df = pd.read_table(bed_path, header=None, sep="\t", compression=compression)
    # chrom = df[0].to_numpy().astype(str)  # TODO: unused variable
    start = df[1].to_numpy()
    end = df[2].to_numpy()

    # calculate center point and create dataframe
    middle = np.round((start + end) / 2).astype(int)
    left_window = np.round(window / 2).astype(int)
    right_window = window - left_window

    # calculate new start and end points
    start = middle - left_window
    end = middle + right_window

    # create dictionary for dataframe
    data = {}
    for i in range(len(df.columns)):
        data[i] = df[i].to_numpy()
    data[1] = start
    data[2] = end
    for i in range(3, df.shape[1]):
        data[i] = df[i].to_numpy()

    # create new dataframe
    df_new = pd.DataFrame(data)

    # save dataframe with fixed width window size to a bed file
    df_new.to_csv(output_path, sep="\t", header=None, index=False)



def get_bed_coords(bed_path):
    assert os.path.exists(bed_path), "No such bed file."

    # set up the compression argument
    if bed_path.split(".")[-1] == "gz" or bed_path.split(".")[-1] == "gzip":
        compression = "gzip"
    else:
        compression = None

    # load bed file
    df = pd.read_table(bed_path, header=None, sep="\t", compression=compression)
    chrom = df[0].to_numpy().astype(str)  # TODO: unused variable
    start = df[1].to_numpy().astype(str)
    end = df[2].to_numpy().astype(str)
    strand = df[5].to_numpy().astype(str)

    coords = []
    for a in zip(chrom, start, end, strand):
        coords.append("%s:%s-%s(%s)" % (a[0], a[1], a[2], a[3]))
    return np.array(coords)


def filter_nonsense_sequences(sequences):
    """Filter sequences with N.
    Parameters
    ----------
    sequences : <numpy.ndarray>
        A numpy vector of sequence strings.
    Returns
    -------
    filter_sequences : <numpy.ndarray>
        The parsed sequences from the input fasta file as a numpy array of sequences.
    good_index : <numpy.ndarray>
        A numpy array of indices corresponding to sequences without nonsense 'N'
        entries.
    Example
    -------
    >>> print(sequences)
    ['GGCTGAAATGGCCACTGGAA' 'ACGCTCTCTCATCAAGTGGT' 'GCAGAANANCGAACACCAAC'
    'NNCNNCANCNACNGGGGAAC' 'GCCTAGTCCAGACATAATTC']
    >>> print(filter_nonsense_sequences(sequences))
    (array(['GGCTGAAATGGCCACTGGAA', 'ACGCTCTCTCATCAAGTGGT',
            'GCCTAGTCCAGACATAATTC'], dtype='<U20'), array([0, 1, 4]))
    """

    # filter sequences if contains at least one 'N' character
    good_index = []
    filter_sequences = []
    for i, seq in enumerate(sequences):
        if "N" not in seq.upper():
            good_index.append(i)
            filter_sequences.append(seq)
    return np.array(filter_sequences), np.array(good_index)



def _is_gzipped(filepath):
    """Return `True` if the file is gzip-compressed.
    This function does not depend on the suffix. Instead, the magic number of the file
    is compared to the GZIP magic number `1f 8b`. 
    """
    with open(filepath, "rb") as f:
        return f.read(2) == b"\x1f\x8b"


def parse_fasta(filepath):
    """Parse FASTA file into arrays of descriptions and sequence data.
    Parameters
    ----------
    filepath : Path-like
        FASTA file to parse. Can be gzip-compressed.
    Returns
    -------
    tuple of two numpy arrays
        The first array contains the sequences, and the second array contains the name
        of each sequence. These arrays have the same length in the first dimension.
    """
    # FASTA format described here.
    # https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=BlastHelp
    descriptions = []
    sequences = []
    prev_line_was_sequence = False

    with open(filepath, "rt") as f:  
        for line in f:
            line = line.strip()
            # handle blank lines
            if not line:
                continue
            is_description = line.startswith(">")
            if is_description:
                description = line[1:].strip()  # prune ">" char
                descriptions.append(description)
                prev_line_was_sequence = False
            else:  # is sequence data
                sequence = line.upper()
                if prev_line_was_sequence:
                    # This accounts for sequences that span multiple lines.
                    sequences[-1] += sequence
                else:
                    sequences.append(sequence)
                prev_line_was_sequence = True
    return np.array(sequences), np.array(descriptions)

    

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


def match_gc_content(pos_one_hot, neg_one_hot, neg_pos_ratio=1):

    N, L, A = pos_one_hot.shape
    gc_pos = np.sum(np.sum(pos_one_hot[:,:,[1,2]], axis=2), axis=1)/L
    gc_neg = np.sum(np.sum(neg_one_hot[:,:,[1,2]], axis=2), axis=1)/L
    print('    Average GC content for positive sequences: %.3f'%(np.mean(gc_pos)))
    print('    Average GC content for negative sequences: %.3f'%(np.mean(gc_neg)))

    pos_index = np.argsort(gc_pos)
    neg_index = np.argsort(gc_neg)
    num_neg = len(neg_index)
    num_pos = len(pos_index)

    match_index = []
    if num_neg > num_pos:
        k = 0
        status = True
        for i in pos_index:
            for j in range(k, num_neg):
                if gc_pos[i] < gc_neg[neg_index[j]]:
                    if k > num_neg:
                        status = False
                        break
                    else:
                        # print("%.2f vs %.2f"%(gc_pos[i], gc_neg[neg_index[j]]))
                        match_index.append(neg_index[j])
                        k = j+1
                        break
            if not status:
                break

    remainder = int(num_pos*neg_pos_ratio) - len(match_index)
    print('    Found %d GC-matched sequences.'%(len(match_index)))
    if remainder > 0:
        print('    Adding %d more random negative sequences.'%(remainder))
        remain_index = np.array(list(set(range(num_neg)) - set(match_index)))
        index = np.random.permutation(len(remain_index))[:remainder]     
        # index = np.argsort(gc_neg[remain_index])[::-1]
        for n in remain_index[index[:remainder]]:
            match_index.append(n)
                
    match_index = np.array(match_index)
    print('    Average GC content for sub-sampled negative sequences: %.3f'%(np.mean(gc_neg[match_index])))

    return neg_one_hot[match_index], match_index

def split_dataset_random(x, y, coords, valid_frac=0.1, test_frac=0.2):
    """split dataset into training, cross-validation, and test set"""

    def split_index(num_data, valid_frac, test_frac):
        # split training, cross-validation, and test sets

        train_frac = 1 - valid_frac - test_frac
        cum_index = np.array(
            np.cumsum([0, train_frac, valid_frac, test_frac]) * num_data
        ).astype(int)
        shuffle = np.random.permutation(num_data)
        train_index = shuffle[cum_index[0] : cum_index[1]]
        valid_index = shuffle[cum_index[1] : cum_index[2]]
        test_index = shuffle[cum_index[2] : cum_index[3]]

        return train_index, valid_index, test_index

    # split training, cross-validation, and test sets
    train_index, valid_index, test_index = split_index(len(x), valid_frac, test_frac)

    # split dataset
    train = (x[train_index], y[train_index, :], coords[train_index])
    valid = (x[valid_index], y[valid_index, :], coords[valid_index])
    test = (x[test_index], y[test_index, :], coords[test_index])

    return train, valid, test


def split_dataset_chrom(x, y, coords, valid_chr, test_chr):

    def parse_held_out_chrom(x, y, coords, split_chroms):
        """extract data according to chromosome"""
        # get list of chromosomes for each data entry
        chrom_list = np.array([n.split(':')[0] for n in coords])

        # parse held out chrom data
        index = []
        for chrom_index in split_chroms:
            index.append(np.where(chrom_index == chrom_list)[0])
        index = np.concatenate(index)
        return x[index], y[index], coords[index]

    def remove_held_out_chrom(x, y, coords, remove_chroms):
        """remove data according to chromosome"""

        # get list of chromosomes for each data entry
        chrom_list = np.array([n.split(':')[0] for n in coords])

        # filter out held out chroms in data
        chrom_filt = np.copy(chrom_list)
        x_filt = np.copy(x)
        y_filt = np.copy(y)
        coords_filt = np.copy(coords)
        for chrom_index in remove_chroms:
            index = np.where(chrom_index == chrom_filt)[0]
            chrom_filt = np.delete(chrom_filt, index, axis=0)
            x_filt = np.delete(x_filt, index, axis=0)
            y_filt = np.delete(y_filt, index, axis=0)
            coords_filt = np.delete(coords_filt, index, axis=0)
        return x_filt, y_filt, coords_filt

    # get list of held out chromosomes 
    valid_chr = ['chr'+chrom_index for chrom_index in valid_chr.split(',')]
    test_chr  = ['chr'+chrom_index for chrom_index in test_chr.split(',')]

    # generate dataset
    test  = parse_held_out_chrom(x, y, coords, test_chr)
    valid = parse_held_out_chrom(x, y, coords, valid_chr)
    train = remove_held_out_chrom(x, y, coords, valid_chr+test_chr)

    return train, valid, test


#------------------------------------------------------------
# bedtools functions
#------------------------------------------------------------



def bedtools_getfasta(
    bed_path, genome_path, output_path, strand=True, bedtools_exe="bedtools"
):
    """Extract DNA sequences from a fasta file based on feature coordinates.
    Wrapper around `bedtools getfasta`. This function was made to
    work with bedtools version 2.27.1. It is not guaranteed to work
    with other versions. It is not even guaranteed to work with version 2.27.1, but
    it could and probably will.
    Parameters
    ----------
    genome_path : str, Path-like
        path to reference genome in fasta format.
    output_path : str, Path-like
        Output FASTA file.
    bed_path : str, Path-like
        BED/GFF/VCF file of ranges to extract from `input_fasta`.
    strand : bool
        Force strandedness. If the feature occupies the antisense
        strand, the squence will be reverse complemented.
    exe_call : Path-like
        The path to the `bedtools` executable. By default, uses `bedtools` in `$PATH`.
    Returns
    -------
    Instance of `subprocess.CompletedProcess`.
    """
    args = [str(bedtools_exe), "getfasta"]
    if strand:
        args.append("-s")
    args.extend(
        ["-fi", str(genome_path), "-bed", str(bed_path), "-fo", str(output_path)]
    )
    try:
        return subprocess.run(
            args, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
    except subprocess.CalledProcessError as e:
        raise subprocess.SubprocessError(e.stderr.decode()) from e



def bedtools_intersect(
    a, b, output_path, write_a=True, nonoverlap=True, bedtools_exe="bedtools",
):
    """Report overlaps between two feature files.
    This is an incomplete wrapper around `bedtools intersect` version 2.27.1.
    The set of arguments here does not include all of the command-line arguments.
    Parameters
    ----------
    a : Path-like
        First feature file <bed/gff/vcf/bam>.
    b : Path-like
        Second feature file <bed/gff/vcf/bam>.
    output_bedfile : Path-like
        Name of output file. Can be compressed (`.bed.gz`).
    write_a : bool
        Write the original entry in `a` for each overlap.
    write_b : bool
        Write the original entry in `b` for each overlap.
    invert_match : bool
        Only report those entries in `a` that have no overlaps with `b`.
    bedtools_exe : Path-like
        The path to the `bedtools` executable. By default, uses `bedtools` in `$PATH`.
    Returns
    -------
    Instance of `subprocess.CompletedProcess`.
    """
    args = [str(bedtools_exe), "intersect"]
    if write_a:
        args.append("-wa")
    if nonoverlap:
        args.append("-v")
    args.extend(["-a", str(a), "-b", str(b)])

    try:
        process = subprocess.run(
            args, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        if not process.stdout:
            raise subprocess.SubprocessError(
                f"empty stdout, aborting. stderr is {process.stderr.decode()}"
            )
        with open(output_path, mode="wb") as f:  # type: ignore
            f.write(process.stdout)
        return process
    except subprocess.CalledProcessError as e:
        raise subprocess.SubprocessError(e.stderr.decode()) from e





if __name__ == "__main__":

    main()

