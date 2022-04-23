import os
import h5py
import pandas as pd
import numpy as np

#-------------------------------------------------------------------------------
# File handling
#-------------------------------------------------------------------------------

def make_directory(path, foldername, verbose=1):
  """make a directory"""

  if not os.path.isdir(path):
    os.mkdir(path)
    print("making directory: " + path)

  outdir = os.path.join(path, foldername)
  if not os.path.isdir(outdir):
    os.mkdir(outdir)
    print("making directory: " + outdir)
  return outdir



def is_gzipped(filepath):
  """Return `True` if the file is gzip-compressed.
  This function does not depend on the suffix. Instead, the magic number of the file
  is compared to the GZIP magic number `1f 8b`. 
  """
  with open(filepath, "rb") as f:
    return f.read(2) == b"\x1f\x8b"


def save_dataset_hdf5(filepath, train, valid, test, coords=False):
  """Saves an h5 dataset based on train--> x, y, coords"""
  with h5py.File(filepath, "w") as fout:
    fout.create_dataset("x_train", data=train[0], compression="gzip")
    fout.create_dataset("y_train", data=train[1], compression="gzip")
    fout.create_dataset("x_valid", data=valid[0], compression="gzip")
    fout.create_dataset("y_valid", data=valid[1], compression="gzip")
    fout.create_dataset("x_test", data=test[0], compression="gzip")
    fout.create_dataset("y_test", data=test[1], compression="gzip")
    if coords:
      fout.create_dataset("coords_train", data=train[2].astype("S"), compression="gzip")
      fout.create_dataset("coords_valid", data=valid[2].astype("S"), compression="gzip")
      fout.create_dataset("coords_test", data=test[2].astype("S"), compression="gzip")



#-------------------------------------------------------------------------------
# Parsing
#-------------------------------------------------------------------------------
    


def get_chrom_sizes(chrom_sizes_path):
  """Load chrom sizes in dictionary,
  Parameters
    ----------
    chrom_sizes_path : str
        path to chrom.sizes 
  """

  df = pd.read_csv(chrom_sizes_path, delimiter='\t', header=None)
  chroms = df[0].values
  sizes = df[1].values
  chrom_sizes = OrderedDict()
  for chrom, size in zip(df[0].values, df[1].values):
    chrom_sizes[chrom] = int(size)
  return chrom_sizes




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




#-------------------------------------------------------------------------------
# Munging
#-------------------------------------------------------------------------------



def filter_chr_list(chrom_sizes, ignore_chr_y, ignore_auxiliary_chr):
  """determine which chromosomes to keep for dataset"""

  # regular expression for standard chromosomes
  primary_re = re.compile("chr\\d+$")

  if ignore_auxiliary_chr:
    keep_chr = []
    for chr in chrom_sizes.keys():
      if primary_re.match(chr):
        keep_chr.append(chr)
    keep_chr.append('chrX')
    keep_chr.append('chrY')
  else:
    keep_chr = list(chrom_sizes.keys())

  if ignore_chr_y:
    keep_chr.remove('chrY')

  return keep_chr


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



def convert_fasta_to_onehot(fasta_file, alphabet='ACGT', uncertain_N=True):
  """open fasta and save one-hot in dictionary with coordinate as key"""

  seq_vecs = OrderedDict()
  for line in open(fasta_file):
    if line[0] == ">":
      name = line[1:].rstrip()
    else:
      seq_vecs[name] = convert_one_hot([line.rstrip().upper()], alphabet, uncertain_N)

  return seq_vecs



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


