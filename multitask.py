import os, sys, gzip, re, h5py
from collections import OrderedDict
import numpy as np
import subprocess
import utils



def process_data(
  exp_list_path,
  prefix_save_path,
  genome_path,
  chrom_sizes_path,
  seq_len=1000,
  merge_overlap_len=200,
  valid_chr='12,13',
  test_chr='7,8',
  ignore_chr_y=False,
  ignore_auxiliary_chr=False,
  uncertain_N=True,
  alphabet='ACGT',
  blacklist_path=None,        # path to blacklist bed file
  standard_coord=False,
  save_tfr=False,
):

  # generate a merged multitask bed file and activity table (outputs: prefix_merged.bed and prefix_activity.tsv)
  print("Merging bed files from: %s"%(exp_list_path))
  generate_multitask_bed(
    exp_list_path=exp_list_path,
    prefix_save_path=prefix_save_path,
    seq_len=seq_len,
    merge_overlap_len=merge_overlap_len,
    chrom_sizes_path=chrom_sizes_path,
    ignore_chr_y=ignore_chr_y,
    ignore_auxiliary_chr=ignore_auxiliary_chr,
    standard_coord=standard_coord,
  )
  print("Saved to: %s_merged.bed and %s_activity.tsv"%(prefix_save_path, prefix_save_path))


  # remove ENCODE blacklist sites if provided
  if blacklist_path:
    print("removing unmappable regions defined by: %s"%(blacklist_path))

    bed_utils.bedtools_intersect(
      a=prefix_save_path+"_merged.bed",
      b=blacklist_path,
      output_path= prefix_save_path+"_merged.bed", 
      write_a=False, nonoverlap=True
    )

    bed_utils.bedtools_intersect(
      a=prefix_save_path+"_activity.tsv",
      b=blacklist_path,
      output_path= prefix_save_path+"_activity.tsv", 
      write_a=False, nonoverlap=True
    )

  # convert merged bed file into fasta file (eg. prefix.fa)
  print("Extracting sequences from reference genome: %s"%(genome_path))
  bed_utils.bedtools_getfasta(
    bed_path=prefix_save_path+"_merged.bed", 
    genome_path=genome_path,
    output_path=prefix_save_path+".fa",
  )

  # load sequences from fasta file
  print("Processing files -- converting to one-hot and generating labels")
  seq_vecs = utils.convert_fasta_to_onehot(prefix_save_path+".fa", alphabet, uncertain_N)
  os.remove(prefix_save_path+".fa")

  # load scores
  seq_targets, target_labels = parse_targets(prefix_save_path+"_activity.tsv")

  # align and construct input matrix
  seqs, targets, coords = align_seqs_targets(seq_vecs, seq_targets, sort=False)

  # split dataset by chromosome
  print("Splittig dataset into:")
  print("    validation: %s"%(valid_chr))
  print("    test: %s"%(test_chr))
  train, valid, test = utils.split_dataset_chrom(seqs, targets, coords, valid_chr, test_chr)

  # save data
  if save_tfr:
    filepath = prefix_save_path+'_train.tfr'
    print("Saving to: %s"%(filepath))
    utils.write_TFRecord_dataset(filepath, train[0], train[1], train[2])
    filepath = prefix_save_path+'_valid.tfr'
    print("Saving to: %s"%(filepath))
    utils.write_TFRecord_dataset(filepath, valid[0], valid[1], valid[2])
    filepath = prefix_save_path+'_test.tfr'
    print("Saving to: %s"%(filepath))
    utils.write_TFRecord_dataset(filepath, test[0], test[1], test[2])
  else:    
    filepath = prefix_save_path+'.h5'
    print("Saving to: %s"%(filepath))
    utils.save_dataset_hdf5(filepath, train, valid, test, coords=True)
  print('%d train sequences'%(len(train[0])))
  print('%d valid sequences'%(len(valid[0])))
  print('%d test sequences '%(len(test[0])))



def generate_multitask_bed(
  exp_list_path,
  prefix_save_path=None,
  seq_len=1000,
  merge_overlap_len=200,
  chrom_sizes_path=None,
  ignore_chr_y=False,
  ignore_auxiliary_chr=False,
  standard_coord=False,
):

  """Merge multiple bed files to select sample sequence regions with at least one
  peak.
  This function outputs a .bed file in the specified directory containing seven
  columns: chromosome, sequence start, sequence end, name, score, strand, and indexes
  of experiments that have a peak in this region.
  Parameters
  ----------
  exp_list_path: str
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
  standard_coord: bool, optional
      if True, then saves chr:start-end(strand), else it's chr_start_end
  Returns
  -------
  None
  Examples
  --------
  >>> multitask_bed_generation(
      example_file,chr_lens_file='/data/hg38.chrom.size',
      seq_len=1000,merge_overlap_len=200,prefix_save_path='/data/bed')
  """

  if not exp_list_path:
    raise Exception(
      "Must provide file labeling the targets and providing BED file paths."
    )

  # read in targets and assign them indexes into the db
  target_names = []
  target_bedpaths = []
  target_index = []
  for line in open(exp_list_path):
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
  if chrom_sizes_path is not None:
    chr_lens = utils.get_chrom_sizes(chrom_sizes_path)
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
    bed_utils.bedtools_sort(chr_file_dict[chr_key], chrom_sbed)
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
      peak_act = activity_set(a[6])
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
          mpeaks = merge_peaks(
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
      mpeaks = merge_peaks(
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
    if standard:
      peak_id = "%s:%s-%s(%s)" % (a[0], a[1], a[2], a[5])
    else:
      peak_id = "%s_%s_%s" % (a[0], a[1], a[2])

    # construct full activity vector
    peak_act = [0] * len(target_names)
    for ai in a[6].split(","):
      if ai != ".":
        peak_act[int(ai)] = 1

    # print line
    cols = [peak_id] + peak_act
    print("\t".join([str(c) for c in cols]), file=final_act_out)
  final_act_out.close()





#------------------------------------------------------------
# Useful functions/classes for multi-task
#------------------------------------------------------------



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
    peak_mids = [find_midpoint(self.start, self.end)]
    peak_mids.append(find_midpoint(peak2.start, peak2.end))

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


def find_midpoint(start, end):
  """ Find the midpoint coordinate between start and end """
  mid = (start + end) / 2
  return int(mid)


def merge_peaks(peaks, peak_size, merge_overlap_len, chr_len):
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


def activity_set(act_cs):
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


def parse_targets(targets_file):
  """Load activity file"""
  seq_targets = OrderedDict()
  for line in open(targets_file):
    a = line.split()
    try:
      seq_targets[a[0]] = np.array([int(a[i]) for i in range(1, len(a))])
    except Exception:
      target_labels = a
  return seq_targets, target_labels


def align_seqs_targets(seq_vecs, seq_targets, sort=True):
  """Align and construct input matrix"""
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
    
