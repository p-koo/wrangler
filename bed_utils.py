import subprocess

def bed_getfasta(
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



def bed_intersect(
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



def bed_sort(
  bed_path, output_path, bedtools_exe="bedtools"
):
  """Sorts bed file.
  Parameters
  ----------
  bed_path : str, Path-like
      Input bed filepath.
  output_path : str, Path-like
      Sorted bed filepath.
  exe_call : Path-like
      The path to the `bedtools` executable. By default, uses `bedtools` in `$PATH`.
  Returns
  -------
  Instance of `subprocess.CompletedProcess`.
  """
  args = [str(bedtools_exe), "sortBed"]
  args.extend(
    ["-i", str(bed_path), ">", str(output_path)]
  )
  try:
    return subprocess.run(
      args, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
  except subprocess.CalledProcessError as e:
    raise subprocess.SubprocessError(e.stderr.decode()) from e



def filter_blacklist(input_bed_path, output_bed_path, blacklist_bed_path):
"""filter out blacklisted (i.e. unmappable regions)"""

  subprocess.call(
    "bedtools intersect -a {} -b {} -v > {}".format(
      input_bed_path, blacklist_bed_path, output_bed_path
    ),
    shell=True,
  )


