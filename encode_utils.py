import os 
import numpy as np 
import pandas as pd 
import subprocess
import utils

def generate_exp_list_from_metadata(
  data_path,
  metadata_path,
  save_path=None,
  criteria={},
  label_list=['Assay', 'Experiment target', 'Biosample term name', 'Experiment accession'],
  download=True,
  exp_accession_list=[],
):
  """Generate subset of metadata table and sample file for further processing.
  Parameters
  ----------
  data_path : str
      dataset directory with files
  metadata_path : str
      metadata containing experiments of interest
  save_path : str
      path where the bed file will be saved
  criteria : dict
      dictionary of column, value pairs to use in making a selection
  exp_accession_list : list
      list of experiments to select, if empty select all in the metadata table
  """

  def make_label(df_entry, label_list):
    """Generate a unique label"""
    items = [
      str(c.values[0])
      for c in [df_entry[entry] for entry in label_list]
    ]
    return "_".join(items).replace(" ", "-")


  # load meta data
  metadata = pd.read_csv(metadata_path, sep="\t")

  # get list of experiment accessions from metadata if not given
  if not exp_accession_list:
    print("Generating experiment list file from all of the metadata table")
    exp_accession_list = list(set(metadata["Experiment accession"]))

  # get the name and filepath of experiments that satisfy criteria
  summary = []
  for i, exp_accession in enumerate(exp_accession_list):

    # get experiment details
    exp_df = metadata[(metadata["Experiment accession"] == exp_accession)]

    # generate name
    exp_name = make_label(exp_df, label_list)
   
    # filter by criteria
    if criteria:
      for name, value in criteria.items():
        exp_df = exp_df[(exp_df[name] == value)]

    if len(exp_df) == 0:
      print("    Warning: criteria not met for " + exp_accession)
    else:
      print("    Processed: " + exp_accession)
      exp_df = exp_df.iloc[[0]]
      
      # check to see if file exists
      if 'bed' in criteria['File format']:
        ext = 'bed.gz'
      elif 'bigWig' == criteria['File format']:
        ext = 'bigWig'
      filepath = check_filepath(data_path, exp_df["File accession"].values[0], ext)
      print(filepath)
      if not filepath: 
        if download:  # download file from link provided in metatable
          print('  Downloading %s: %s'%(exp_name, exp_df['File download URL'].values[0]))
          download_path = os.path.join(data_path, exp_df["File accession"].values[0] + ext)
          cmd = 'wget -O ' + download_path + ' ' + exp_df['File download URL'].values[0]
          subprocess.call(cmd, shell='True')
          filepath = os.path.join(data_path, exp_df["File accession"].values[0], ext)
          #if 'bed' in ext:
          #  cmd = 'gunzip ' + download_path
          #  subprocess.call(cmd, shell ='True')

      # store results
      if filepath:
        summary.append([make_label(exp_df, label_list), filepath])

  if save_path:
    # save to file
    with open(save_path, "w") as fout:
      for line in summary:
        fout.write("{}\t{}\n".format(line[0], line[1]))

  print('%d out of %d experiments processed in: %s'%(len(summary), len(exp_accession_list), save_path))
  return summary




def generate_exp_list_from_dir(
    data_path,
    save_path,
    ext='bed', # 'bed' or 'gz'
):
  """Generate a sample file for further processing based on files in a directory
  Parameters
  ----------
  data_path : str
      dataset directory with files
  save_path : str
      path where the bed file will be saved
  ext : str
      extension of files to be included in experiment list
  """

  filepaths = []
  names = []
  for f in os.listdir(data_path):
    splits = f.split('.')
    if ext in splits:
      names.append(splits[0])
      filepaths.append(os.path.join(data_path, f))

  if save_path:
    # save to file
    with open(save_path, "w") as fout:
      for name, filepath in zip(names, filepaths):
        fout.write("{}\t{}\n".format(name, filepath))

  print('%d experiments processed in: %s'%(len(names), save_path))
  return [names, filepaths]



def get_path_from_metadata(
  data_path,
  criteria={},
  metadata_name='metadata.tsv',
  ext='bed',
):
  """Generate subset of metadata table and sample file for further processing.
  Parameters
  ----------
  data_path : str
      dataset directory with files
  metadata_path : str
      metadata containing experiments of interest
  criteria : dict
      dictionary of column, value pairs to use in making a selection
  """

  # load meta data
  metadata = pd.read_csv(os.path.join(data_path, metadata_name), sep="\t")
  exp_accession_list = list(set(metadata["Experiment accession"]))
  
  # get the name and filepath of experiments that satisfy criteria
  names = []
  filepaths = []
  for i, exp_accession in enumerate(exp_accession_list):
    exp_df = metadata[(metadata["Experiment accession"] == exp_accession)]

    # filter by criteria
    if criteria:
      for name, value in criteria.items():
        exp_df = exp_df[(exp_df[name] == value)]

    if len(exp_df) != 0:
      sizes = []
      for j in range(len(exp_df)):
        name = exp_df.iloc[[j]]["File accession"].values[0]
        filepath = check_filepath(data_path, name, ext)
        # store results
        if filepath:
          names.append(name)
          filepaths.append(filepath)
          sizes.append(int(exp_df.iloc[[j]]["Size"].values[0]))
      sizes = np.array(sizes)
      break

  if names:
    index = np.argsort(sizes)[::-1]
    return np.array(names)[index], np.array(filepaths)[index]
  else:
    print("No file matching criteria found.")
    return None, None



def check_filepath(data_path, filename, ext='bed'):
  """Generate path where the file is found.
  Parameters
  ----------
  data_path : str
      dataset path to file
  filename : str
      name of file
  """

  filepath = os.path.join(data_path, filename + '.' + ext)
  if os.path.isfile(filepath):
    return filepath

  elif os.path.isfile(filepath + ".gz"):
    return filepath + ".gz"
  else:
    return ""

