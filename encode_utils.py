import os
import pandas as pd
import numpy as np


def generate_exp_bed_list(
    data_path,
    metadata_path,
    exp_bed_list_path=None,
    criteria={},
    label_list=['Assay', 'Experiment target', 'Biosample term name', 'Experiment accession'],
    exp_accession_list=[],
):
    """Generate subset of metadata table and sample file for further processing
    Parameters
    ----------
    data_path : str
        dataset directory with files
    metadata_path : str
        metadata containing experiments of interest
    exp_bed_list_path : str
        path where the bed file will be saved
    criteria : dict
        dictionary of column, value pairs to use in making a selection
    exp_accession_list : list
        list of experiments to select, if empty select all in the metadata table
    """

    def make_label(df_entry):
        """Generate a unique label"""
        items = [
            str(c.values[0])
            for c in [df_entry[entry] for entry in label_list]
        ]
        return "_".join(items).replace(" ", "-")



    def check_filepath(data_path, filename):
        """Generate path where the file is found"""
        filepath = os.path.join(data_path, filename + ".bed")
        if os.path.isfile(filepath):
            return filepath

        elif os.path.isfile(filepath + ".gz"):
            return filepath + ".gz"
        else:
            return ""

    # load meta data
    metadata = pd.read_csv(metadata_path, sep="\t")

    # get list of experiment accessions from metadata if not given
    if not exp_accession_list:
        print("Generating experiment list file from all of the metadata table")
        exp_accession_list = list(set(metadata["Experiment accession"]))
    
    # get the name and filepath of experiments that satisfy criteria
    summary = []
    for i, exp_accession in enumerate(exp_accession_list):
        exp_df = metadata[(metadata["Experiment accession"] == exp_accession)]

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
            filepath = check_filepath(data_path, exp_df["File accession"].values[0])

            # store results
            if filepath:
                summary.append([make_label(exp_df), filepath])
    
    if exp_bed_list_path:
        # save to file
        with open(exp_bed_list_path, "w") as fout:
            for line in summary:
                fout.write("{}\t{}\n".format(line[0], line[1]))

    print('%d out of %d experiments processed in: %s'%(len(summary), len(exp_accession_list), exp_bed_list_path))
    return summary

