import os
import pandas as pd
import argparse
import yaml
import glob

def _load_metadata(input_filename):
    # look at extension
    if input_filename.endswith(".tsv"):
        input_df = pd.read_csv(input_filename, sep="\t")
    elif input_filename.endswith(".csv"):
        input_df = pd.read_csv(input_filename, sep=",")
    elif input_filename.endswith(".xlsx"):
        input_df = pd.read_excel(input_filename)
    else:
        input_df = pd.read_csv(input_filename, sep=None)

    # Removing rows where nan is in the filename
    input_df = input_df.dropna(subset=["filename"])

    # Lets make sure the filename is good
    input_df["filename"] = input_df["filename"].apply(lambda x: os.path.basename(x))

    return input_df

def _filter_metadata(input_df, filelist):
    # Making sure is basename
    filelist = set([os.path.basename(x) for x in filelist])

    # Making sure filename is basename
    input_df["filename"] = input_df["filename"].apply(lambda x: os.path.basename(x))

    # Filtering the metadata
    filtered_input_df = input_df[input_df["filename"].isin(filelist)]

    # We can handle a special case where there is no extensions and the result is empty
    if len(filtered_input_df) == 0:
        try:
            # create a mapping of the filelist without extension to with extension
            filelist_mapping = {os.path.splitext(x)[0]: x for x in filelist}

            # lets try mapping the filename column using this filelist mapping
            remapped_filtered_input_df = input_df.copy()
            remapped_filtered_input_df["filename"] = remapped_filtered_input_df["filename"].apply(lambda x: filelist_mapping.get(x, x))
            filtered_input_df = remapped_filtered_input_df[remapped_filtered_input_df["filename"].isin(filelist)]
        except Exception as e:
            print(f"Error during filename remapping: {e}")
            pass

    return filtered_input_df

def load_usi_list(input_filename):
    # Checking the file extension
    if input_filename.endswith(".yaml"):
        # Loading yaml file
        parameters = yaml.load(open(input_filename), Loader=yaml.SafeLoader)
        usi_list = parameters["usi"].split("\n")
    elif input_filename.endswith(".tsv"):
        df = pd.read_csv(input_filename, sep="\t")
        usi_list = df["usi"].tolist()

    # Cleaning USI list
    usi_list = [usi.lstrip().rstrip() for usi in usi_list]

    return usi_list

def get_redu_metadata():
    url = "https://redu.gnps2.org/dump"

    metadata_df = pd.read_csv(url, sep="\t")

    return metadata_df

def match_usi_to_redu_metadata(usi_list, redu_df):

    usi_df = pd.DataFrame({"USI": usi_list})

    # Merge tables
    merged_df = pd.merge(usi_df, redu_df, on="USI", how="inner")
    
    # get all columns from redu metadata other than USI as list in one line
    redu_columns = [col for col in merged_df.columns if col not in ["USI", "filename", "ATTRIBUTE_DatasetAccession"]]

    for column in redu_columns:
        # rename the column with ATTRIBUTE_ prefix using rename mechanism
        merged_df = merged_df.rename(columns={column: "ATTRIBUTE_" + column})

    # Drop unnecessary columns
    merged_df = merged_df.drop(columns=["filename"])

    # Extract filename from USI
    merged_df["filename"] = merged_df["USI"].apply(lambda x: x.split(":")[2])
    merged_df["filename"] = merged_df["filename"].apply(lambda x: os.path.basename(x))

    return merged_df

def main():
    # Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('input_metadata')
    parser.add_argument('input_usi_information')
    parser.add_argument('merged_metadata')

    parser.add_argument('--include_redu', default="No", help='Include redu metadata integration')
    parser.add_argument('--per_file_grouping', default="No", help='Treat each file as a group')
    parser.add_argument('--spectra_folder', default=None, help='Input Folder of Spectra Files')

    args = parser.parse_args()

    # Load USI list and extract filenames
    usi_list = load_usi_list(args.input_usi_information)
    usi_filenames = [os.path.basename(usi.split(":")[2]) for usi in usi_list]
    
    #remove USIs which we have more than once keeping the last occurance (which should be the one we keep after download)
    last_occurrence_index = {filename: idx for idx, filename in enumerate(usi_filenames)}
    usi_list = [usi_list[idx] for filename, idx in last_occurrence_index.items()]
    usi_filenames = [usi_filenames[idx] for filename, idx in last_occurrence_index.items()]

    # Getting all the input files
    all_input_files = usi_filenames.copy()

    if args.spectra_folder is not None:
        spectra_files = []
        spectra_files += glob.glob(args.spectra_folder + "/*.mzML")
        spectra_files += glob.glob(args.spectra_folder + "/*.mzml")
        spectra_files += glob.glob(args.spectra_folder + "/*.mzXML")
        spectra_files += glob.glob(args.spectra_folder + "/*.mzxml")
        spectra_files += glob.glob(args.spectra_folder + "/*.mgf")
        spectra_files += glob.glob(args.spectra_folder + "/*.MGF")

        # getting basename for all input files
        spectra_filenames = [os.path.basename(x) for x in spectra_files]

        # TODO: figure out a way we can take care of duplicated filenames between public and private data. Will likely have to be done already durng download since spectra_folder private files are likely replaced by downlaoded public files

        all_input_files += spectra_filenames


    if not os.path.exists(args.input_metadata):
        # This is likely not valid, lets skip it
        input_metadata = pd.DataFrame()
    else:
        input_metadata = _load_metadata(args.input_metadata)

    usi_list_with_redu_data_df = pd.DataFrame()

    if args.include_redu == "Yes" and len(usi_list) > 0:

        # Downloading the ReDU Metadata
        redu_df = get_redu_metadata()

        # We'll match up the data from the dataset accession and the filename
        usi_list_with_redu_data_df = match_usi_to_redu_metadata(usi_list, redu_df)


    # combine private metadata for public data with redu metadata for public files
    if 'filename' in input_metadata.columns and len(usi_list_with_redu_data_df) > 0:

        if input_metadata['filename'].isin(usi_list_with_redu_data_df['filename']).any():

            # merge the input_metadata with usi_list_with_redu_data_df
            usi_list_with_redu_data_df = pd.merge(input_metadata.drop(columns=input_metadata.columns.intersection(usi_list_with_redu_data_df.columns).drop('filename')),
                                                  usi_list_with_redu_data_df, on="filename", how="right")

            # remove rows with filenames from input_metadata that have matches in usi_list_with_redu_data_df
            input_metadata = input_metadata[~input_metadata['filename'].isin(usi_list_with_redu_data_df['filename'])]

    # combine private metadata for private raw data with metadata of public files
    if len(usi_list_with_redu_data_df) > 0 and len(input_metadata) > 0:

        input_metadata['ATTRIBUTE_DataSource'] = 'private'
        output_metadata = pd.concat([input_metadata, usi_list_with_redu_data_df], ignore_index=True, sort=False)

    # if we have only metadata for public files just return it (no private raw data-metadat to add)
    elif len(usi_list_with_redu_data_df) > 0:
        output_metadata = usi_list_with_redu_data_df

    # if we only have private metadata but no metadata for public files return this
    elif len(input_metadata) > 0:
        output_metadata = input_metadata

    # if we have no metadata return empty dataframe
    else:
        output_metadata = pd.DataFrame()

    # If we want per file grouping
    if args.per_file_grouping == "Yes":
        new_metadata = pd.DataFrame()

        new_metadata["filename"] = all_input_files
        new_metadata["ATTRIBUTE_PerFileGrouping"] = new_metadata["filename"]

        # joining it with the input_metadata
        if "filename" in output_metadata:
            output_metadata = pd.merge(output_metadata, new_metadata, on="filename", how="left")
        else:
            output_metadata = new_metadata

    # Filtering the output metadata
    try:
        output_metadata = _filter_metadata(output_metadata, all_input_files)
    except:
        pass

    # outputing metadata
    output_metadata.to_csv(args.merged_metadata, sep="\t", index=False)


if __name__ == '__main__':
    main()