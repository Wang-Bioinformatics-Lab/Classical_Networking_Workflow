import os
import pandas as pd
import argparse
import yaml
import glob
import numpy as np

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

    try:
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

    except Exception as e:
        print("Error while attempting to read USIs:", e)
        usi_list = []

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
    parser.add_argument('--private_spectra_file_list', default=None, help='Input TSV file of private spectra files')

    args = parser.parse_args()

    # Load USI list and extract filenames
    usi_list = load_usi_list(args.input_usi_information)
    usi_filenames = []

    try:
        usi_filenames = [os.path.basename(usi.split(":")[2]) for usi in usi_list]
        
        # Remove USIs which we have more than once, keeping the first occurrence
        first_occurrence_index = {filename: idx for idx, filename in enumerate(usi_filenames) if filename not in locals().get('first_occurrence_index', {})}
        usi_list = [usi_list[idx] for filename, idx in first_occurrence_index.items()]
        usi_filenames = [usi_filenames[idx] for filename, idx in first_occurrence_index.items()]

        # Remove USIs of files which have identical names to private spectra files
        private_spectra = pd.read_csv(args.private_spectra_file_list, sep='\t')
        private_spectra['Filename'] = private_spectra['Filename'].apply(os.path.basename)
        private_spectra_file_names = private_spectra['Filename'].tolist() if not private_spectra.empty else []

        usi_list, usi_filenames = map(list, zip(*[
            (usi, filename) for usi, filename in zip(usi_list, usi_filenames) 
            if filename not in private_spectra_file_names
        ]))

    except Exception as e:
        print("Error during USI processing:", e)
        usi_list = []
        usi_filenames = []
        pass

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

        # Convert back to lists after filtering
        usi_list = list(usi_list)
        usi_filenames = list(usi_filenames)

        # Getting all the input files
        all_input_files = usi_filenames.copy()

        all_input_files += spectra_filenames


    if not os.path.exists(args.input_metadata):
        # This is likely not valid, lets skip it
        input_metadata = pd.DataFrame()
    else:
        input_metadata = _load_metadata(args.input_metadata)

    usi_list_with_redu_data_df = pd.DataFrame()

    if args.include_redu == "Yes" and len(usi_list) > 0:
        try:
            # Downloading the ReDU Metadata
            redu_df = get_redu_metadata()

            # We'll match up the data from the dataset accession and the filename
            usi_list_with_redu_data_df = match_usi_to_redu_metadata(usi_list, redu_df)

            # If public files are not in ReDU DataSource public but not in redu as info
            # Find filenames in usi_filenames that are not present in the 'filename' column of usi_list_with_redu_data_df
            missing_filenames = [filename for filename in usi_filenames if filename not in usi_list_with_redu_data_df['filename'].values]
            # If there are any missing filenames, add them with all values being NaN
            if missing_filenames:
                # Create a DataFrame for the missing entries with NaN for other columns
                missing_data = {col: np.nan for col in usi_list_with_redu_data_df.columns}
                missing_data['filename'] = missing_filenames
                missing_data['ATTRIBUTE_DataSource'] = 'public_no_ReDU'
                
                # Convert the dictionary to a DataFrame, repeating NaNs for each missing filename
                new_entries_df = pd.DataFrame(missing_data)
                
                # Concatenate the new entries with the original DataFrame
                usi_list_with_redu_data_df = pd.concat([usi_list_with_redu_data_df, new_entries_df], ignore_index=True)


        except Exception as e:
            print("There was a problem with downloading or matching the ReDU metadata:", e)


    # combine private metadata for public data with redu metadata for public files
    ##############################################################################
    if 'filename' in input_metadata.columns and len(usi_list_with_redu_data_df) > 0:
        try:

            input_metadata_redu_enrichment = input_metadata[~input_metadata['filename'].isin(private_spectra_file_names)]
            input_metadata_redu_enrichment = input_metadata_redu_enrichment.drop(columns=input_metadata_redu_enrichment.columns.intersection(usi_list_with_redu_data_df.columns).drop('filename'))

            if input_metadata_redu_enrichment['filename'].isin(usi_list_with_redu_data_df['filename']).any():


                # merge the input_metadata with usi_list_with_redu_data_df
                usi_list_with_redu_data_df = pd.merge(input_metadata_redu_enrichment,
                                                    usi_list_with_redu_data_df, on="filename", how="right")
                
                # remove rows with filenames from input_metadata that have matches in usi_list_with_redu_data_df
                input_metadata = input_metadata[input_metadata['filename'].isin(private_spectra_file_names)]

        except Exception as e:
            print("There was a problem with merging private and redu metadata for public data:", e)

    # combine private metadata for private raw data with metadata of public files
    ##############################################################################
    if len(usi_list_with_redu_data_df) > 0 and len(private_spectra_file_names) > 0:

        try:

            if not 'filename' in input_metadata.columns:
                print('No metadata file with filename column found.')
                input_metadata = pd.DataFrame({'filename': []})

            # add rows for private spectra files that are not in metadata table
            existing_filenames = set(input_metadata['filename'])
            new_filenames = [fname for fname in private_spectra_file_names if fname not in existing_filenames]

            if len(new_filenames) > 0:
                print('Some spectra files not present in metadata files.')

                missing_data = {col: np.nan for col in input_metadata.columns}
                missing_data['filename'] = new_filenames

                new_entries_df = pd.DataFrame(missing_data)
                input_metadata = pd.concat([input_metadata, new_entries_df], ignore_index=True)

            input_metadata['ATTRIBUTE_DataSource'] = 'private'
            output_metadata = pd.concat([input_metadata, usi_list_with_redu_data_df], ignore_index=True, sort=False)

        except Exception as e:
            print("There was a problem with merging private metadata for private raw data with metadata of public files:", e)

    # if we have only metadata for public files just return it (no private raw data-metadat to add)
    ##############################################################################
    elif len(usi_list_with_redu_data_df) > 0:
        output_metadata = usi_list_with_redu_data_df

    # if we only have private metadata but no metadata for public files return this
    ##############################################################################
    elif len(input_metadata) > 0:
        output_metadata = input_metadata

    # if we have no metadata return empty dataframe
    else:
        output_metadata = pd.DataFrame()

    # If we want per file grouping
    if args.per_file_grouping == "Yes":

        try:
            new_metadata = pd.DataFrame()

            new_metadata["filename"] = all_input_files
            new_metadata["ATTRIBUTE_PerFileGrouping"] = new_metadata["filename"]

            # joining it with the input_metadata
            if "filename" in output_metadata:
                output_metadata = pd.merge(output_metadata, new_metadata, on="filename", how="left")
            else:
                output_metadata = new_metadata

        except Exception as e:
            print("There was a problem with per file grouping:", e)


    # Remove all columns where we have no data
    if len(output_metadata) > 0:
        try:
            output_metadata = output_metadata.loc[:, ~output_metadata.apply(lambda col: col.isna() | (col == 'missing value')).all()]
        except Exception as e:
            print("There was a problem with removing columns with no data:", e)


    # Filtering the output metadata
    try:
        output_metadata = _filter_metadata(output_metadata, all_input_files)
    except:
        pass

    # outputing metadata
    output_metadata.to_csv(args.merged_metadata, sep="\t", index=False)


if __name__ == '__main__':
    main()