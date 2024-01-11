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

    # Lets make sure the filename is good
    input_df["filename"] = input_df["filename"].apply(lambda x: os.path.basename(x))
    input_df = input_df[input_df['filename'].notna() & input_df['filename'].ne('')]

    return input_df

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
    usi_df = pd.DataFrame({"usi": usi_list})

    if not 'usi' in list(redu_df.columns):
        redu_df['usi'] = redu_df['filename']
        redu_df['usi'] = redu_df['usi'].str[2:]
        redu_df['usi'] = 'mzspec:' + redu_df['usi']
        redu_df['usi'] = redu_df['usi'].apply(lambda x: x.replace('/', ':', 1))

    merged_df = pd.merge(usi_df, redu_df, on="usi", how="left")

    merged_df.loc[merged_df['filename'].isna() | merged_df['filename'].eq(''), 'filename'] = merged_df['usi'].apply(lambda x: x.split(":")[2])
    merged_df["filename"] = merged_df['usi'].apply(lambda x: x.split(":")[2])
        
    # TODO: make sure appropriate columns have attribute in redu metdata
    all_columns = list(merged_df.columns)
    # Include these columns SampleType SampleTypeSub1 NCBITaxonomy YearOfAnalysis SampleCollectionMethod   SampleExtractionMethod   MassSpectrometer  IonizationSourceAndPolarity    ChromatographyAndPhase  BiologicalSex  UBERONBodyPartName    HealthStatus   DOIDCommonName  Country  HumanPopulationDensity 
    redu_columns = ["SampleType", "SampleTypeSub1", "NCBITaxonomy", "YearOfAnalysis", "SampleCollectionMethod", "SampleExtractionMethod", "MassSpectrometer", "IonizationSourceAndPolarity", "ChromatographyAndPhase", "BiologicalSex", "UBERONBodyPartName", "HealthStatus", "DOIDCommonName", "Country", "HumanPopulationDensity", "DataSource"]

    for column in all_columns:
        if column in redu_columns:
            # rename the column with ATTRIBUTE_ prefix using rename mechanism
            merged_df = merged_df.rename(columns={column: "ATTRIBUTE_" + column})

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

    if not os.path.exists(args.input_metadata):
        # This is likely not valid, lets skip it
        input_metadata = pd.DataFrame()
    else:
        input_metadata = _load_metadata(args.input_metadata)

    if args.include_redu == "Yes":
        # we need to load the usi list
        usi_list = load_usi_list(args.input_usi_information)

        # Downloading the ReDU Metadata
        redu_df = get_redu_metadata()

        # We'll match up the data from the dataset accession and the filename
        merged_df = match_usi_to_redu_metadata(usi_list, redu_df)

        # Merging the metadata with the input
        if "filename" in input_metadata:
            if "ATTRIBUTE_DataSource" in merged_df:
                input_metadata["ATTRIBUTE_DataSource"] = 'user input data'
            input_metadata = pd.concat([input_metadata, merged_df], ignore_index=True)
        else:
            input_metadata = merged_df

    # If we want per file grouping
    if args.per_file_grouping == "Yes":
        all_input_files = []
        
        all_input_files += glob.glob(args.spectra_folder + "/*.mzML")
        all_input_files += glob.glob(args.spectra_folder + "/*.mzml")
        all_input_files += glob.glob(args.spectra_folder + "/*.mzXML")
        all_input_files += glob.glob(args.spectra_folder + "/*.mzxml")
        all_input_files += glob.glob(args.spectra_folder + "/*.mgf")
        all_input_files += glob.glob(args.spectra_folder + "/*.MGF")

        new_metadata = pd.DataFrame()

        # getting basename for all input files
        all_input_files = [os.path.basename(x) for x in all_input_files]

        new_metadata["filename"] = all_input_files
        new_metadata["ATTRIBUTE_PerFileGrouping"] = new_metadata["filename"]

        # joining it with the input_metadata
        if "filename" in input_metadata:
            input_metadata = pd.merge(input_metadata, new_metadata, on="filename", how="left")
        else:
            input_metadata = new_metadata

    # outputing metadata
    input_metadata.to_csv(args.merged_metadata, sep="\t", index=False)


if __name__ == '__main__':
    main()