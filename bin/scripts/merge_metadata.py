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
        except:
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
    usi_df = pd.DataFrame({"usi": usi_list})
    usi_df["filename"] = usi_df["usi"].apply(lambda x: x.split(":")[2])
    usi_df["basename"] = usi_df["filename"].apply(lambda x: os.path.basename(x))
    usi_df["ATTRIBUTE_DatasetAccession"] = usi_df["usi"].apply(lambda x: x.split(":")[1])
    usi_df["merge_column"] = usi_df["ATTRIBUTE_DatasetAccession"] + ":" + usi_df["basename"]

    redu_df["basename"] = redu_df["filename"].apply(lambda x: os.path.basename(x))
    redu_df["merge_column"] = redu_df["ATTRIBUTE_DatasetAccession"] + ":" + redu_df["basename"]

    merged_df = pd.merge(usi_df, redu_df, on="merge_column", how="left")

    # we should cleanup filenames
    try:
        merged_df["filename"] = merged_df["filename_x"]

        # lets drop filename_x and filename_y
        merged_df = merged_df.drop(columns=["filename_x", "filename_y"])
    except:
        pass
    
    # Clean up accession
    try:
        merged_df["ATTRIBUTE_DatasetAccession"] = merged_df["ATTRIBUTE_DatasetAccession_x"]
        merged_df = merged_df.drop(columns=["ATTRIBUTE_DatasetAccession_x", "ATTRIBUTE_DatasetAccession_y"])
    except:
        pass
    
    # TODO: make sure appropriate columns have attribute in redu metdata
    all_columns = list(merged_df.columns)
    # Include these columns SampleType SampleTypeSub1 NCBITaxonomy YearOfAnalysis SampleCollectionMethod   SampleExtractionMethod   MassSpectrometer  IonizationSourceAndPolarity    ChromatographyAndPhase  BiologicalSex  UBERONBodyPartName    HealthStatus   DOIDCommonName  Country  HumanPopulationDensity 
    redu_columns = ["SampleType", "SampleTypeSub1", "NCBITaxonomy", "YearOfAnalysis", "SampleCollectionMethod", "SampleExtractionMethod", "MassSpectrometer", "IonizationSourceAndPolarity", "ChromatographyAndPhase", "BiologicalSex", "UBERONBodyPartName", "HealthStatus", "DOIDCommonName", "Country", "HumanPopulationDensity"]

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


    # Getting all the input files is available
    all_input_files = None

    if args.spectra_folder is not None:
        all_input_files = []
        
        all_input_files += glob.glob(args.spectra_folder + "/*.mzML")
        all_input_files += glob.glob(args.spectra_folder + "/*.mzml")
        all_input_files += glob.glob(args.spectra_folder + "/*.mzXML")
        all_input_files += glob.glob(args.spectra_folder + "/*.mzxml")
        all_input_files += glob.glob(args.spectra_folder + "/*.mgf")
        all_input_files += glob.glob(args.spectra_folder + "/*.MGF")

        # getting basename for all input files
        all_input_files = [os.path.basename(x) for x in all_input_files]

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
            input_metadata = pd.merge(input_metadata, merged_df, on="filename", how="left")
        else:
            input_metadata = merged_df

    # If we want per file grouping
    if args.per_file_grouping == "Yes":
        new_metadata = pd.DataFrame()

        new_metadata["filename"] = all_input_files
        new_metadata["ATTRIBUTE_PerFileGrouping"] = new_metadata["filename"]

        # joining it with the input_metadata
        if "filename" in input_metadata:
            input_metadata = pd.merge(input_metadata, new_metadata, on="filename", how="left")
        else:
            input_metadata = new_metadata

    # Filtering the output metadata
    try:
        input_metadata = _filter_metadata(input_metadata, all_input_files)
    except:
        pass

    # outputing metadata
    input_metadata.to_csv(args.merged_metadata, sep="\t", index=False)


if __name__ == '__main__':
    main()