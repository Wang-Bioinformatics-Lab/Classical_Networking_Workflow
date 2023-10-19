import os
import pandas as pd
import argparse
import yaml


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
    
    # TODO: make sure appropriate columns have attribute in redu metdata
    all_columns = list(merged_df.columns)
    redu_columns = ["SampleType"]

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
            input_metadata = pd.merge(input_metadata, merged_df, on="filename", how="left")
        else:
            input_metadata = merged_df

    # outputing metadata
    input_metadata.to_csv(args.merged_metadata, sep="\t", index=False)


if __name__ == '__main__':
    main()