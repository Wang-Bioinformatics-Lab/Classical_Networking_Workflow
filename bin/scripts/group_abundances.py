#!/usr/bin/python


import sys
import getopt
import os
import argparse
import statistics
import pandas as pd
from tqdm import tqdm

from collections import defaultdict

def create_attribute_group_list(metadata_df):
    # Determining all the groups we want to calculate over columsn with prefix ATTRIBUTE_
    all_attributes = [x for x in metadata_df.columns if x.startswith("ATTRIBUTE_")]

    all_attribute_groups = []

    for attribute in all_attributes:
        # Getting the groups in each attribute and creating a list of them
        attribute_groups = metadata_df[attribute].unique().tolist()

        # Creating a dictionary for each attribute group
        for attribute_group in attribute_groups:
            attribute_group_dict = {}
            attribute_group_dict["attribute"] = attribute
            attribute_group_dict["group"] = attribute_group
            all_attribute_groups.append(attribute_group_dict)

    return all_attributes, all_attribute_groups


# This function calculates all the group counds for all the relevant columns
def calculate_groups_metadata(clustersummary_df, clusterinfo_df, metadata_df):
    # Cleaning the filenames
    clusterinfo_df["#Filename"] = clusterinfo_df["#Filename"].apply(lambda x: os.path.basename(x))
    metadata_df["filename"] = metadata_df["filename"].apply(lambda x: os.path.basename(x))

    # Folding in the metadata into clusterinfo
    clusterinfo_df = clusterinfo_df.merge(metadata_df, left_on="#Filename", right_on="filename", how="left")

    # First lets group the cluster info by cluster index
    grouped_clusterinfo_df = clusterinfo_df.groupby("#ClusterIdx")

    # loop through all the clusters
    cluster_summary_list = clustersummary_df.to_dict(orient="records")

    # Getting the attributes
    all_attributes, all_attribute_groups = create_attribute_group_list(metadata_df)

    for cluster in tqdm(cluster_summary_list):
        # filter for the grouped list
        cluster_index = cluster["cluster index"]
        clusterinfo_per_group_df = grouped_clusterinfo_df.get_group(cluster_index)

        # TODO: We can likely speed this up with pandas operations

        # Attribute_group
        for attribute_group in all_attribute_groups:
            #print(attribute_group)

            # filtering the data
            group_count = len(clusterinfo_per_group_df[clusterinfo_per_group_df[attribute_group["attribute"]] == attribute_group["group"]])
            group_column = "{}:GNPSGROUP:{}".format(attribute_group["attribute"], attribute_group["group"])

            cluster[group_column] = group_count
    
        # Adding the cluster information for which group membership it is a part of
        for attribute in all_attributes:
            # Finding all groups in the attribute
            all_groups = set(clusterinfo_per_group_df[attribute])

            # Converting to string
            all_groups = [str(x) for x in all_groups]

            cluster[attribute] = ",".join(all_groups)

    return pd.DataFrame(cluster_summary_list)


def calculate_attribute_metadata(clustersummary_df, metadata_df):
    # Getting the attributes
    all_attributes, all_attribute_groups = create_attribute_group_list(metadata_df)

    cluster_summary_list = clustersummary_df.to_dict(orient="records")

    for cluster in tqdm(cluster_summary_list):
        # filter for the grouped list
        cluster_index = cluster["cluster index"]

        for attribute in all_attributes:
            # Finding all groups in the attribute via list comprehension
            current_attribute_groups = [x for x in all_attribute_groups if x["attribute"] == attribute]

            positive_groups = []

            for attribute_group in current_attribute_groups:
                group_column = "{}:GNPSGROUP:{}".format(attribute_group["attribute"], attribute_group["group"])
                if cluster[group_column] > 0:
                    positive_groups.append(attribute_group["group"])

            # Converting to string
            positive_groups = [str(x) for x in positive_groups]

            cluster[attribute] = ",".join(positive_groups)

    return pd.DataFrame(cluster_summary_list)



def main():
    parser = argparse.ArgumentParser(description='Creates enriched cluster info summary')
    parser.add_argument('input_clusterinfo_file', help='input_clusterinfo_file')
    parser.add_argument('input_clusterinfosummary_file', help='input_clusterinfosummary_file')
    parser.add_argument('input_metadata', help='input_group_mapping_filename')
    parser.add_argument('output_clusterinfosummary_filename', help='output_clusterinfosummary_filename')
    args = parser.parse_args()

    # Loading Data
    clustersummary_df = pd.read_csv(args.input_clusterinfosummary_file, sep="\t")
    clustersinfo_df = pd.read_csv(args.input_clusterinfo_file, sep="\t")

    try:
        metadata_df = pd.read_csv(args.input_metadata, sep="\t")
    except:
        metadata_df = pd.DataFrame()

    if len(metadata_df) > 0:
        # Enriching metadata group counts
        clustersummary_df = calculate_groups_metadata(clustersummary_df, clustersinfo_df, metadata_df)

        # Now we will write out the attribute column
        clustersummary_df = calculate_attribute_metadata(clustersummary_df, metadata_df)

    # Writing out the file
    clustersummary_df.to_csv(args.output_clusterinfosummary_filename, sep="\t", index=False)

    exit(0)


    """Loading group filenames"""
    group_to_files, files_to_groups = load_group_mapping(args.input_group_mapping_filename)
    print("Loaded Group Mapping")
    
    
    print("Loaded Cluster Summary")

    attribute_to_groups = load_attribute_mapping(args.input_attribute_mapping_filename)

    #Calculating the spectrum counts per group
    cluster_to_group_counts = defaultdict(lambda: defaultdict(lambda: 0))
    cluster_to_files = defaultdict(set)
    cluster_to_RT = defaultdict(list)
    line_count = 0

    



    # for line in open(args.input_clusterinfo_file):
    #     line_count += 1
    #     if line_count == 1:
    #         continue
    #     if line_count % 10000 == 0:
    #         print(line_count)

    #     splits = line.rstrip().split("\t")
    #     cluster_index = splits[0]
    #     filename = os.path.basename(splits[1])
    #     rt = float(splits[6])

    #     group_membership = files_to_groups[filename]
    #     cluster_to_files[cluster_index].add(filename)
    #     cluster_to_RT[cluster_index].append(rt)

    #     for group in group_membership:
    #         cluster_to_group_counts[cluster_index][group] += 1

    # print(len(cluster_summary_list))

    

    print("Setting up grouping", len(group_to_files.keys()))
    for cluster_summary_object in cluster_summary_list:
        cluster_index = cluster_summary_object["cluster index"]
        for group in group_to_files:
            group_count = 0
            if group in cluster_to_group_counts[cluster_index]:
                group_count = cluster_to_group_counts[cluster_index][group]
            cluster_summary_object[group] = group_count

        for attribute in attribute_to_groups:
            groups_to_include = []
            for group in attribute_to_groups[attribute]:
                if group in cluster_summary_object:
                    if cluster_summary_object[group] > 0:
                        groups_to_include.append(group)

            cluster_summary_object[attribute] = ",".join(groups_to_include).replace("GNPSGROUP:", "")


    print("Default Attributes")
    calculate_default_attributes(cluster_summary_list, group_to_files.keys())

    # print("calculate_cluster_file_stats")
    # calculate_cluster_file_stats(cluster_summary_list, cluster_to_files, mangled_mapping)

    # print("rt stats")
    # calculate_rt_stats(cluster_summary_list, cluster_to_RT)

    # print("populate_network_component")
    # populate_network_component(cluster_summary_list, args.input_networking_pairs)

    # print("calculate_ancillary_information")
    # calculate_ancillary_information(cluster_summary_list, params_object["task"][0])    

    # print("populate_network_identifications")
    # populate_network_identifications(cluster_summary_list, args.input_library_search)

    # Write pandas
    df = pd.DataFrame(cluster_summary_list)
    df.to_csv(args.output_clusterinfosummary_filename, sep="\t", index=False)




if __name__ == "__main__":
    main()
