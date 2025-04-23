import os
import pandas as pd
import argparse
import numpy as np

def rewrite_falcon_mgf(mgf_input, mgf_output):
    with open(mgf_input, 'r') as infile, open(mgf_output, 'w') as outfile:
        for line in infile:
            if line.startswith("CLUSTER="):
                try:
                    cluster_num = int(line.strip().split("=")[1])
                    outfile.write(f"SCANS={cluster_num + 1}\n")
                except ValueError:
                    outfile.write(line)
            else:
                outfile.write(line)

def main():
    parser = argparse.ArgumentParser(description='Summarizing Falcon Results')
    parser.add_argument('falcon_clusters', help='falcon_clusters')
    parser.add_argument('falcon_mgf', help='falcon_mgf')
    #parser.add_argument('output_summary_folder', help='output_summary_folder')
    args = parser.parse_args()

    clusterinfo_df = pd.read_csv(args.falcon_clusters, sep=',', comment='#')

    print(args)
    print(clusterinfo_df)

    clusterinfo_df = clusterinfo_df.sort_values(by='cluster', key=lambda x: x.replace(-1, np.inf))

    # Filtering out not in clusters data
    clusterinfo_df = clusterinfo_df[clusterinfo_df["cluster"] != -1]


    # Grouping by cluster
    grouped_cluster_df = clusterinfo_df.groupby(["cluster"])
    cluster_summary_list = []
    for cluster, cluster_group_df in grouped_cluster_df:
        #TODO :Read these from mgf, as the representative is a medoid

        cluster_count = len(cluster_group_df)
        cluster_mz = cluster_group_df["precursor_mz"].mean()
        cluster_rt = cluster_group_df["retention_time"].mean()
        cluster_charge = cluster_group_df["precursor_charge"].mean()
        # adjust the col name to map the classical MN wokflow
        output_dict = {}
        output_dict["number of spectra"] = cluster_count
        output_dict["parent mass"] = cluster_mz
        output_dict["RTMean"] = cluster_rt
        output_dict["precursor charge"] = cluster_charge
        output_dict["cluster index"] = cluster[0] + 1

        cluster_summary_list.append(output_dict)

    # Creating a cluster summary
    cluster_summary_df = pd.DataFrame(cluster_summary_list)
    cluster_summary_df.to_csv("clustersummary.tsv", sep='\t', index=False)

    # Creating cluster info
    clusterinfo_df["filename"] = clusterinfo_df["identifier"].apply(lambda x: x.split(":")[2] + ".mzML")
    clusterinfo_df["scan"] = clusterinfo_df["identifier"].apply(lambda x: x.split(":")[-1])

    # Rename relevant columns to MS-Cluster format
    clusterinfo_df = clusterinfo_df.rename(columns={
        "filename": "#Filename",
        "cluster": "#ClusterIdx",
        "scan": "#Scan",
        "precursor_mz": "#ParentMass",
        "precursor_charge": "#Charge",
        "retention_time": "#RetTime"
    })

    # Just to prevent other processes in the workflow raise error
    clusterinfo_df["#PrecIntensity"] = 0

    # Select required columns
    clusterinfo_df = clusterinfo_df[[
        "#ClusterIdx", "#Filename", "#Scan", "#ParentMass", "#Charge", "#RetTime",  "#PrecIntensity"
    ]]
    clusterinfo_df.to_csv("clusterinfo.tsv", sep='\t', index=False)

    #TODO: Rewriting MGF files
    rewrite_falcon_mgf(args.falcon_mgf, "specs_ms.mgf")
    # TODO: Maybe make this compatible with FBMN, since its already clustered.

if __name__ == "__main__":
    main()