#!/usr/bin/python


import sys
import getopt
import os
import json
import argparse
import molecular_network_filtering_library
from collections import defaultdict
import pandas as pd

def main():
    parser = argparse.ArgumentParser(description='Add Component Number to network')
    parser.add_argument('clustersummary', help='clustersummary')
    parser.add_argument('networking_pairs_results_file', help='networking_pairs_results_file')
    parser.add_argument('output_network_edges', help='output_network_edges')

    args = parser.parse_args()

    # Reading the cluster summary
    cluster_summary_df = pd.read_csv(args.clustersummary, sep="\t")

    # Reading the pairs file
    pairs_df = pd.read_csv(args.networking_pairs_results_file, sep="\t")

    # Adding in other columsn that may have been dropped
    if "DeltaMZ" not in pairs_df.columns:
        # Adding delta m/z
        
        pairs_list = pairs_df.to_dict(orient="records")

        for pair in pairs_list:
            cluster1 = pair["CLUSTERID1"]
            cluster2 = pair["CLUSTERID2"]

            mz1 = cluster_summary_df[cluster_summary_df["cluster index"] == cluster1]["parent mass"].values[0]
            mz2 = cluster_summary_df[cluster_summary_df["cluster index"] == cluster2]["parent mass"].values[0]

            deltamz = mz1 - mz2

            pair["DeltaMZ"] = deltamz

        pairs_df = pd.DataFrame(pairs_list)

    if "ComponentIndex" in pairs_df.columns:
        pairs_df.to_csv(args.output_network_edges, sep="\t", index=False)
        exit(0)

    # Saving a copy in the output
    pairs_df.to_csv(args.output_network_edges, sep="\t", index=False)

    # Loading and calculating
    G = molecular_network_filtering_library.loading_network(args.output_network_edges, hasHeaders=True)
    
    if G == None:
        exit(0)

    molecular_network_filtering_library.output_graph_with_headers(G, args.output_network_edges)


if __name__ == "__main__":
    main()