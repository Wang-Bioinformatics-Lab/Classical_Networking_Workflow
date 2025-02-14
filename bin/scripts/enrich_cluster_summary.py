#!/usr/bin/python

import sys
import getopt
import os
import json
import argparse
import uuid
from collections import defaultdict
import csv
import re
import pandas as pd
import glob

def main():
    parser = argparse.ArgumentParser(description='Running library search parallel')
    parser.add_argument('input_clustersummary', help='spectrafolder')
    parser.add_argument('input_filtered_pairs', help='output folder for parameters')
    parser.add_argument('input_library_matches', help='output folder for parameters')
    parser.add_argument('clustersummary_with_network', help='output folder for parameters')
    args = parser.parse_args()

    # Loading Cluster Summary
    df_cluster_summary = pd.read_csv(args.input_clustersummary, sep='\t')

    # Loading Filtered Pairs
    df_filtered_pairs = pd.read_csv(args.input_filtered_pairs, sep='\t')

    # Loading Library Matches
    try:
        df_library_matches = pd.read_csv(args.input_library_matches, sep='\t')
    except:
        df_library_matches = pd.DataFrame(columns=['#Scan#', 'Compound_Name', 'MQScore'])

    # Get mapping from node  to component
    node_to_component = {}
    for index, row in df_filtered_pairs.iterrows():
        component = int(row['ComponentIndex'])

        node1 = int(row['CLUSTERID1'])
        node2 = int(row['CLUSTERID2'])

        node_to_component[node1] = component
        node_to_component[node2] = component

    # Adding component to cluster_summary
    df_cluster_summary['component'] = df_cluster_summary['cluster index'].map(node_to_component)
    df_cluster_summary['component'] = df_cluster_summary['component'].fillna(-1)
    # Making sure column is an integer
    df_cluster_summary['component'] = df_cluster_summary['component'].astype(int)

    # Sort library matches by score
    df_library_matches = df_library_matches.sort_values(by=['MQScore'], ascending=False)

    # Grouping by the scan
    df_library_matches = df_library_matches.groupby('#Scan#').head(1)

    # Merging in library matches
    df_cluster_summary = pd.merge(df_cluster_summary, df_library_matches, how='left', left_on='cluster index', right_on='#Scan#')

    #df_cluster_summary['library matches'] = df_cluster_summary['cluster index'].map(df_library_matches.set_index('#Scan#')['Compound_Name'])

    # Exporting summary with components
    df_cluster_summary.to_csv(args.clustersummary_with_network, sep='\t', index=False)

if __name__ == "__main__":
    main()
