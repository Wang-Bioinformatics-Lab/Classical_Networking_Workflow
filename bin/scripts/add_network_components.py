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
    parser.add_argument('networking_pairs_results_file', help='networking_pairs_results_file')
    parser.add_argument('output_network_edges', help='output_network_edges')

    args = parser.parse_args()

    # Early Exit
    pairs_df = pd.read_csv(args.networking_pairs_results_file, sep="\t")
    if "ComponentIndex" in pairs_df.columns:
        pairs_df.to_csv(args.output_network_edges, sep="\t", index=False)
        exit(0)

    # Loading and calculating
    G = molecular_network_filtering_library.loading_network(args.networking_pairs_results_file, hasHeaders=True)
    if G == None:
        exit(0)

    molecular_network_filtering_library.output_graph_with_headers(G, args.output_network_edges)


if __name__ == "__main__":
    main()