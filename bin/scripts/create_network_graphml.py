#!/usr/bin/python


import sys
import getopt
import os
import molecular_network_filtering_library
import networkx as nx
import argparse

def main():
    parser = argparse.ArgumentParser(description='Creating Clustering Info Summary')
    parser.add_argument('input_clusterinfo_summary', help='input_clusterinfo_summary')
    parser.add_argument('input_pairs', help='input_pairs')
    parser.add_argument('input_library_matches', help='input_library_matches')
    parser.add_argument('output_graphml', help='output_graphml')

    args = parser.parse_args()

    #Doing other filtering
    G = molecular_network_filtering_library.loading_network(args.input_pairs, hasHeaders=True)
    molecular_network_filtering_library.add_clusterinfo_summary_to_graph(G, args.input_clusterinfo_summary)
    molecular_network_filtering_library.add_library_search_results_to_graph(G, args.input_library_matches)

    nx.write_graphml(G, sys.argv[4], infer_numeric_types=True)




if __name__ == "__main__":
    main()
