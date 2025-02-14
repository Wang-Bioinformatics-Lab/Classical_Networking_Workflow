#!/usr/bin/python

import sys
import os
import glob
import argparse
import ming_spectrum_library
import pandas as pd

def write_params(input_spectra_folder, tool_path, output_params_filename,
                 min_cluster_size=2,
                 pm_tolerance=2.0, 
                 fragment_tolerance=0.5, 
                 min_peak_intensity=0.0,
                 window_filter=1,
                 precursor_filter=1):

    # listing all files
    all_mgf_files = glob.glob(os.path.join(input_spectra_folder, "*.mgf"))
    all_mzxml_files = glob.glob(os.path.join(input_spectra_folder, "*.mzXML"))
    all_mzml_files = glob.glob(os.path.join(input_spectra_folder, "*.mzML"))

    all_spectrum_files = all_mgf_files + all_mzxml_files + all_mzml_files

    # Sort these filenames
    all_spectrum_files.sort()

    with open(output_params_filename, "w") as params_file:
        params_file.write("CLUSTER_MODEL=LTQ_TRYP\n")
        params_file.write("CLUST_RANK_FILTER=6\n")
        params_file.write("CORRECT_PM=no\n")
        params_file.write("GUESS_CHARGE=no\n")
        params_file.write("MIN_SPECTRUM_QUALITY=0.0\n")

        # Basic Tolerances
        params_file.write("TOLERANCE_PEAK={}\n".format(fragment_tolerance))
        params_file.write("TOLERANCE_PM={}\n".format(pm_tolerance))
        
        # Window Filtering
        if int(window_filter) == 1:
            params_file.write("CLUST_RANK_FILTER=6\n")
            params_file.write("RANK_FILTER=6\n")
            params_file.write("RANK_FILTER_RADIUS=50.0\n")
        
        params_file.write("MIN_PEAK_INT={}\n".format(min_peak_intensity))

        # Filtering
        params_file.write("FILTER_PRECURSOR_WINDOW={}\n".format(precursor_filter))
        params_file.write("FILTER_STDDEV_PEAK_INT=0.0\n")


        params_file.write("CLUSTER_MIN_SIZE={}\n".format(min_cluster_size))

        params_file.write("EXE_DIR={}\n".format(tool_path))
        params_file.write("INPUT_SPECS_MS={}\n".format(";".join(all_spectrum_files)))
    
def main():
    # Parse arguments
    parser = argparse.ArgumentParser(description='MSCluster Wrapper')
    parser.add_argument('input_spectra_folder', help='Input Spectra Folder')
    parser.add_argument('tool_dir', help='tool_dir')
    parser.add_argument('output_spectra_folder', help='Output Spectra Folder')
    parser.add_argument('final_output_folder', help='final_output_folder')

    parser.add_argument('--min_cluster_size', default="2", help='min_cluster_size')

    parser.add_argument('--pm_tolerance', default="2.0", help='pm_tolerance')
    parser.add_argument('--fragment_tolerance', default="0.5", help='fragment_tolerance')

    # Filters
    parser.add_argument('--min_peak_intensity', default="0.0", help='min_peak_intensity')
    parser.add_argument('--window_filter', default="1", help='window_filter')
    parser.add_argument('--precursor_filter', default="1", help='precursor_filter')
    
    args = parser.parse_args()

    # Writing the parameters
    parameters_filename = "mscluster.params"
    write_params(args.input_spectra_folder, args.tool_dir, parameters_filename, 
                 min_cluster_size=args.min_cluster_size,
                 pm_tolerance=args.pm_tolerance, 
                 fragment_tolerance=args.fragment_tolerance, 
                 min_peak_intensity=args.min_peak_intensity,
                 window_filter=args.window_filter,
                 precursor_filter=args.precursor_filter)

    # Running the data
    mainspecnets_binary = os.path.join(args.tool_dir, "main_specnets")

    cmd = "{} {} -ll 0 -f mscluster".format(mainspecnets_binary, parameters_filename)
    ret_code = os.system(cmd)

    if ret_code != 0:
        exit(1)

    # Perform a rewrite of the mgf file, this does not need to be done, because original ordering is fine
    specs_mgf_filename = os.path.join(args.output_spectra_folder, "specs_ms.mgf")
    spectrum_collection = ming_spectrum_library.SpectrumCollection(specs_mgf_filename)
    spectrum_collection.load_from_mgf()

    # TODO: we need to make sure that there are empty spectra
    output_mgf_filename = os.path.join(args.final_output_folder, "specs_ms.mgf")

    # Making sure the scans are continuous
    spectrum_collection.make_scans_contiguous()
    
    # Saving out the MGF
    spectrum_collection.save_to_mgf(open(output_mgf_filename, "w"), renumber_scans=False)

    # Creating the clusterinfo file
    path_to_clusterinfo = os.path.join(args.tool_dir, "clusterinfo")
    clusterinfo_file = os.path.join(args.final_output_folder, "clusterinfo.tsv")
    clustersummary_file = os.path.join(args.final_output_folder, "clustersummary.tsv")
    cmd = "{} --outfile {} --out-summary-file {}".format(path_to_clusterinfo, clusterinfo_file, clustersummary_file)
    os.system(cmd)

    # Loading data to be refined
    cluster_summary_df = pd.read_csv(clustersummary_file, sep="\t")
    clusterinfo_df = pd.read_csv(clusterinfo_file, sep="\t")

    # rewrite output to only include the cluster summary that have size greater than the parameter
    cluster_summary_df = cluster_summary_df[cluster_summary_df["number of spectra"] >= int(args.min_cluster_size)]
    included_clusters_set = set(cluster_summary_df["cluster index"].tolist())

    # Selecting only the clusters in the summary file
    clusterinfo_df = clusterinfo_df[clusterinfo_df["#ClusterIdx"].isin(included_clusters_set)]

    # Reading the cluster info file and figuring out the mean RT
    clusterinfo_filtered_df = clusterinfo_df[["#ClusterIdx", "#RetTime"]]
    clusterinfo_filtered_df = clusterinfo_filtered_df.rename(columns={"#ClusterIdx": "cluster index", "#RetTime": "retention time"})
    
    # Make sure RT is in minutes
    clusterinfo_filtered_df["retention time"] = clusterinfo_filtered_df["retention time"] / 60.0

    clusters_mean_df = clusterinfo_filtered_df.groupby("cluster index").mean()
    clusters_mean_df = clusters_mean_df.reset_index()
    # rename to RTMean
    clusters_mean_df = clusters_mean_df.rename(columns={"retention time": "RTMean"})

    # Merging this into the clustersummary
    cluster_summary_df = cluster_summary_df.merge(clusters_mean_df, on="cluster index", how="left")

    # Outputting again to be used downstream
    cluster_summary_df.to_csv(clustersummary_file, sep="\t", index=False)

    # Outputting clusterinfo file
    clusterinfo_df["#RetTime"] = clusterinfo_df["#RetTime"] / 60.0
    clusterinfo_df.to_csv(clusterinfo_file, sep="\t", index=False)

    # Do clean up out output spectra folder
    all_pklbin_files = glob.glob(os.path.join(args.output_spectra_folder, "specs_ms_*.pklbin"))

    """Disabling Removing Files because they are needed in a later step"""
    for filetoremove in all_pklbin_files:
       print("Removing ", filetoremove)
       os.remove(filetoremove)

if __name__ == "__main__":
    main()
