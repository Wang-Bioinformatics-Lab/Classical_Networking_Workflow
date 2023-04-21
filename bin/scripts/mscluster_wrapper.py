#!/usr/bin/python

import sys
import os
import glob
import argparse
import ming_spectrum_library

def write_params(input_spectra_folder, tool_path, params_filename):

    # listing all files
    all_mgf_files = glob.glob(os.path.join(input_spectra_folder, "*.mgf"))
    all_mzxml_files = glob.glob(os.path.join(input_spectra_folder, "*.mzXML"))
    all_mzml_files = glob.glob(os.path.join(input_spectra_folder, "*.mzML"))

    all_spectrum_files = all_mgf_files + all_mzxml_files + all_mzml_files

    with open(params_filename, "w") as params_file:
        params_file.write("CLUSTER_MIN_SIZE=2\n")
        params_file.write("EXE_DIR={}\n".format(tool_path))
        params_file.write("INPUT_SPECS_MS={}\n".format(";".join(all_spectrum_files)))

    
def main():
    # Parse arguments
    parser = argparse.ArgumentParser(description='MSCluster Wrapper')
    parser.add_argument('input_spectra_folder', help='Input Spectra Folder')
    parser.add_argument('tool_dir', help='tool_dir')
    parser.add_argument('output_spectra_folder', help='Output Spectra Folder')
    
    args = parser.parse_args()

    # Writing the parameters
    parameters_filename = "mscluster.params"
    write_params(args.input_spectra_folder, args.tool_dir, parameters_filename)

    # Running the data
    mainspecnets_binary = os.path.join(args.tool_dir, "main_specnets")

    cmd = "{} {} -ll 9 -f mscluster".format(mainspecnets_binary, parameters_filename)
    ret_code = os.system(cmd)

    if ret_code != 0:
        exit(1)

    # Perform a rewrite of the mgf file, this does not need to be done, because original ordering is fine
    specs_mgf_filename = os.path.join(args.output_spectra_folder, "specs_ms.mgf")
    spectrum_collection = ming_spectrum_library.SpectrumCollection(specs_mgf_filename)
    spectrum_collection.load_from_mgf()

    spectrum_collection.save_to_mgf(open(specs_mgf_filename, "w"), renumber_scans=False)

    # Do clean up out output spectra folder
    # all_pklbin_files = glob.glob(os.path.join(output_spectra_folder, "specs_ms_*.pklbin"))



    """Disabling Removing Files because they are needed in a later step"""
    #TODO: Move clusterinfo into this step so we can get rid of these files. 
    #for filetoremove in all_pklbin_files:
    #    print("Removing ", filetoremove)
    #    os.remove(filetoremove)

if __name__ == "__main__":
    main()
