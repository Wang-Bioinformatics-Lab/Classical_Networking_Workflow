import argparse
import glob
import os
import csv

# Function to find files with specified extensions
def find_spectra_files(directory):
    spectra_files = []
    # Append files with specified extensions to the spectra_files list
    spectra_files += glob.glob(os.path.join(directory, "*.mzML"))
    spectra_files += glob.glob(os.path.join(directory, "*.mzml"))
    spectra_files += glob.glob(os.path.join(directory, "*.mzXML"))
    spectra_files += glob.glob(os.path.join(directory, "*.mzxml"))
    spectra_files += glob.glob(os.path.join(directory, "*.mgf"))
    spectra_files += glob.glob(os.path.join(directory, "*.MGF"))
    return spectra_files

# Main function to handle argument parsing and file output
def main():
    # Setup argument parser
    parser = argparse.ArgumentParser(description="Find all spectra files in a directory and save filenames to a TSV file.")
    parser.add_argument("directory", help="Directory to search for spectra files")
    parser.add_argument("output", help="Output TSV file to save the filenames")
    args = parser.parse_args()
    
    # Find spectra files in the specified directory
    spectra_files = find_spectra_files(args.directory)
    
    # Write filenames to the output TSV file
    with open(args.output, "w", newline='') as tsv_file:
        writer = csv.writer(tsv_file, delimiter='\t')
        writer.writerow(["Filename"])  # Header
        for file in spectra_files:
            writer.writerow([file])  # Write each filename in a new row

# Ensure the script runs only when executed directly
if __name__ == "__main__":
    main()
