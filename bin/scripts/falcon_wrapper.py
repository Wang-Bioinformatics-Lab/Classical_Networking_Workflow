#!/usr/bin/python

import sys
import os
import glob
import argparse
import subprocess
import pandas as pd
import ming_spectrum_library

def run_falcon(input_spectra_path, output_prefix="falcon", 
               precursor_tol="20 ppm", fragment_tol=0.05,
               min_mz_range=0, min_mz=0, max_mz=30000, eps=0.1):
    """
    Runs Falcon with specified parameters.
    input_spectra_path can be a file or a folder.
    """
    # Check if input is a file or folder
    if os.path.isfile(input_spectra_path):
        # Single file input
        all_spectrum_files = [input_spectra_path]
    elif os.path.isdir(input_spectra_path):
        # Folder input - list all spectrum files
        all_mgf_files = glob.glob(os.path.join(input_spectra_path, "*.mgf"))
        all_mzxml_files = glob.glob(os.path.join(input_spectra_path, "*.mzXML"))
        all_mzml_files = glob.glob(os.path.join(input_spectra_path, "*.mzML"))
        
        all_spectrum_files = all_mgf_files + all_mzxml_files + all_mzml_files
        
        if len(all_spectrum_files) == 0:
            print(f"ERROR: No spectrum files found in {input_spectra_path}")
            sys.exit(1)
        
        # Sort these filenames
        all_spectrum_files.sort()
    else:
        print(f"ERROR: Input path does not exist: {input_spectra_path}")
        sys.exit(1)
    
    # Create pattern for falcon (it accepts wildcards or file list)
    if len(all_spectrum_files) == 1:
        mzml_pattern = all_spectrum_files[0]
    else:
        # Falcon can accept multiple files, we'll pass them as space-separated
        mzml_pattern = " ".join(all_spectrum_files)
    
    # Parse precursor_tol - falcon expects two arguments
    # Format can be "2.0" (just number) or "20 ppm" (value and unit)
    # Falcon needs two arguments, so we need to handle both cases
    precursor_tol_str = str(precursor_tol).strip()
    if " " in precursor_tol_str:
        # Already has format like "20 ppm" or "2.0 Da" - split into two args
        parts = precursor_tol_str.split(None, 1)  # Split on first space
        if len(parts) == 2:
            precursor_tol_args = f"{parts[0]} {parts[1]}"
        else:
            # Fallback: use value twice
            precursor_tol_args = f"{parts[0]} {parts[0]}"
    else:
        # Just a number, assume Da and use twice (min and max tolerance)
        try:
            tol_value = float(precursor_tol_str)
            precursor_tol_args = f"{tol_value} {tol_value}"
        except ValueError:
            # If can't parse, use default
            precursor_tol_args = "2.0 2.0"
    
    command = (
        f"falcon {mzml_pattern} {output_prefix} "
        f"--export_representatives "
        f"--precursor_tol {precursor_tol_args} "
        f"--fragment_tol {fragment_tol} "
        f"--min_mz_range {min_mz_range} "
        f"--min_mz {min_mz} --max_mz {max_mz} "
        f"--eps {eps} "
        f"--hash_len 400 "
        f"--n_neighbors_ann 64 "
        f"--n_probe 16 "
        f"--batch_size 32768 "
    )
    print(f"[run_falcon] Running: {command}")
    process = subprocess.Popen(command, shell=True)
    retcode = process.wait()
    if retcode != 0:
        print(f"ERROR: Falcon failed with exit code {retcode}")
        sys.exit(1)


def main():
    # Parse arguments
    parser = argparse.ArgumentParser(description='Falcon Clustering Wrapper')
    parser.add_argument('input_spectra_folder', help='Input Spectra Folder')
    parser.add_argument('output_spectra_folder', help='Output Spectra Folder')
    parser.add_argument('final_output_folder', help='final_output_folder')
    
    parser.add_argument('--min_cluster_size', default="2", help='min_cluster_size (not used by falcon directly)')
    
    parser.add_argument('--pm_tolerance', default="20 ppm", help='pm_tolerance (precursor tolerance)')
    parser.add_argument('--fragment_tolerance', default="0.05", help='fragment_tolerance')
    
    # Filters (not all are used by falcon)
    parser.add_argument('--min_peak_intensity', default="0.0", help='min_peak_intensity (not used by falcon)')
    parser.add_argument('--window_filter', default="1", help='window_filter (not used by falcon)')
    parser.add_argument('--precursor_filter', default="1", help='precursor_filter (not used by falcon)')
    
    # Falcon-specific parameters
    parser.add_argument('--eps', default="0.1", help='Falcon eps parameter')
    parser.add_argument('--min_mz', default="0", help='Falcon min_mz parameter')
    parser.add_argument('--max_mz', default="30000", help='Falcon max_mz parameter')
    
    args = parser.parse_args()
    
    # Check if input is a file or folder (Nextflow may pass a file)
    input_spectra_path = args.input_spectra_folder
    if os.path.isfile(args.input_spectra_folder):
        # If it's a single file, we need to use its directory
        input_spectra_path = os.path.dirname(args.input_spectra_folder)
        if not input_spectra_path:
            input_spectra_path = "."
    
    # Running falcon
    output_prefix = "falcon"
    run_falcon(input_spectra_path, output_prefix,
               precursor_tol=args.pm_tolerance,
               fragment_tol=float(args.fragment_tolerance),
               min_mz_range=0,
               min_mz=int(args.min_mz),
               max_mz=int(args.max_mz),
               eps=float(args.eps))
    
    # Falcon outputs:
    # - falcon.csv (cluster assignments)
    # - falcon.mgf (representative spectra)
    
    falcon_csv = f"{output_prefix}.csv"
    falcon_mgf = f"{output_prefix}.mgf"
    
    if not os.path.exists(falcon_csv):
        print(f"ERROR: Falcon output {falcon_csv} not found")
        sys.exit(1)
    
    if not os.path.exists(falcon_mgf):
        print(f"ERROR: Falcon output {falcon_mgf} not found")
        sys.exit(1)
    
    # Convert falcon output to mscluster format first
    script_dir = os.path.dirname(os.path.abspath(__file__))
    convert_script = os.path.join(script_dir, "convert_falcon_to_mscluster_format.py")
    
    clusterinfo_file = os.path.join(args.final_output_folder, "clusterinfo.tsv")
    clustersummary_file = os.path.join(args.final_output_folder, "clustersummary.tsv")
    
    # For conversion, we need the original input path to get spectrum info
    convert_input_path = input_spectra_path
    
    convert_cmd = (
        f"python {convert_script} "
        f"{falcon_csv} "
        f"{convert_input_path} "
        f"{clusterinfo_file} "
        f"{clustersummary_file} "
        f"--min_cluster_size {args.min_cluster_size}"
    )
    
    print(f"Converting falcon output to mscluster format...")
    print(f"Running: {convert_cmd}")
    ret_code = os.system(convert_cmd)
    
    if ret_code != 0:
        print(f"ERROR: Conversion failed with exit code {ret_code}")
        sys.exit(1)
    
    # Now update MGF file with correct scan numbers based on clusterinfo.tsv
    # The scan number in MGF should match the #ClusterIdx in clusterinfo.tsv
    # Read clusterinfo to get cluster indices
    clusterinfo_df = pd.read_csv(clusterinfo_file, sep='\t')
    
    # Get unique cluster indices (these will be the scan numbers in MGF)
    unique_clusters = sorted(clusterinfo_df['#ClusterIdx'].unique())
    
    # Read falcon MGF file - parse manually to access cluster field
    output_mgf_filename = os.path.join(args.final_output_folder, "specs_ms.mgf")
    os.makedirs(args.final_output_folder, exist_ok=True)
    
    # Parse falcon MGF to get cluster information and update scan numbers
    mgf_spectra = []
    with open(falcon_mgf, 'r') as f:
        spec = None
        for line in f:
            line = line.strip()
            if line == 'BEGIN IONS':
                spec = {'peaks': [], 'headers': {}}
            elif line == 'END IONS':
                if spec is not None:
                    mgf_spectra.append(spec)
                spec = None
            elif spec is not None:
                if '=' in line:
                    key, val = line.split('=', 1)
                    spec['headers'][key.upper()] = val
                else:
                    # Peak data
                    parts = line.split()
                    if len(parts) == 2:
                        try:
                            mz, intensity = float(parts[0]), float(parts[1])
                            spec['peaks'].append((mz, intensity))
                        except:
                            pass
    
    # Create mapping: falcon cluster (0-based) -> mscluster cluster index (1-based, which is #ClusterIdx)
    # falcon cluster IDs in MGF are 0-based, but we converted to 1-based in clusterinfo
    cluster_to_scan = {}
    for cluster_idx in unique_clusters:
        # falcon uses 0-based, mscluster uses 1-based
        falcon_cluster = cluster_idx - 1
        cluster_to_scan[falcon_cluster] = cluster_idx
    
    # Write MGF with correct scan numbers
    # Only include spectra that correspond to clusters in clusterinfo.tsv
    # This ensures MGF and clusterinfo.tsv are consistent
    skipped_count = 0
    with open(output_mgf_filename, 'w') as out_mgf:
        for spec in mgf_spectra:
            # Get cluster ID from falcon MGF
            falcon_cluster = -1
            if 'CLUSTER' in spec['headers']:
                try:
                    falcon_cluster = int(spec['headers']['CLUSTER'])
                except (ValueError, TypeError):
                    pass
            
            # Only include spectra that are in clusterinfo.tsv
            # Skip clusters that were filtered out by min_cluster_size
            if falcon_cluster in cluster_to_scan:
                scan_number = cluster_to_scan[falcon_cluster]
                
                # Write MGF entry with correct scan number
                out_mgf.write("BEGIN IONS\n")
                
                # Write headers, updating SCANS= line
                for key, val in spec['headers'].items():
                    if key == 'SCANS':
                        out_mgf.write(f"SCANS={scan_number}\n")
                    else:
                        out_mgf.write(f"{key}={val}\n")
                
                # Add SCANS if it wasn't in the original
                if 'SCANS' not in spec['headers']:
                    out_mgf.write(f"SCANS={scan_number}\n")
                
                # Write peaks
                for mz, intensity in spec['peaks']:
                    out_mgf.write(f"{mz} {intensity}\n")
                
                out_mgf.write("END IONS\n\n")
            else:
                # Skip clusters that are not in clusterinfo.tsv (filtered out by min_cluster_size)
                skipped_count += 1
                if falcon_cluster != -1:  # Only warn for non-singleton clusters
                    print(f"INFO: Skipping falcon cluster {falcon_cluster} (filtered out by min_cluster_size)")
    
    if skipped_count > 0:
        print(f"INFO: Skipped {skipped_count} spectra that were filtered out by min_cluster_size")
    
    print(f"Falcon clustering completed. Output files:")
    print(f"  - MGF: {output_mgf_filename}")
    print(f"  - ClusterInfo: {clusterinfo_file}")
    print(f"  - ClusterSummary: {clustersummary_file}")


if __name__ == "__main__":
    main()
