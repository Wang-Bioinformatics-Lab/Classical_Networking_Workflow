#!/usr/bin/python

import sys
import os
import argparse
import pandas as pd


def _read_mzml_precursor_intensities(path):
    """(scan -> precursor intensity) for MS2 spectra in an mzML file.

    Reads pyteomics directly so a missing optional CV term on one spectrum
    (e.g. 'collision energy') does not abort parsing of the whole file.
    """
    from pyteomics import mzml as _pmzml
    out = {}
    for spectrum in _pmzml.read(path):
        if spectrum.get("ms level") != 2:
            continue
        scan = -1
        for tok in spectrum.get("id", "").split():
            if tok.startswith("scan="):
                try:
                    scan = int(tok.split("=", 1)[1])
                except ValueError:
                    pass
            elif tok.startswith("scanId="):
                try:
                    scan = int(tok.split("=", 1)[1])
                except ValueError:
                    pass
        if scan < 0:
            continue
        intensity = 0.0
        try:
            selected = (spectrum["precursorList"]["precursor"][0]
                                ["selectedIonList"]["selectedIon"][0])
            intensity = float(selected.get("peak intensity", 0.0))
        except (KeyError, IndexError, TypeError, ValueError):
            intensity = 0.0
        out[scan] = intensity
    return out


def _read_mzxml_precursor_intensities(path):
    """(scan -> precursor intensity) for MS2 spectra in an mzXML file."""
    from pyteomics import mzxml as _pmzxml
    out = {}
    for spectrum in _pmzxml.read(path):
        if int(spectrum.get("msLevel", 0)) != 2:
            continue
        try:
            scan = int(spectrum.get("num", -1))
        except (ValueError, TypeError):
            scan = -1
        if scan < 0:
            continue
        intensity = 0.0
        try:
            prec_list = spectrum.get("precursorMz") or []
            if prec_list:
                intensity = float(prec_list[0].get("precursorIntensity", 0.0))
        except (KeyError, IndexError, TypeError, ValueError):
            intensity = 0.0
        out[scan] = intensity
    return out


def build_precursor_intensity_lookup(input_spectra_folder, needed_basenames):
    """Map (basename, scan) -> precursor intensity by re-reading the input spectra.

    Falcon's CSV doesn't carry per-spectrum precursor intensity, so without
    this lookup #PrecIntensity stays 0 and every cell of the
    precursor-intensity feature quant table is 0. MSCluster's binary reads the
    value from the input mzML directly; we do the same here, extracting the
    'peak intensity' CV term of the selected precursor ion (or
    @precursorIntensity for mzXML). MGF inputs are skipped because standard
    MGF does not carry precursor intensity.
    """
    intensity_lookup = {}
    if not input_spectra_folder or not os.path.isdir(input_spectra_folder):
        print(f"WARNING: input_spectra_folder '{input_spectra_folder}' is not a directory; "
              f"precursor intensities will default to 0")
        return intensity_lookup

    for basename in needed_basenames:
        spectra_path = os.path.join(input_spectra_folder, basename)
        if not os.path.isfile(spectra_path):
            print(f"WARNING: input spectrum file not found: {spectra_path}; "
                  f"precursor intensities for this file will default to 0")
            continue

        ext = os.path.splitext(basename)[1].lower()
        try:
            if ext == ".mzml":
                per_file = _read_mzml_precursor_intensities(spectra_path)
            elif ext == ".mzxml":
                per_file = _read_mzxml_precursor_intensities(spectra_path)
            else:
                # MGF and unknown extensions: no per-spectrum precursor
                # intensity available, default to 0.
                continue
        except Exception as e:
            print(f"WARNING: failed to read precursor intensities from {spectra_path}: {e}")
            continue

        for scan_int, intensity in per_file.items():
            intensity_lookup[(basename, scan_int)] = float(intensity)

    return intensity_lookup


def convert_falcon_to_mscluster_format(falcon_csv, input_spectra_folder, output_clusterinfo, output_clustersummary, min_cluster_size=2):
    """
    Convert falcon output to mscluster format.

    Falcon format: cluster, filename, scan, precursor_mz, retention_time, new_batch
    MSCluster format: #ClusterIdx, #Filename, #SpecIdx, #Scan, #ParentMass, #Charge, #RetTime, #PrecIntensity

    Precursor intensity is not present in falcon's CSV, so it is looked up from
    the original input spectra files (mzML/mzXML) using ming_spectrum_library.
    Charge falls back to falcon's columns if available, otherwise defaults to 0.
    """
    # Load falcon CSV
    clusterinfo_df = pd.read_csv(falcon_csv, sep=',', comment='#')
    
    print(f"Loaded {len(clusterinfo_df)} rows from falcon CSV")
    print(f"Columns: {clusterinfo_df.columns.tolist()}")
    
    
    if 'identifier' in clusterinfo_df.columns:

        clusterinfo_df["filename"] = clusterinfo_df["identifier"].apply(
            lambda x: x.split(":")[2] + ".mzML" if len(x.split(":")) > 2 else (x.split(":")[-2] + ".mzML" if len(x.split(":")) > 1 else "unknown.mzML")
        )
        clusterinfo_df["scan"] = clusterinfo_df["identifier"].apply(
            lambda x: int(x.split(":")[-1]) if x.split(":")[-1].isdigit() else 0
        )
    
    # Ensure we have the required columns
    required_cols = ['cluster', 'filename', 'scan', 'precursor_mz', 'retention_time']
    missing_cols = [col for col in required_cols if col not in clusterinfo_df.columns]
    
    if missing_cols:
        print(f"ERROR: Required columns not found in falcon CSV: {missing_cols}")
        print(f"Available columns: {clusterinfo_df.columns.tolist()}")
        print(f"First few rows:")
        print(clusterinfo_df.head())
        sys.exit(1)
    
    if min_cluster_size > 1:
        clusterinfo_df = clusterinfo_df[clusterinfo_df['cluster'] != -1]

    # Build (basename, scan) -> precursor intensity lookup from the original
    # input spectra files. Falcon does not write precursor intensities into its
    # CSV, so without this every value in the precursor-intensity feature quant
    # table ends up as 0.
    needed_basenames = sorted({os.path.basename(str(fn)) for fn in clusterinfo_df['filename'].unique()})
    precursor_intensity_lookup = build_precursor_intensity_lookup(input_spectra_folder, needed_basenames)
    print(f"Loaded precursor intensities for {len(precursor_intensity_lookup)} spectra "
          f"across {len(needed_basenames)} input files")

    # Convert to mscluster format
    mscluster_rows = []
    spec_idx_counter = 0

    for idx, row in clusterinfo_df.iterrows():
        cluster_idx = int(row['cluster'])
        
        # Skip singletons if min_cluster_size > 1
        if cluster_idx == -1 and min_cluster_size > 1:
            continue
        
        # Handle cluster indexing: falcon uses 0-based, mscluster uses 1-based
        # But we also need to handle -1 (singletons)
        if cluster_idx == -1:
            # Singletons: use a large number or handle separately
            cluster_idx = 999999  # Use a large number for singletons
        else:
            cluster_idx = cluster_idx + 1  # Convert to 1-based
        
        filename = str(row['filename'])
        scan = int(row['scan'])
        precursor_mz = float(row['precursor_mz'])
        retention_time = float(row['retention_time'])
        
        # Use default values for fields that falcon doesn't provide
        # Don't fetch from MGF files - just use defaults
        charge = 0
        precursor_intensity = 0.0
        
   
        if 'precursor_charge' in row and pd.notna(row['precursor_charge']):
            try:
                charge = int(row['precursor_charge'])
            except (ValueError, TypeError):
                charge = 0
        elif 'charge' in row and pd.notna(row['charge']):
            try:
                charge = int(row['charge'])
            except (ValueError, TypeError):
                charge = 0
        
        if 'precursor_intensity' in row and pd.notna(row['precursor_intensity']):
            try:
                precursor_intensity = float(row['precursor_intensity'])
            except (ValueError, TypeError):
                precursor_intensity = 0.0
        else:
            precursor_intensity = precursor_intensity_lookup.get(
                (os.path.basename(filename), scan), 0.0
            )


        # Convert retention time to seconds if it's in minutes
        # Falcon typically outputs RT in minutes, mscluster uses seconds
        if retention_time > 0 and retention_time < 1000:  # Likely in minutes
            retention_time = retention_time * 60.0
        
        if not filename.startswith('input_spectra/'):
            filename = f"input_spectra/{filename}"
        
        spec_idx = spec_idx_counter
        spec_idx_counter += 1
        
        mscluster_row = {
            '#ClusterIdx': cluster_idx,
            '#Filename': filename,
            '#SpecIdx': spec_idx,
            '#Scan': scan,
            '#ParentMass': precursor_mz,
            '#Charge': charge,
            '#RetTime': retention_time,
            '#PrecIntensity': precursor_intensity
        }
        mscluster_rows.append(mscluster_row)
    
    # Create DataFrame
    mscluster_df = pd.DataFrame(mscluster_rows)
    
    # IMPORTANT: Before filtering, save the original falcon cluster ID (0-based) for each row
    # This will help us match falcon MGF clusters to clusterinfo clusters later
    # The original falcon cluster ID is stored in the 'cluster' column from falcon.csv
    # We need to track: original_falcon_cluster (0-based) -> sequential_index (1-based)
    
    # First, add a column to track original falcon cluster (0-based) before filtering
    # We'll need to reconstruct this from the cluster_idx we converted
    # cluster_idx was converted from falcon's 0-based to 1-based, so original = cluster_idx - 1
    mscluster_df['_original_falcon_cluster'] = mscluster_df['#ClusterIdx'] - 1
    
    # Filter by min_cluster_size
    if min_cluster_size > 1:
        # Count spectra per cluster
        cluster_counts = mscluster_df['#ClusterIdx'].value_counts()
        valid_clusters = cluster_counts[cluster_counts >= min_cluster_size].index
        mscluster_df = mscluster_df[mscluster_df['#ClusterIdx'].isin(valid_clusters)]
    
    # IMPORTANT: Remap cluster indices to sequential (1, 2, 3, ...) for consistency
    # This ensures clusterinfo.tsv, clustersummary.tsv, and MGF SCANS all use the same sequential indices
    # This is necessary for compatibility with ExecMolecularParallelPairs which uses index-based CLUSTERID
    original_clusters = sorted(mscluster_df['#ClusterIdx'].unique())
    cluster_remap = {orig: new for new, orig in enumerate(original_clusters, start=1)}
    
    # Create mapping: original falcon cluster (0-based) -> sequential index (1-based)
    # This mapping will be used in falcon_wrapper.py to match MGF clusters
    falcon_cluster_to_sequential = {}
    for orig_cluster_idx in original_clusters:
        original_falcon_cluster = orig_cluster_idx - 1  # Convert back to 0-based
        sequential_idx = cluster_remap[orig_cluster_idx]
        falcon_cluster_to_sequential[original_falcon_cluster] = sequential_idx
    
    # Remap #ClusterIdx in mscluster_df
    mscluster_df['#ClusterIdx'] = mscluster_df['#ClusterIdx'].map(cluster_remap)
    
    # Remove the temporary column
    mscluster_df = mscluster_df.drop(columns=['_original_falcon_cluster'])
    
    # Create cluster summary with remapped indices
    cluster_summary_rows = []
    for cluster_idx in sorted(mscluster_df['#ClusterIdx'].unique()):
        cluster_data = mscluster_df[mscluster_df['#ClusterIdx'] == cluster_idx]
        num_spectra = len(cluster_data)
        
        # Calculate mean RT (in minutes)
        mean_rt = cluster_data['#RetTime'].mean() / 60.0
        
        # Calculate parent mass (mean of #ParentMass)
        parent_mass = cluster_data['#ParentMass'].mean()
        
        # Calculate precursor mass (parent mass / charge, then take mean)
        # If charge is 0, use parent mass directly
        precursor_masses = []
        for idx, row in cluster_data.iterrows():
            if row['#Charge'] > 0:
                precursor_masses.append(row['#ParentMass'] / row['#Charge'])
            else:
                precursor_masses.append(row['#ParentMass'])
        precursor_mass = sum(precursor_masses) / len(precursor_masses) if precursor_masses else parent_mass
        
        # Calculate precursor charge (most common charge, or mean if no clear mode)
        charges = cluster_data['#Charge'].values
        charges_nonzero = charges[charges > 0]
        if len(charges_nonzero) > 0:
            # Use mode (most common charge) using pandas
            charge_series = pd.Series(charges_nonzero)
            mode_values = charge_series.mode()
            if len(mode_values) > 0:
                precursor_charge = int(mode_values.iloc[0])
            else:
                precursor_charge = int(charges_nonzero.mean())
        else:
            precursor_charge = 0
        
        # Calculate sum of precursor intensity
        sum_precursor_intensity = cluster_data['#PrecIntensity'].sum()
        
        cluster_summary_row = {
            'cluster index': cluster_idx,  # Already remapped to sequential
            'number of spectra': num_spectra,
            'parent mass': parent_mass,
            'precursor charge': precursor_charge,
            'precursor mass': precursor_mass,
            'sum(precursor intensity)': sum_precursor_intensity,
            'RTMean': mean_rt
        }
        cluster_summary_rows.append(cluster_summary_row)
    
    cluster_summary_df = pd.DataFrame(cluster_summary_rows)
    cluster_summary_df = cluster_summary_df.sort_values('cluster index')
    
    # Ensure cluster index is string type to match network graph node types
    # Network graph nodes from pairs file (CLUSTERID1/CLUSTERID2) are typically strings
    # This ensures type matching when add_clusterinfo_summary_to_graph checks "if cluster_index in G"
    cluster_summary_df['cluster index'] = cluster_summary_df['cluster index'].astype(str)
    
    # Save outputs
    mscluster_df.to_csv(output_clusterinfo, sep='\t', index=False)
    cluster_summary_df.to_csv(output_clustersummary, sep='\t', index=False)
    
    print(f"Converted {len(mscluster_df)} spectra in {len(cluster_summary_df)} clusters")
    print(f"Saved clusterinfo to {output_clusterinfo}")
    print(f"Saved clustersummary to {output_clustersummary}")
    print(f"Note: Cluster indices have been remapped to sequential (1, 2, 3, ...) for consistency")


def main():
    parser = argparse.ArgumentParser(description='Convert Falcon output to MSCluster format')
    parser.add_argument('falcon_csv', help='Falcon CSV output file')
    parser.add_argument('input_spectra_folder', help='Input spectra folder (for reference, not used for fetching data)')
    parser.add_argument('output_clusterinfo', help='Output clusterinfo.tsv file')
    parser.add_argument('output_clustersummary', help='Output clustersummary.tsv file')
    parser.add_argument('--min_cluster_size', type=int, default=2, help='Minimum cluster size')
    
    args = parser.parse_args()
    
    convert_falcon_to_mscluster_format(
        args.falcon_csv,
        args.input_spectra_folder,
        args.output_clusterinfo,
        args.output_clustersummary,
        args.min_cluster_size
    )


if __name__ == "__main__":
    main()
