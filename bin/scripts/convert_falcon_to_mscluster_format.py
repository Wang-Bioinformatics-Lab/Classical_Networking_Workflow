#!/usr/bin/python

import sys
import os
import argparse
import pandas as pd



def convert_falcon_to_mscluster_format(falcon_csv, input_spectra_folder, output_clusterinfo, output_clustersummary, min_cluster_size=2):
    """
    Convert falcon output to mscluster format.
    
    Falcon format: cluster, filename, scan, precursor_mz, retention_time, new_batch
    MSCluster format: #ClusterIdx, #Filename, #SpecIdx, #Scan, #ParentMass, #Charge, #RetTime, #PrecIntensity
    
    Note: If falcon doesn't have certain fields (like charge, intensity), we use default values (0)
    instead of fetching from MGF files.
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
    
    # Filter by min_cluster_size
    if min_cluster_size > 1:
        # Count spectra per cluster
        cluster_counts = mscluster_df['#ClusterIdx'].value_counts()
        valid_clusters = cluster_counts[cluster_counts >= min_cluster_size].index
        mscluster_df = mscluster_df[mscluster_df['#ClusterIdx'].isin(valid_clusters)]
    
    # Create cluster summary
    cluster_summary_rows = []
    for cluster_idx in mscluster_df['#ClusterIdx'].unique():
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
            'cluster index': cluster_idx,
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
    
    # Save outputs
    mscluster_df.to_csv(output_clusterinfo, sep='\t', index=False)
    cluster_summary_df.to_csv(output_clustersummary, sep='\t', index=False)
    
    print(f"Converted {len(mscluster_df)} spectra in {len(cluster_summary_df)} clusters")
    print(f"Saved clusterinfo to {output_clusterinfo}")
    print(f"Saved clustersummary to {output_clustersummary}")


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
