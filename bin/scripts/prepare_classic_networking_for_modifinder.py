import sys
import argparse
import pandas as pd

from datetime import datetime

def getCurrentDatatimeAsString():
   
    # Get current date and time
    current_datetime = datetime.now()

    # Convert datetime object to string
    current_datetime_str = current_datetime.strftime('%Y-%m-%d-%H-%M-%S')
    
    return current_datetime_str

def dataframe_to_dict(df, key_column, value_column):
    # Initialize an empty dictionary
    result_dict = {}

    # Iterate over each row in the DataFrame
    for index, row in df.iterrows():
        # Retrieve key and value from specified columns
        key_value = row[key_column]
        value_value = row[value_column]
        
        # Add key-value pair to the dictionary
        result_dict[key_value] = value_value

    return result_dict
def scan2lib_row_dict(merged_df):
    scan2row_dict = {}
    for index, row in merged_df.iterrows():
        scan = int(row['#Scan#'])
        scan2row_dict[scan] = row

    return scan2row_dict

def process_pair(librow, clusterid2, gnps_taskid, addepairs, outputdata):
    clusterid1 = int(librow['#Scan#'])
    
    pair1 = f"{clusterid1}_{clusterid2}"
    pair2 = f"{clusterid2}_{clusterid1}"
    
    if pair1 not in addepairs:
    
        
        smiles = librow['Smiles']
        spectrumID = librow['SpectrumID']
        compoundName = librow['Compound_Name']

        usi1 = f"mzspec:GNPS:GNPS-LIBRARY:accession:{spectrumID}" 
        usi2 = f"mzspec:GNPS2:TASK-{gnps_taskid}-nf_output/clustering/specs_ms.mgf:scan:{clusterid2}"        

        outputdata.append({'USI1': usi1, 'USI2': usi2, 'SMILES1': smiles, 'SpectrumID': spectrumID, 'Compound_Name': compoundName, 'CLUSTERID1': clusterid1, 'CLUSTERID2': clusterid2 })

        addepairs.add(pair1)
        addepairs.add(pair2)

def process_bothIDed_pair(librow1, librow2, addepairs, outputdata):
    
    clusterid1 = int(librow1["#Scan#"])
    clusterid2 = int(librow2["#Scan#"])

    pair1 = f"{clusterid1}_{clusterid2}"
    pair2 = f"{clusterid2}_{clusterid1}"

    if pair1 not in addepairs:
    
        scan = int(librow1['#Scan#'])
        scan2 = int(librow2['#Scan#'])
        smiles = librow1['Smiles']
        smiles2 = librow2['Smiles']
        spectrumID = librow1['SpectrumID']
        spectrumID2 = librow2['SpectrumID']
        compoundName = librow1['Compound_Name']

        usi1 = f"mzspec:GNPS:GNPS-LIBRARY:accession:{spectrumID}"
        usi2 = f"mzspec:GNPS:GNPS-LIBRARY:accession:{spectrumID2}" 
        #usi2 = f"mzspec:GNPS2:TASK-{gnps_taskid}-nf_output/clustering/specs_ms.mgf:scan:{clusterid2}"        

        outputdata.append({'USI1': usi1, 'USI2': usi2, 'SMILES1': smiles, 'SpectrumID': spectrumID, 'Compound_Name': compoundName, 'CLUSTERID1': clusterid1, 'CLUSTERID2': clusterid2, 'SMILES2': smiles2})

        addepairs.add(pair1)
        addepairs.add(pair2)
        
def main():
    # Create the parser
    parser = argparse.ArgumentParser(description="Process some input files and a task ID.")

    # Add the arguments
    parser.add_argument('library_file', type=str, help='The merged results file with GNPS (TSV format)')
    parser.add_argument('filtered_pairs_file', type=str, help='The filtered pairs file (TSV format)')
    parser.add_argument('gnps2_taskid', type=str, help='The GNPS task ID')

    # output
    parser.add_argument('output_file', type=str, help='The output file (CSV format)')

    # Parse the arguments
    args = parser.parse_args()

    # Access the arguments
    merged_file = args.library_file
    filtered_pairs_file = args.filtered_pairs_file
    gnps2_taskid = args.gnps2_taskid
    
    # Load merged_results_with_gnps.tsv into a DataFrame
    try:
        merged_df = pd.read_csv(merged_file, sep='\t')
    except FileNotFoundError:
        print(f"Error: File '{merged_file}' not found.")
        sys.exit(1)
    
    #ignore hits with no smiles structure
    merged_df = merged_df.dropna(subset=['Smiles'])
    
    merged_df = merged_df[merged_df['Smiles'].str.strip() != '']

    # Load filtered_pairs.tsv into a DataFrame
    try:
        filtered_pairs_df = pd.read_csv(filtered_pairs_file, sep='\t')
    except FileNotFoundError:
        print(f"Error: File '{filtered_pairs_file}' not found.")
        sys.exit(1)
    
    #scan2smiles_dict = dataframe_to_dict(merged_df, '#Scan#', 'Smiles')
    scan2librowdict = scan2lib_row_dict(merged_df)

    # Prepare output data
    output_data = []
    
    addedPairs = set()
  
    # Check if scan is present in either CLUSTERID1 or CLUSTERID2 columns of filtered_pairs_df
    for index2, row2 in filtered_pairs_df.iterrows():
        suspect_id = 0
        clusterid1 = int(row2['CLUSTERID1'])
        clusterid2 = int(row2['CLUSTERID2'])
        librow = 0

        if clusterid1 in scan2librowdict and clusterid2 in scan2librowdict:
            #both scans are identified and have smiles
            librow1 = scan2librowdict[clusterid1]
            librow2 = scan2librowdict[clusterid2]
            # process_bothIDed_pair(librow1, librow2, addedPairs, output_data)

        elif clusterid1 in scan2librowdict:
            #clusterid1 is identified
            librow = scan2librowdict[clusterid1]
            process_pair(librow, clusterid2, gnps2_taskid, addedPairs, output_data)

        elif clusterid2 in scan2librowdict:
            #clusterid2 is identified
            librow = scan2librowdict[clusterid2]
            process_pair(librow, clusterid1, gnps2_taskid, addedPairs, output_data)
                
    if output_data:
        # Create output DataFrame
        output_df = pd.DataFrame(output_data)
        
        # Write output to file
        output_df.to_csv(args.output_file, sep=',', index=False)
        print(f"Output written to {args.output_file}")
        
    else:
        output_df = pd.DataFrame(columns=['USI1', 'USI2', 'SMILES1'])
        output_df.to_csv(args.output_file, sep=',', index=False)
        print("No matches found.")
    

if __name__ == "__main__":
    main()