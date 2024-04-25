import os
import pandas as pd
import argparse

def create_feature_table(clusterinfo_df, type="presence"):
    # applying basename
    clusterinfo_df['#Filename'] = clusterinfo_df['#Filename'].apply(lambda x: os.path.basename(x))
    all_filenames = clusterinfo_df['#Filename'].unique()

    print(all_filenames)

    # grouping by filename
    grouped_df = clusterinfo_df.groupby('#ClusterIdx')

    # iterate through groups
    feature_table = []
    for name, group in grouped_df:
        feature_dict = {}

        # by default lets calculate the mean of the m/z and retention time
        mean_mz = group['#ParentMass'].mean()
        mean_rt = group['#RetTime'].mean()

        feature_dict["row ID"] = name
        feature_dict["row m/z"] = mean_mz
        feature_dict["row retention time"] = mean_rt

        for filename in all_filenames:
            # checking if these filenames are in here
            if type == "presence":
                if filename in group['#Filename'].values:
                    feature_dict[filename + " Peak area"] = 1
                else:
                    feature_dict[filename + " Peak area"] = 0
            elif type == "spectrumcount":
                if filename in group['#Filename'].values:
                    # count the number of times it appears
                    feature_dict[filename + " Peak area"] = len(group[group['#Filename'] == filename])
                else:
                    feature_dict[filename + " Peak area"] = 0
            elif type == "precursorintensity":
                if filename in group['#Filename'].values:
                    # count the number of times it appears
                    feature_dict[filename + " Peak area"] = max(group[group['#Filename'] == filename]["#PrecIntensity"])
                else:
                    feature_dict[filename + " Peak area"] = 0

        feature_table.append(feature_dict)

    # create a dataframe
    feature_table_df = pd.DataFrame(feature_table)

    return feature_table_df

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('input_clusterinfo')
    parser.add_argument('output_featuretable_presence')
    parser.add_argument('output_featuretable_spectrumcount')
    parser.add_argument('output_featuretable_precursorintensity')
    args = parser.parse_args()

    clusterinfo_df = pd.read_csv(args.input_clusterinfo, sep='\t')

    print(clusterinfo_df.head())

    # presence abscence
    feature_table_df = create_feature_table(clusterinfo_df, type="presence")
    feature_table_df.to_csv(args.output_featuretable_presence, sep=',', index=False)

    # precursor intensity
    feature_table_df = create_feature_table(clusterinfo_df, type="spectrumcount")
    feature_table_df.to_csv(args.output_featuretable_spectrumcount, sep=',', index=False)

    # precursor intensity
    feature_table_df = create_feature_table(clusterinfo_df, type="precursorintensity")
    feature_table_df.to_csv(args.output_featuretable_precursorintensity, sep=',', index=False)






    

    

if __name__ == '__main__':
    main()