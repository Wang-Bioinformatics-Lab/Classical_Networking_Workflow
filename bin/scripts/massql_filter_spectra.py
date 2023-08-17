import argparse
import shutil

from massql import msql_parser
from massql import msql_engine

def main():
    parser = argparse.ArgumentParser(description='Filtering of spectra')
    parser.add_argument('input_mgf', help='Input file')
    parser.add_argument('output_mgf', help='Output file')
    parser.add_argument('--massql', help='Query used for filtering')

    args = parser.parse_args()

    if args.massql == "None":
        # copy the inptu to the output and call it a day
        shutil.copyfile(args.input_mgf, args.output_mgf)
        exit(0)

    # Here we'll actually do the filtering

    # Lets parse via the API
    parsed_version = msql_parser.parse_msql(args.massql)

    results_df = msql_engine.process_query(args.massql, args.input_mgf)

    print(results_df)

    # now what we can do here is get the full ms2 data, and then filter only the ones in the original
    ms1_df, ms2_df = msql_fileloading.load_data(args.input_mgf)

    # getting values unique to the ms2_df
    unique_to_original = ms2_df.merge(results_df, how='left', indicator=True).query('_merge == "left_only"').drop('_merge', axis=1)

    raise Exception()



    # This is debug no-op
    shutil.copyfile(args.input_mgf, args.output_mgf)


if __name__ == '__main__':
    main()