import argparse
import shutil

def main():
    parser = argparse.ArgumentParser(description='Filtering of spectra')
    parser.add_argument('input_mgf', help='Input file')
    parser.add_argument('output_mgf', help='Output file')
    parser.add_argument('massql', help='Query used for filtering')

    args = parser.parse_args()

    if args.massql == "None":
        # copy the inptu to the output and call it a day
        shutil.copyfile(args.input_mgf, args.output_mgf)
        exit(0)

    # Here we'll actually do the filtering

    # Lets parse via the API
    url = "https://massql.gnps2.org/parse?query=".format(args.massql)


    # This is debug no-op
    shutil.copyfile(args.input_mgf, args.output_mgf)


if __name__ == '__main__':
    main()