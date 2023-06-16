#!/usr/bin/python
import sys
import os
import json
import argparse
from collections import defaultdict
import csv
import pandas as pd
import glob
import shutil
import requests
import yaml

def main():
    parser = argparse.ArgumentParser(description='Running library search parallel')
    parser.add_argument('input_download_file', help='input_download_file')
    parser.add_argument('output_folder', help='output_folder')
    parser.add_argument('--existing_data', default=None, help='folder of existing data')
    args = parser.parse_args()

    print(args)

    # checking the input file exists
    if not os.path.isfile(args.input_download_file):
        print("Input file does not exist")
        exit(0)

    # Checking the file extension
    if args.input_download_file.endswith(".yaml"):
        # Loading yaml file
        parameters = yaml.load(open(args.input_download_file), Loader=yaml.SafeLoader)
        usi_list = parameters["usi"].split("\n")
    elif args.input_download_file.endswith(".tsv"):
        df = pd.read_csv(args.input_download_file, sep="\t")
        usi_list = df["usi"].tolist()

    # Lets download these files
    for usi in usi_list:
        # Getting the path to the original file
        url = "https://dashboard.gnps2.org/downloadlink"
        params = {"usi": usi}
        r = requests.get(url, params=params)

        if r.status_code == 200:
            download_url = r.text

            # download in chunks using requests
            r = requests.get(download_url, stream=True)
            with open(os.path.join(args.output_folder, os.path.basename(download_url)), 'wb') as fd:
                for chunk in r.iter_content(chunk_size=128):
                    fd.write(chunk)

    # copying all files from existing data into the output folder
    if args.existing_data is not None:
        src = args.existing_data
        dst = os.path.join(args.output_folder, os.path.basename(args.existing_data))
        if os.path.islink(src):
            linkto = os.readlink(src)
            os.symlink(linkto, dst)
        else:
            shutil.copytree(src, dst)

        
        



if __name__ == "__main__":
    main()
