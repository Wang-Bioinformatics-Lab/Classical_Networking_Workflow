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
import uuid

def main():
    parser = argparse.ArgumentParser(description='Running library search parallel')
    parser.add_argument('input_download_file', help='input_download_file')
    parser.add_argument('output_folder', help='output_folder')
    parser.add_argument('--cache_directory', default=None, help='folder of existing data')
    args = parser.parse_args()

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

            target_filename = os.path.basename(download_url)
            target_path = os.path.join(args.output_folder, target_filename)

            # Checking the cache
            if args.cache_directory is not None and os.path.exists(args.cache_directory):
                
                namespace = uuid.UUID('6ba7b810-9dad-11d1-80b4-00c04fd430c8')
                hashed_id = str(uuid.uuid3(namespace, usi)).replace("-", "")

                cache_path = os.path.join(args.cache_directory, hashed_id)
                cache_path = os.path.realpath(cache_path)
                
                cache_filename = cache_path + "-" + target_filename[-50:]

                # If we find it, we can create a link to it
                if os.path.exists(cache_filename):
                    print("Found in cache", cache_path)
                    os.symlink(cache_filename, target_path)

                    continue

                # Saving file to cache if we don't
                r = requests.get(download_url, stream=True)
                try:
                    with open(os.path.join(args.output_folder, cache_filename), 'wb') as fd:
                        for chunk in r.iter_content(chunk_size=128):
                            fd.write(chunk)
                        
                    # Creating symlink
                    os.symlink(cache_filename, target_path)
                except:
                    # We are likely writing to read only file system for the cache
                    with open(target_path, 'wb') as fd:
                        for chunk in r.iter_content(chunk_size=128):
                            fd.write(chunk)
            else:
                # download in chunks using requests
                r = requests.get(download_url, stream=True)
                with open(target_path, 'wb') as fd:
                    for chunk in r.iter_content(chunk_size=128):
                        fd.write(chunk)

if __name__ == "__main__":
    main()
