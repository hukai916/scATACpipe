#!/usr/bin/env python

"""
Split bed files into clusters given cluster file and fragment file.

Usage:
python split_bed.py Cluster_xxx.tsv Fragment.bed.gz

Cluster_xxx.tsv file must consist of two columns (first column is cell barcode, second is the assigned cluster), separated by tab.
Fragment.bed.gz must be in .bed or .bed.gz format, and the last second column must be cell barcode.
"""

import sys
import csv
import gzip
import os
import shutil
from pathlib import Path

cluster_file    = sys.argv[1]
fragment_file   = sys.argv[2]

cluster_dict    = {}

with open(cluster_file) as tsv_file:
    tsv = csv.reader(tsv_file, delimiter='\t')
    for row in tsv:
        cluster_dict[row[0]] = row[1]

clusters = set(x for x in cluster_dict.values())

# create output folder using cluser_file name:
out_dir = "split_" + Path(cluster_file).stem
if os.path.exists(out_dir):
    shutil.rmtree(out_dir)
os.mkdir(out_dir)

def split_bed(f, cluster_dict, output_cluster_dict):
    for i, line in enumerate(f):
        if i % 1000000 == 0:
            print(str(int(i / 1000000) * 1000000) + " reads processed ...")
        if not line.startswith("#"): # in case some software generated bed file starts with comment lines
            try: # some barcode are not in the cluster_dict
                barcode = line.split()[-2]
                if cluster_dict[barcode] in output_cluster_dict:
                    output_cluster_dict[cluster_dict[barcode]].append(line)
                else:
                    output_cluster_dict[cluster_dict[barcode]] = [line]
            except:
                continue

def save_file(out_dir, output_cluster_dict):
    print("Saving to output ...")
    for cluster in output_cluster_dict:
        output_file = out_dir + "/cluster_" + cluster + ".txt"
        with open(output_file, "w") as res:
            res.write("".join(output_cluster_dict[cluster]))

# loop over fragment file and output to corresponding cluster file.
output_cluster_dict = {}

if fragment_file.endswith(".gz"):
    with gzip.open(fragment_file, mode='rt') as f:
        split_bed(f, cluster_dict, output_cluster_dict)
    save_file(out_dir, output_cluster_dict)
elif fragment_file.endswith(".bed"):
    with open(fragment_file, mode='rt') as f:
        split_bed(f, cluster_dict, output_cluster_dict)
    save_file(out_dir, output_cluster_dict)
else:
    print("Pls supply either .bed or .bed.gz file.")
    exit()
print("Done!")
