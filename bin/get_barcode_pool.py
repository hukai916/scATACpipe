#!/usr/bin/env python

"""
Given full whitelist and index fastq, determine the valid barcode pool and its count.

Usage:
python get_barcode_pool.py whitelist.txt.gz index.fastq.gz read_count_cutoff outfile.txt

read_count_cutoff is an integer indicating the minimum number of reads to include certain barcode into the valid barcode pool.
"""

import sys
import gzip
import pysam

whitelist_file = sys.argv[1]
index_fastq    = sys.argv[2]
read_count_cutoff   = int(sys.argv[3])
outfile_name   = sys.argv[4]

# read in whitelist and make a dict:
whitelist_dict = {}
if whitelist_file.endswith(".gz"):
    with gzip.open(whitelist_file, "rt") as f:
        for line in f:
            if not line.strip() in whitelist_dict:
                whitelist_dict[line.strip()] = 0
else:
    with open(whitelist_file) as f:
        for line in f:
            if not line.strip() in whitelist_dict:
                whitelist_dict[line.strip()] = 0

# count barcode from index fastq file:
not_in_whitelist = 0 # count how many barcodes are not in whitelist
with pysam.FastxFile(index_fastq) as f:
    for read in f:
        if read.sequence in whitelist_dict:
            whitelist_dict[read.sequence] += 1
        else:
            not_in_whitelist += 1

total_valid_barcode = sum([x for x in list(whitelist_dict.values()) if x > read_count_cutoff])
total_unique_valid_barcode = sum([1 for x in list(whitelist_dict.values()) if x > read_count_cutoff])

# store output:
if outfile_name.endswith(".gz"):
    out_file = gzip.open(outfile_name, "wt")
else:
    out_file = open(outfile_name, "wt")

for barcode in whitelist_dict:
    if whitelist_dict[barcode] > read_count_cutoff:
        # out_file.write("\t".join([barcode, str(whitelist_dict[barcode] / total_valid_barcode)]) + "\n") # for frequency
        out_file.write("\t".join([barcode, str(whitelist_dict[barcode])]) + "\n") # for raw counts
out_file.close()

# output some statistics:
print("total valid bacodes: ", total_valid_barcode)
print("total unique valid barcodes: ", total_unique_valid_barcode)
print("not_in_whitelist / total valid barcodes: ", not_in_whitelist/total_valid_barcode)
print("If above value > 0.25, you may have used the wrong whitelist.")
print("Done!")
