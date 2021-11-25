#!/usr/bin/env python

"""
Given whitelist_barcode folder, index_fastq, outfile_prefix, select the correct whitelist with the highest overlapping fraction. If fraction less than 0.5, exits with error msg.
"""

import sys
import gzip
import os
from shutil import copy2

whitelist_barcode_folder = sys.argv[1]
index_fastq    = sys.argv[2]
outfile_prefix = sys.argv[3]

def get_fraction(index_fastq, barcode_file):
    """
    Given index_fastq and barcode_file, calculate the percentage of index_fastq that overlaps with barcode.
    """

    whitelist_barcodes = {}
    if barcode_file.endswith(".gz"):
        with gzip.open(barcode_file, "rt") as f:
            for line in f:
                whitelist_barcodes[line.strip()] = 1
    else:
        with open(barcode_file) as f:
            for line in f:
                whitelist_barcodes[line.strip()] = 1

    reads_in, reads_not_in = 0, 0
    if index_fastq.endswith(".gz"):
        with gzip.open(index_fastq, "rt") as f:
            for i, line in enumerate(f):
                if i%4 == 1:
                    if line.strip() in whitelist_barcodes:
                        reads_in += 1
                    else:
                        reads_not_in += 1
    else:
        with open(index_fastq) as f:
            for i, line in enumerate(f):
                if i%4 == 1:
                    if line.strip() in whitelist_barcodes:
                        reads_in += 1
                    else:
                        reads_not_in += 1

    if (reads_in + reads_not_in) == 0:
        sys.exit("Error: no reads detected in index fastq file!")
    else:
        return reads_in / (reads_in + reads_not_in)

barcode_files = [os.path.join(whitelist_barcode_folder, file) for file in os.listdir(whitelist_barcode_folder)]
fractions     = [get_fraction(index_fastq, barcode_file) for barcode_file in barcode_files]
max_frac      = fractions.index(max(fractions))

if max(fractions) < 0.5:
    sys.exit("Valid barcode ratio is less than 0.5 across all supplied whitelist barcode files, must used the wrong whitelist barcode file!")
else:
    copy2(barcode_files[max_frac], os.path.join("./" + outfile_prefix + "_" + os.path.basename(barcode_files[max_frac])))
    print("Selected whitelist barcode file: " + barcode_files[max_frac])
    print("Valid barcode fraction against this whitelist: " + str(max(fractions)))
