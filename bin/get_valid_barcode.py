#!/usr/bin/env python

"""
Given valid_barcode_frequency.txt, whitelist_barcode.txt, outfile name, output valid barcode by overlapping whitelist and valid_barcode.
"""

import sys
import gzip
import os
from shutil import copy2

valid_barcode_frequency = sys.argv[1]
whitelist_barcode       = sys.argv[2]
outfile_name            = sys.argv[3]
outfile_frequency       = sys.argv[4]

dict_valid_barcode = {}
if valid_barcode_frequency.endswith(".gz"):
    with gzip.open(valid_barcode_frequency, "rt") as f:
        for line in f:
            tem = line.split()
            if not tem[0] in dict_valid_barcode:
                dict_valid_barcode[tem[0]] = tem[1]
else:
    with open(valid_barcode_frequency) as f:
        for line in f:
            tem = line.split()
            if not tem[0] in dict_valid_barcode:
                dict_valid_barcode[tem[0]] = tem[1]

# output valid_barcode:
with open(outfile_name, "w") as fout:
    if whitelist_barcode.endswith(".gz"):
        with gzip.open(whitelist_barcode, "rt") as fin:
            for line in fin:
                tem = line.split()[0]
                if tem in dict_valid_barcode:
                    fout.write(tem + "\n")
    else:
        with open(whitelist_barcode) as fin:
            for line in fin:
                tem = line.split()[0]
                if tem in dict_valid_barcode:
                    fout.write(tem + "\n")

# output valid_barcode_frequency:
with open(outfile_frequency, "w") as fout:
    if valid_barcode_frequency.endswith(".gz"):
        with gzip.open(valid_barcode_frequency, "rt") as fin:
            for line in fin:
                tem = line.split()[0]
                if tem in dict_valid_barcode:
                    fout.write(line, end = "")
    else:
        with open(valid_barcode_frequency) as fin:
            for line in fin:
                tem = line.split()[0]
                if tem in dict_valid_barcode:
                    fout.write(line, end = "")
