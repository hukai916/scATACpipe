#!/usr/bin/env python

"""
Filter fragment file based on valid barcode list.
"""

import sys
import gzip
import os

fragment        = sys.argv[1]
valid_barcode   = sys.argv[2]

dict_valid_barcode = {}
if valid_barcode.endswith(".gz"):
    with gzip.open(valid_barcode, "rt") as f:
        for line in f:
            tem = line.split()
            if not tem[0] in dict_valid_barcode:
                dict_valid_barcode[tem[0]] = 1
else:
    with open(valid_barcode) as f:
        for line in f:
            tem = line.split()
            if not tem[0] in dict_valid_barcode:
                dict_valid_barcode[tem[0]] = 1

if fragment.endswith(".gz"):
    with gzip.open(fragment, "rt") as fin:
        for line in fin:
            tem = line.split()[0]
            if tem in dict_valid_barcode:
                print(tem + "\n")
else:
    with open(fragment) as fin:
        for line in fin:
            tem = line.split()[0]
            if tem in dict_valid_barcode:
                print(tem + "\n")
