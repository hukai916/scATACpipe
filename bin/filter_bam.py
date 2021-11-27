#!/usr/bin/env python

"""
Filter bam file based on valid barcode list.
"""

import sys
import gzip
import os
import pysam

bam             = sys.argv[1]
valid_barcode   = sys.argv[2]
outbam          = sys.argv[3]

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

samfile = pysam.AlignmentFile(bam, "rb")
outbam  = pysam.AlignmentFile(outbam, "wb", template=samfile)

for read in samfile.fetch():
    cr = read.get_tag("CR")
    if cr in dict_valid_barcode:
        outbam.write(read)
outbam.close()
samfile.close()
