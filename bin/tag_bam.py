#!/usr/bin/env python

"""
Add tag CB to bam file via a tagfile: for tag pair (raw-corrected) in the tag file, assign the corrected tag to the CB field of input bam where the raw matches the one in the name line of the input bam.

Usage:
python tag_bam.py in.bam tag.txt outname.bam

"""

import sys
import pysam
import re

bam             = sys.argv[1]
tagfile         = sys.argv[2]
outname         = sys.argv[3]

dict_tag = {}
with open(tagfile, "r") as f:
    for line in open(f):
        tem = line.split()
        raw_barcode = tem[0]
        corrected_barcode = tem[1]
        if not corrected_barcode == "undetermined":
            if not raw_barcode in dict_tag:
                dict_tag[raw_barcode] = corrected_barcode

inbam   = pysam.AlignmentFile(bam, "rb")
outbam  = pysam.AlignmentFile(outname, "wb", template = inbam)

for read in inbam.fetch():
    raw_barcode = re.search('[^:]*', read.query_name).group()
    if raw_barcode in dict_tag:
        read.set_tag("CB", dict_tag[raw_barcode])
        outbam.write(read)

inbam.close()
outbam.close()
