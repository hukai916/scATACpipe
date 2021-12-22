#!/usr/bin/env python

"""
Extract paried reads from BAM.

Usage:
python extract_pair.py inbam outbam_name.bam

"""

import sys
import pysam
import os

inname    = sys.argv[1]
outname   = sys.argv[2]

if not os.path.exists(inname + ".bai"):
    pysam.index(inname)
inbam   = pysam.AlignmentFile(inname, "rb", check_sq = False)
outbam  = pysam.AlignmentFile(outname, "wb", template = inbam)

read_dict = {}
for read in inbam.fetch():
    if not read.query_name in read_dict:
        read_dict[read.query_name] = 1
    else:
        read_dict[read.query_name] += 1
read_dict = { key:read_dict[key] for key in read_dict if read_dict[key] == 2 }

for read in inbam.fetch():
    if read.query_name in read_dict:
        outbam.write(read)

inbam.close()
outbam.close()
