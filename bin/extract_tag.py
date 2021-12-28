#!/usr/bin/env python

"""
Extract specified tag fields from BAM.

Usage:
python extract_tag.py inbam tag1,tag2 outfile_name

"""

import sys
import pysam

bam     = sys.argv[1]
tags    = sys.argv[2]
outfile = sys.argv[3]

tag_list = tags.split(",")

inbam = pysam.AlignmentFile(bam, "rb", check_sq = False)
with open(outfile, "w") as f:
    for read in inbam:
        line = "\t".join([read.get_tag(tag) for tag in tag_list])
        line = read.get_tag(tag_list[0]) + ":" + read.query_name + "\t" + line # add raw_barcode:query_name as first column
        f.write(line + "\n")
inbam.close()
