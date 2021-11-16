#!/usr/bin/env python

"""
Match bed to GTF col1:

Usage:
python extract_bed.py fragment.bed annotation.gtf

"""

import sys
import gzip

bed = sys.argv[1]
gtf = sys.argv[2]

chrList = []
for line in open(gtf):
    tem = line.split()[0]
    if not tem in chrList:
        chrList.append(tem)

for line in open(bed):
    tem = line.split()[0]
    if tem in chrList:
        print(line, end = "")
