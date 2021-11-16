#!/usr/bin/env python

"""
Convert gff3 to gtf.

Usage:
python gff3_to_gtf.py xxx.gff3 (output will be xxx.gtf)

"""

import sys
from bioinfokit.analys import gff

infile_gff3 = sys.argv[1]

gff.gff_to_gtf(file=infile_gff3)
