#!/usr/bin/env python

"""
Match genome.fa config to GTF annotation:
GTF may contain annotations that dont have matching config in the genome.fa, which is problematic to cellranger index command.

Usage:
python extract_gtf.py genome.fa xxx.gtf

"""

import sys
import re
from Bio import SeqIO
import gzip

genome = sys.argv[1]
gtf    = sys.argv[2]

if genome.endswith(".gz"):
    fasta = SeqIO.parse(gzip.open(genome, mode = 'rt'),'fasta')
else:
    fasta = SeqIO.parse(open(genome),'fasta')
contig  = [record.id for record in fasta]

# check if chrM or chrMT:
if "chrM" in contig:
    chr_mito = "chrM"
elif "chrMT" in contig:
    chr_mito = "chrMT"
else:
    chr_mito = ""

for line in open(gtf):
    chr = line.split()[0]
    if chr in contig:
        print(line, end = "")
    else:
        if chr == "chrM":
            print(re.sub("^chrM", chr_mito, line), end = "")
        elif chr == "chrMT":
            print(re.sub("^chrMT", chr_mito, line), end = "")
