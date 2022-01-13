#!/usr/bin/env python

"""
Filter out features in GTF that go beyond the boundaries of its genomic contig.
Usage:
python extract_gtf.py genome.fa XXX.gtf

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
contig_length = dict([(record.id, len(record.seq)) for record in fasta])

for line in open(gtf):
    chr = line.split()[0]
    end = int(line.split()[4])
    assert chr in contig_length, "Feature contig not found in genome!"
    if end <= contig_length[chr]:
        print(line, end = "")
