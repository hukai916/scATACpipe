#!/usr/bin/env python

"""
Add tag CB to bam file via a tagfile: for tag pair (raw-corrected) in the tag file, assign the corrected tag to the CB field of input bam where the raw matches the one in the name line of the input bam.

Usage:
python tag_bam.py in.bam tag.txt outname.bam
The tag.txt follows:
query_name raw_barcode corrected_barcode

"""

import sys
import pysam
import re
import subprocess
import math
import os
from multiprocessing import Pool
import functools
from sinto import utils
from os.path import exists

bam             = sys.argv[1]
tagfile         = sys.argv[2]
outname         = sys.argv[3]
nproc           = int(sys.argv[4])

dict_tag = {}
with open(tagfile, "r") as f:
    for line in f:
        tem = line.split()
        query_name      = tem[0]
        raw_barcode     = tem[1]
        corrected_barcode = tem[2]
        if not corrected_barcode == "undetermined":
            dict_tag[query_name] = [raw_barcode, corrected_barcode]

def split_bam(infile, prefix, nproc):
    inbam = pysam.AlignmentFile(infile, "rb")

    # initialise all output bam files
    outbams = [pysam.AlignmentFile("".join([prefix, "_", str(i), ".bam"]), "wb", template = inbam) for i in range(int(nproc))]

    # iterate through inbam and output to outbam according to line number
    line_number = 0
    for read in inbam.fetch():
        outbam = outbams[line_number % int(nproc)]
        outbam.write(read)
        line_number += 1

    # close up
    inbam.close()
    [outbam.close() for outbam in outbams]

    return(outbams)

def set_tag_chunk(chunk, dict_tag, tag):
    chunk_name = os.path.basename(chunk)

    print("Processing ", "chunk: " + chunk_name, " ...")
    outname = os.path.join("tagged_" + chunk_name + ".bam")
    outbam  = pysam.AlignmentFile(outname, "wb", template = chunk)
    for read in chunk.fetch():
        if read.query_name in dict_tag:
            read.set_tag(tag, dict_tag[read.query_name][1])
            outbam.write(read)
    outbam.close()

    # index for later merging:
    outname_sorted = os.path.join("srt_" + outname)
    pysam.sort("-o", outname_sorted, "-m", "4G", "-@", "1", outname)
    pysam.index(outname_sorted)
    try:
        os.remove(outname)
    except:
        pass

    return(outname_sorted)

chunks  = split_bam(bam, "tmp_chunk", nproc)

with Pool(nproc) as p:
    chunk_bam_lists = p.map_async(
        functools.partial(set_tag_chunk,
                          dict_tag = dict_tag,
                          tag      = "CB"
                          ),
                          chunks
                          ).get()

merge_param = ["-f", "-@", str(nproc), outname] +  chunk_bam_lists
pysam.merge(*merge_param)

# clean intermediate files:
for bam in chunk_bam_lists:
    try:
        os.remove(bam)
        os.remove(bam + ".bai")
    except:
        pass
