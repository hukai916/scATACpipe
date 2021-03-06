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
from multiprocessing import Pool, Manager
import functools
from sinto import utils
from os.path import exists

bam             = sys.argv[1]
tagfile         = sys.argv[2]
outname         = sys.argv[3]
nproc           = int(sys.argv[4])

dict_tag = Manager().dict() # to share the dictionary across processes to save memory usage
with open(tagfile, "r") as f:
    for line in f:
        query_name, raw_barcode, corrected_barcode = line.strip().split("\t")
        if not corrected_barcode == "undetermined":
            dict_tag[query_name] = [raw_barcode, corrected_barcode]

def split_bam(infile, tag, prefix, nproc):
    inbam = pysam.AlignmentFile(infile, "rb")

    # initialise all output bam files
    outbam_files = ["".join([prefix + "_", str(i), ".bam"]) for i in range(int(nproc))]
    outbams      = [pysam.AlignmentFile(outbam_files[i], "wb", template = inbam) for i in range(int(nproc))]
    tag_dicts    = [{} for i in range(int(nproc))]

    # iterate through inbam and output to outbam according to line number
    line_number = 0
    for read in inbam.fetch():
        _pos   = line_number % int(nproc)
        outbam = outbams[_pos]
        tag_dict = tag_dicts[_pos]
        outbam.write(read)

        line_number += 1
    print("Splitting step done ...")
    # close up
    inbam.close()
    print("inbam closed")
    [outbam.close() for outbam in outbams]
    print("outbam closed")
    return(outbam_files)

def set_tag_chunk(chunk, dict_tag, tag):
    chunk_name = os.path.basename(chunk)

    print("Processing ", "chunk: " + chunk_name, " ...")
    outname = os.path.join(chunk_name + ".bam")
    outname_sorted = os.path.join("srt_" + outname)
    outname_final  = os.path.join("tagged_" + outname)

    # index to use fetch()
    pysam.sort("-o", outname_sorted, "-m", "4G", "-@", "1", chunk)
    pysam.index(outname_sorted)

    chunk   = pysam.AlignmentFile(outname_sorted, "rb")
    outbam  = pysam.AlignmentFile(outname_final, "wb", template = chunk)
    for read in chunk.fetch():
        if read.query_name in dict_tag:
            read.set_tag(tag, dict_tag[read.query_name][1])
            outbam.write(read)
    chunk.close()
    outbam.close()

    # index for later merging:
    pysam.index(outname_final)

    try:
        os.remove(outname_sorted)
        os.remove(outname_sorted + ".bai")
        os.remove(chunk_name)
    except:
        pass

    return(outname_final)

print("Splitting BAM into " + str(nproc) + " chunks ...")
chunks  = split_bam(bam, "tmp_chunk", nproc)
print("Splitted!")

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
