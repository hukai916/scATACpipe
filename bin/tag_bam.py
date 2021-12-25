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
            dict_tag[query_name] = [raw_barcode, correct_barcode]

# If using custom_chunk_bam, need to add the header too, below won't work due to lack of header: samtools view fails
# def custom_chunk_bam(bam, nproc):
#     if not exists(bam + ".bai"): # assuming position sorted, otherwise index won't work.
#         pysam.index(bam)
#     total_reads = int(pysam.idxstats(bam).split()[-1])
#     chunk_reads = math.ceil(total_reads / nproc)
#     command = "samtools view " + bam + " | split --lines=" + chunk_reads + " --filter='samtools view -o \${FILE}.bam' - tmp_"
#     subprocess.call(command, shell = True)
#     bam_chunks = [f for f in os.listdir('.') if re.match(r'tmp_*\.bam', f)]
#     return(bam_chunks)

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

# def set_tag_interval(intervals, inbam, dict_tag, tag):
#     inbam   = pysam.AlignmentFile(inbam, "rb")
#     temp_name_start = "_".join([str(intervals[0][i]) for i in [0, 1]])
#     temp_name_end   = "_".join([str(intervals[-1][i]) for i in [0, 2]])
#     temp_name       = temp_name_start + "_to_" + temp_name_end
#
#     print("Processing ", "chunk_" + temp_name, " ...")
#     prefix  = re.sub(".bam$", "", os.path.basename(bam))
#     outname = os.path.join("tagged_" + prefix + "_chunk_" + temp_name + ".bam")
#     outbam  = pysam.AlignmentFile(outname, "wb", template = inbam)
#
#     for i in intervals:
#         for read in inbam.fetch(i[0], i[1], i[2]):
#             if read.query_name in dict_tag:
#                 read.set_tag(tag, dict_tag[read.query_name][1])
#                 outbam.write(read)
#                 # raw_barcode = re.search('[^:]*', read.query_name).group()
#                 # if raw_barcode in dict_tag:
#                 #     read.set_tag(tag, dict_tag[raw_barcode])
#                 #     outbam.write(read)
#             # else:
#             #     print(raw_barcode)
#
#     inbam.close()
#     outbam.close()
#     # Index for later merging:
#     outname_sorted = os.path.join("tagged_" + prefix + "_chunk_" + temp_name + ".srt.bam")
#     pysam.sort("-o", outname_sorted, "-m", "4G", "-@", "1", outname)
#     pysam.index(outname_sorted)
#     try:
#         os.remove(outname)
#     except:
#         pass
#
#     return(outname_sorted)

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

inbam   = pysam.AlignmentFile(bam, "rb")
chunks  = split_bam(inbam, "tmp_chunk", nproc)
inbam.close()

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
