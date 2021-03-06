#!/usr/bin/env python

"""
Add tag CB to bam file via a tagfile: for tag pair (raw-corrected) in the tag file, assign the corrected tag to the CB field of input bam where the raw matches the one in the name line of the input bam.

Usage:
python tag_bam.py in.bam tag.txt outname.bam

Legacy version: split bam by read chunk

"""

import sys
import pysam
import re
import subprocess
import math
import os
from multiprocessing import Pool

bam             = sys.argv[1]
tagfile         = sys.argv[2]
outname         = sys.argv[3]
nproc           = sys.argv[4]

dict_tag = {}
with open(tagfile, "r") as f:
    for line in f:
        tem = line.split()
        raw_barcode = tem[0]
        corrected_barcode = tem[1]
        if not corrected_barcode == "undetermined":
            if not raw_barcode in dict_tag:
                dict_tag[raw_barcode] = corrected_barcode

def chunk_bam(bam, nproc):
    pysam.index(bam)
    total_reads = int(pysam.idxstats(bam).split()[-1])
    chunk_reads = math.ceil(total_reads / nproc)
    command = "samtools view " + bam + " | split --lines=" + chunk_reads + " --filter='samtools view -o \${FILE}.bam' - tmp_"
    subprocess.call(command, shell = True)
    bam_chunks = [f for f in os.listdir('.') if re.match(r'tmp_*\.bam', f)]
    return(bam_chunks)

def set_tag(bam, dict_tag, tag):
    inbam   = pysam.AlignmentFile(bam, "rb")
    prefix  = re.sub(".bam$", "", os.path.basename(bam))
    outname = os.path.join("tagged_" + prefix + ".bam")
    outbam  = pysam.AlignmentFile(outname, "wb", template = inbam)

    for read in bam.fetch():
        raw_barcode = re.search('[^:]*', read.query_name).group()
        if raw_barcode in dict_tag:
            read.set_tag(tag, dict_tag[raw_barcode])
            outbam.write(read)
    inbam.close()
    outbam.close()
    # Index for later merging:
    outname_sorted = os.path.join("tagged_" + prefix + ".srt.bam")
    pysam.sort("-o", outname_sorted, "-m", "4G", "-@", "1", outname)
    pysam.index(outname_sorted)
    try:
        os.remove(outname)
    except:
        pass

    return(outname_sorted)

# for read in inbam.fetch():
#     raw_barcode = re.search('[^:]*', read.query_name).group()
#     if raw_barcode in dict_tag:
#         read.set_tag("CB", dict_tag[raw_barcode])
#         outbam.write(read)
inbam   = pysam.AlignmentFile(bam, "rb")
outbam  = pysam.AlignmentFile(outname, "wb", template = inbam)

chunk_bams = chunk_bam(inbam, nproc)
with Pool(nproc) as p:
    chunk_bam_lists = p.map_async(
        functools.partial(set_tag,
                          dict_tag = dict_tag,
                          tag      = "CB"
                          ),
                          chunk_bams
                          ).get()
inbam.close()
outbam.close()

merge_param = ["-f", "-@", str(nproc), outname] +  out_bams
pysam.merge(*merge_param)

# clean intermediate files:
for bam in chunk_bam_lists:
    try:
        os.remove(bam)
        os.remove(bam + ".bai")
    except:
        pass
