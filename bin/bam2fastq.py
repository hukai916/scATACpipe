#!/usr/bin/env python

"""
Convert pheniqs output bam into fastq and save a summary statistics file.
Extract fastq only from pheniqs output bam files.

Usage:
python bam2fastq.py XXX.bam outfile_name
"""

import pysam
import sys
import gzip

bamfile = sys.argv[1]
outname = sys.argv[2]
inbam = pysam.AlignmentFile(bamfile, check_sq = False)

# some statistics
valid_read_num   = 0
rescued_read_num = 0
discard_read_num = 0
count            = 0

# output R1 and R2 fastq:
r1 = gzip.open(outname + ".R1.fastq.gz", "wt")
r2 = gzip.open(outname + ".R2.fastq.gz", "wt")
r3 = gzip.open(outname + ".R3.fastq.gz", "wt") # stores corrected index

for read in inbam:
    count += 1
    if count % 100000 == 0:
        print("Processed ", int(count / 100000) * 100000, " reads ...")
    tag_rg = read.get_tag("RG") # stores the corrected barcode from whitelist
    tag_bc = read.get_tag("BC") # stores the sequenced barcode (note BC is intentionally used to store cell barcode here)
    tag_qt = read.get_tag("QT") # stores the index sequencing quality
    if tag_rg == "undetermined":
        discard_read_num += 1
    else:
        if tag_bc == tag_rg:
            valid_read_num += 1
        else:
            rescued_read_num += 1
        line1 = "".join(["@", tag_rg, ":", read.query_name])
        line2 = read.query_sequence
        line3 = "+"
        line4 = read.qual
        if (read.flag & 64): # bitwise operation to see if read is first in pair
            r1.write("\n".join([line1, line2, line3, line4]) + "\n")
        elif (read.flag & 128):
            r2.write("\n".join([line1, line2, line3, line4]) + "\n")
        r3.write("\n".join([line1, tag_rg, "+", tag_qt]) + "\n")
inbam.close()
r1.close()
r2.close()
r3.close()

# store a summary file:
print("Done!\n")
print("Total valid reads: ", valid_read_num)
print("Total rescued reads: ", rescued_read_num)
print("Total discarded reads: ", discard_read_num)

summary = "Summary (R_correct_barcode): " + "total valid: " + str(valid_read_num) + "; total corrected: " + str(rescued_read_num) + "; total discarded: " + str(discard_read_num) + "."
with open("summary_" + outname + ".txt", "w") as f:
    f.write(summary)
