#!/usr/bin/env python

"""
To mimick sinto::addbarcodes() with support of multiprocessing. Assuming the orders of records match in fq1, fq2, and fq3.

Usage:
python addbarcodes_parallel.py fq1 f2 f3

# fq1 refers to barcode fastq in .gz format.
# fq2 refers to R1 fastq in .gz format.
# fq3 refers to R2 fastq in .gz format.

Dev notes:
1. may try convert fastq to bam, and indexing bam for fast range retrival?
"""

import sys
from multiprocessing import Pool
import functools
import gzip
import subprocess
import math

fq1 = sys.argv[1] # barcode fastq
barcode_length = int(sys.argv[2]) # barcode length
fq2 = sys.argv[3] # R1 fastq
fq3 = sys.argv[4] # R2 fastq
nproc = int(sys.argv[5]) # number of cpu

# check if fastq ends with .gz:
assert fq1.endswith(".gz"), "Currently only .gz file supported!"
assert fq2.endswith(".gz"), "Currently only .gz file supported!"
assert fq3.endswith(".gz"), "Currently only .gz file supported!"

# check if lengths match:
n1 = int(subprocess.run("zcat " + fq1 + " | wc -l", shell = True, text = True, capture_output = True).stdout)
n2 = int(subprocess.run("zcat " + fq2 + " | wc -l", shell = True, text = True, capture_output = True).stdout)
n3 = int(subprocess.run("zcat " + fq3 + " | wc -l", shell = True, text = True, capture_output = True).stdout)
assert n1 == n2 == n3, "Fastq read numbers differ!"

# check if line number is divisible by 4:
assert n1 % 4 == 0, "Not valid Fastq format!"

def closest_number(n, m):
    """
    Find the closest number to n that is no less than n and also divisible by m.
    Ref: https://www.geeksforgeeks.org/find-number-closest-n-divisible-m/
    """
    q = int(n / m)
    n1 = m * q
    if((n * m) > 0) :
        n2 = (m * (q + 1))
    else :
        n2 = (m * (q - 1))

    if (abs(n - n1) < abs(n - n2)) :
        if n1 < n:
            n1 += m
        return n1
    if n2 < n:
        n2 += m
    return n2

def get_barcodes(f, interval, bases=12, prefix="", suffix=""):
    """
    Use generator to save memory.
    """
    f_open = gzip.open(f, "rb")
    if f.endswith(".gz"):
        gz = True
    else:
        gz = False

    for index,line in enumerate(f_open):
        if index in interval:
            if (index % 4 == 1):
                if gz:
                    yield prefix + line.decode("utf-8")[:bases] + suffix
                    # cb.append(prefix + line.decode("utf-8")[:bases] + suffix)
                else:
                    yield prefix + line[:bases] + suffix
                    # cb.append(prefix + line[:bases] + suffix)
    f_open.close()
    # return(cb)

def _add_barcode(interval, cb, fq):
    o = "tmp_" + str(interval[0]) + "_" + fq.replace(".fastq.gz", "").replace(".fq.gz", "") + ".barcoded.fastq.gz"
    outfile = gzip.GzipFile(o, mode = "wb")
    with gzip.open(fq, "r") as f:
        for index, line in enumerate(f):
            if index in interval:
                if (index % 4 == 0):
                    rdname = "@" + next(cb) + ":" + line.decode("utf-8")[1:]
                    outfile.write(rdname.encode("utf-8"))
                else:
                    outfile.write(line)
    return o

def add_barcode(intervals, fq1, fq2, fq3, bases = 16):
    # return ["adsf", 'asasdfsdfadsf']
    cb = get_barcodes(f = fq1, interval = intervals, bases = bases)
    R1_barcoded = _add_barcode(intervals, cb, fq2)
    cb = get_barcodes(f = fq1, interval = intervals, bases = bases)
    R2_barcoded = _add_barcode(intervals, cb, fq3)

    return [R1_barcoded, R2_barcoded]

# get fastq intervals
chunk_size = math.ceil(n1 / nproc)
chunk_size = closest_number(chunk_size, 4)

# ensure chunk_size is a fold of 4
def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]
intervals = list(chunks(range(0, n1), chunk_size))
intervals = {i: intervals[i] for i in range(0, len(intervals))}

# add barcode for each chunk:
p = Pool(nproc)
tempfiles = p.map_async(
    functools.partial(
        add_barcode,
        fq1 = fq1,
        fq2 = fq2,
        fq3 = fq3,
        bases = barcode_length),
        intervals.values()).get()

# cat tmp fastq into final fastq.gz and cleanup tmp files:
tmp_R1, tmp_R2 = list(list(zip(*tempfiles))[0]), list(list(zip(*tempfiles))[1])
o_r1 = fq2.replace(".fastq.gz", "").replace(".fq.gz", "") + ".barcoded.fastq.gz"
o_r2 = fq3.replace(".fastq.gz", "").replace(".fq.gz", "") + ".barcoded.fastq.gz"
subprocess.run("cat " + " ".join(tmp_R1) + " > " + o_r1, shell = True, text = True, capture_output = True)
subprocess.run("rm " + " ".join(tmp_R1), shell = True, text = True, capture_output = True)
subprocess.run("cat " + " ".join(tmp_R2) + " > " + o_r2, shell = True, text = True, capture_output = True)
subprocess.run("rm " + " ".join(tmp_R2), shell = True, text = True, capture_output = True)
