#!/usr/bin/env python

import sys
from collections import OrderedDict

filename = sys.argv[1]

dict = OrderedDict()

for line in open(filename):
    if line.startswith("#"):
        print(line, end = "")
    else:
        cols = line.split("\t")
        for x in cols[8].split(";"):
            if "gene_id" in x.lower():
                key = cols[0] + "\t" + x.strip()[9:-1] # better handle the psudo-autosomal genes
                if not key in dict:
                    dict[key] = [line]
                else:
                    dict[key].append(line)
                break

if len(dict) == 0:
    for line in open(filename):
        print(line, end = "")
else:
    fout = open("gene_ranges.tsv", "w")
    for key in dict:
        for x in dict[key]:
            print(x, end = "")
        starts = [int(x.split("\t")[3]) for x in dict[key]]
        ends   = [int(x.split("\t")[4]) for x in dict[key]]
        fout.write(key + "\t" + str(min(starts)) + "\t" + str((max(ends))) + "\n")
        # also output a gene range file:

    fout.close()
