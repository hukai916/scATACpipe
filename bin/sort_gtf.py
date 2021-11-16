#!/usr/bin/env python

import sys
from collections import OrderedDict

filename = sys.argv[1]

dict = OrderedDict()

for line in open(filename):
    if line.startswith("#"):
        print(line, end = "")
    else:
        for x in line.split("\t")[8].split(";"):
            if "gene_id" in x.lower():
                if not x in dict:
                    dict[x] = [line]
                else:
                    dict[x].append(line)
                break

if len(dict) == 0:
    for line in open(filename):
        print(line, end = "")
else:
    for key in dict:
        for x in dict[key]:
            print(x, end = "")
