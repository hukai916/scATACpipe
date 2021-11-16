#!/usr/bin/env python

"""
Given dict_json, genome_name, and what to expect, return the download url.

Usage:
python get_download_url.py dict_json hg19 genome/genome_md5sum/gtf/gtf_md5sum/md5sum
"""

import json
import sys

dict_json = sys.argv[1]
genome_name = sys.argv[2]
what      = sys.argv[3]

f = open(dict_json)
dict_genome = json.load(f)

if genome_name in dict_genome:
    print(dict_genome[genome_name][what])
else:
    print(0)
