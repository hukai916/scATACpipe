#!/usr/bin/env python

"""
Given barcode whitelist, interleaved sample name, number of segments, and template token, and barcode token, output json filename, generate json config.

Usage:
python make_json.py path_to_whitelist sample_A.cram 3 0::,2:: 1::16 sample_A.json

For a detailed explanation of the config, see: https://learn.gencore.bio.nyu.edu/pheniqs/
"""

import json
import sys
import gzip

whitelist      = sys.argv[1]
sample_name    = sys.argv[2]
segment_num    = sys.argv[3]
template_token = sys.argv[4]
barcode_token  = sys.argv[5]
outfile        = sys.argv[6]

config = {}
config['input'] = [val for val in [sample_name] for _ in range(int(segment_num))]
config['template'] = {}
config['template']['transform'] = {}
config['template']['transform']['token'] = template_token.split(",")
config['sample'] = {}
config['sample']['algorithm'] = 'pamld'
config['sample']['confidence threshold'] = 0.99
config['sample']['noise'] = 0.01
config['sample']['transform'] = {}
config['sample']['transform']['token'] = [barcode_token]
config['sample']['codec'] = {}
if whitelist.endswith(".gz"):
    f = gzip.open(whitelist, mode='rt')
else:
    f = open(whitelist)
for line in f:
    barcode = line.split()[0]
    concentration = line.split()[1]
    config['sample']['codec']['@' + barcode] = {}
    config['sample']['codec']['@' + barcode]['barcode'] = [barcode]
    config['sample']['codec']['@' + barcode]['concentration'] = float(concentration)
f.close()

with open(outfile, 'w') as outfile:
    json.dump(config, outfile, indent = 4)
