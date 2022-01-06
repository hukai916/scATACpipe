#!/usr/bin/env python

"""
Convert fragment file to barcode_fragment_frequency file.
Usage:

python frag_to_freq.py fragment_file out_filename
"""

import sys
import gzip

f_frag = sys.argv[1]
f_out  = sys.argv[2]

if f_frag.endswith(".gz"):
    f_in = gzip.open(f_frag, "rb")
else:
    f_in = open(f_frag)

freq_dict = {}
for line in f_in:
    barcode, count = line.strip().split()[-2:]
    if not barcode in freq_dict:
        freq_dict[barcode] = int(count)
    else:
        freq_dict[barcode] += int(count)

freq_dict_sorted = {k: v for k, v in sorted(freq_dict.items(), reverse=True, key=lambda item: item[1])}
with open(f_out, "w") as f:
  for key in freq_dict_sorted:
      print(key, freq_dict_sorted[key])
      f.write(key + "\t" + str(freq_dict_sorted[key]) + "\n")

f_in.close()
