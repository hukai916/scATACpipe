#!/usr/bin/env python

"""
Correct the "gene ranges" for psudo-autosomal genes.
See: https://github.com/gpertea/gffread/issues/86

Usage:
python correct_gtf_gene_range.py XXX.gff3 gene_ranges.tsv

Dev note:
## This is assuming the gff3 is sorted: features belonging to same gene are grouped together.
## Note the GFF3 output by gffread is sorted by position, so the above assumption doesn't hold. Use a gene range file inferred from gtf instead.

"""

import sys

gff3    = sys.argv[1]
gr      = sys.argv[2]


# first iterate through gr to obtain gene ranges
gene_ranges = {}
for line in open(gr):
    cols = line.split("\t")
    unique_id = cols[0] + "\t" + cols[1]
    assert unique_id not in gene_ranges, "Duplicated gene id detected!"
    gene_ranges[unique_id] = [int(cols[2]), int(cols[3])]
# second iterate GFF3 to fix the wrong ranges
for line in open(gff3):
    if line.startswith("#"):
        print(line, end = "")
    else:
        cols = line.split("\t")
        if cols[2] == "gene":
            unique_id = cols[0] + "\t" + cols[8].split("ID=")[1].split(";")[0]
            unique_id = unique_id.strip()
            # print(unique_id)
            if not unique_id in gene_ranges:
                print(unique_id)
                for key in gene_ranges:
                    print(key, key == unique_id, "\"" + key +  "\", \"" + unique_id + "\"" )
                # print(gene_ranges)
                exit()
            assert unique_id in gene_ranges, "Undetected gene_id from GTF!"
            if [int(cols[3]), int(cols[4])] != gene_ranges[unique_id]:
                cols[3], cols[4] = gene_ranges[unique_id]
                print("\t".join([str(x) for x in cols]), end = "")
            else:
                print(line, end = "")
        else:
            print(line, end = "")

### Below are legacy codes.
# first iteration to retrieve real gene ranges for each gene
# gene_ranges = {}
# unique_id = ""
# for line in open(gff3):
#     if not line.startswith("#"):
#         cols = line.split("\t")
#         if cols[2] == "gene":
#             unique_id = cols[0] + "_" + cols[8]
#             assert unique_id not in gene_ranges, "Duplicate gene_ids detected in GFF3 file!"
#             gene_ranges[unique_id] = []
#         else:
#             start = int(cols[3])
#             end   = int(cols[4])
#             assert unique_id != "", "Missing genes in GFF3 file!"
#             if gene_ranges[unique_id] == []:
#                 gene_ranges[unique_id] = [start, end]
#             else:
#                 gene_ranges[unique_id][0] = min(start, gene_ranges[unique_id][0])
#                 gene_ranges[unique_id][1] = max(end, gene_ranges[unique_id][1])

# second iteration to fix the wrong ranges
# for line in open(gff3):
#     if line.startswith("#"):
#         print(line, end = "")
#     else:
#         cols = line.split("\t")
#         if cols[2] == "gene":
#             unique_id = cols[0] + "_" + cols[8]
#             if [int(cols[3]), int(cols[4])] != gene_ranges[unique_id]:
#                 print(line, end = "")
#                 print([int(cols[3]), int(cols[4])])
#                 print(gene_ranges[unique_id], "\n")
