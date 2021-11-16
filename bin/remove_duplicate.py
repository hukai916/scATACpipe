#!/usr/bin/env python

"""
Created on Fri Sep 28 12:08:23 2021

@authors: liuh, huk

Usage:
python remove_duplicate.py -h

"""

import pysam
import sys, argparse, re, os, copy
from statistics import mean

parser = argparse.ArgumentParser(description='Remove PCR duplicates with mapping coordinates and cell barcode.')
parser.add_argument('--inbam', type = argparse.FileType('r'), default = sys.stdin, help = 'Path to paired bam file. (default: stdin)')
parser.add_argument('--outbam', type = argparse.FileType('w'), default = sys.stdout, help = 'Filename for output bam file. (default: stdout)')
parser.add_argument('--shift_forward', type = int, default = 0, help = "Number of bases to shift for the reads mapped to forward strand (+4 for Tn5). (default: %(default)i)")
parser.add_argument('--shift_reverse', type = int, default = 0, help = "Number of bases to shift for the reads mapped to reverse strand (-5 for Tn5). (default: %(default)i)")
parser.add_argument('--barcode_regex', type = str, default = '[^:]*', help = 'Regular expression (must be double quoted) that matches the barcode in name. (default: "%(default)s")')
parser.add_argument('--barcode_tag', type = str, default = 'N/A', help = "Bam tag that stores the cell barcode if available (default: %(default)s)")

args = parser.parse_args()

# read in bam file, close them when done:
pysam.index(args.inbam.name) # for inbam.fetch() to work, must be put before reading in bam.
inbam   = pysam.AlignmentFile(args.inbam, "rb")
outbam  = pysam.AlignmentFile(args.outbam, "wb", template = inbam)

# process on a per chromosome basis:
header_chr = list(inbam.references)
soft_clip_5 = re.compile(r'^(\d+)S')
soft_clip_3 = re.compile(r'(\d+)S$')
soft_clip_num = 0
no_corrected_barcode_num     = 0
not_properly_mapped_num      = 0
total_unique_fragment_num    = 0
total_duplicate_fragment_num = 0

# process in a per chromosome basis:
print("Total contig(s): ", len(header_chr))
for contig in header_chr:
    print("Processing ", contig, " ...")
    read_dict = {} # dict to store reads by name
    unique_fragments_dict = {} # dict to store unique fragments

    # fill in the read_dict dictionary:
    for read in inbam.fetch(contig = contig):
        if not read.query_name in read_dict:
            read_dict[read.query_name] = [read]
        else:
            read_dict[read.query_name].append(read)

    # accounting for soft-clipping and Tn5 shifts:
    for query_name in read_dict.keys():
        if len(read_dict[query_name]) == 2:
            # frag_end_pos = 0
            right_read_check, left_read_check = 0, 0
            not_properly_mapped_check = 0
            for read in read_dict[query_name]:
                cigar = read.cigarstring
                if read.flag in [99, 163]: # reads mapped to forward strand
                    left_read = copy.deepcopy(read) # copy the AlignmentSegment object
                    left_read_reference_end  = left_read.reference_end # reference_end is not writable, changes automatically when reference_start changes
                    left_read_check = 1
                    match_5   = re.search(soft_clip_5, cigar)
                    match_3   = re.search(soft_clip_3, cigar)
                    if match_5:
                        left_read.reference_start = left_read.reference_start - int(match_5.group(1)) + args.shift_forward
                    else:
                        left_read.reference_start = left_read.reference_start + args.shift_forward
                    if match_3:
                        left_read_reference_end = left_read.reference_end + int(match_3.group(1))
                    if match_3 or match_5:
                        soft_clip_num += 1

                    left_read.query_sequence = left_read.query_sequence[args.shift_forward:]
                    left_read.query_qualities = read.query_qualities[args.shift_forward:] # use read since left_read.query_qualities will be None as long as query_sequence is modified.
                    left_read.cigar =((0, left_read.infer_query_length() - args.shift_forward),)  # return cigar contains only 'M'.
                elif read.flag in [147, 83]: # reads mapped to reverse strand
                    right_read = copy.deepcopy(read)
                    right_read_reference_end = right_read.reference_end
                    right_read_check = 1
                    match_3    = re.search(soft_clip_3, cigar)
                    match_5    = re.search(soft_clip_5, cigar)
                    if match_3:
                        # frag_end_pos = right_read.reference_end + int(match_3.group(1)) + args.shift_reverse
                        right_read_reference_end = right_read.reference_end + int(match_3.group(1)) + args.shift_reverse
                    else:
                        # frag_end_pos = right_read.reference_end + args.shift_reverse
                        right_read_reference_end = right_read.reference_end + args.shift_reverse
                    if match_5:
                        right_read.reference_start = right_read.reference_start - int(match_5.group(1))
                    if match_3 or match_5:
                        soft_clip_num += 1

                    if not args.shift_reverse == 0:
                        right_read.query_sequence = right_read.query_sequence[:args.shift_reverse]
                        right_read.query_qualities = read.query_qualities[:args.shift_reverse]
                    # right_read.cigar = ((0, frag_end_pos - right_read.reference_start),)
                    right_read.cigar = ((0, right_read.infer_query_length() + args.shift_reverse),)
                    # infer_query_length should always equal to len(query_sequence)
                else:
                    not_properly_mapped_check = 1
                    if not not_properly_mapped_check:
                        not_properly_mapped_num += 1 # count the number of fragments that are not properly mapped
                    continue
                if left_read_check and right_read_check:
                    left_read.next_reference_start  = right_read.reference_start
                    left_read.template_length       = max(left_read_reference_end, right_read_reference_end) - min(left_read.reference_start, right_read.reference_start)
                    right_read.next_reference_start = left_read.reference_start
                    right_read.template_length      = - left_read.template_length
                    # move barcode to tag if in name:
                    cell_barcode = ''
                    try:
                        if args.barcode_tag == "N/A":
                            cell_barcode = re.search(args.barcode_regex, left_read.query_name).group()
                            left_read.set_tag(tag = "CB", value = cell_barcode, value_type = "Z")
                            right_read.set_tag(tag = "CB", value = cell_barcode, value_type = "Z")
                        else:
                            cell_barcode = left_read.get_tag(args.barcode_tag)
                    except:
                        continue # barcode can't be determined (some reads with poor quality may not have CB tag.)
                    if not cell_barcode == "": # valid fragments
                        fragment_id = ":".join([
                        cell_barcode,
                        left_read.reference_name,
                        str(min(left_read.reference_start, right_read.reference_start)), str(abs(left_read.template_length))
                        ])

                        query_quality = int(mean([x for x in left_read.query_qualities + right_read.query_qualities]))

                        if not fragment_id in unique_fragments_dict:
                            unique_fragments_dict[fragment_id] = [left_read, right_read, query_quality, 1] # last element is the count
                        else:
                            unique_fragments_dict[fragment_id][3] += 1
                            if query_quality > unique_fragments_dict[fragment_id][2]:
                                unique_fragments_dict[fragment_id][:3] = [left_read, right_read, query_quality]
                    else:
                        no_corrected_barcode_num += 1
        else:
            not_properly_mapped_num += 1 # unpaired fragments
            continue
    total_unique_fragment_num += len(unique_fragments_dict)
    total_duplicate_fragment_num += sum([unique_fragments_dict[x][-1] for x in unique_fragments_dict]) - len(unique_fragments_dict)

    for fragment in unique_fragments_dict:
        outbam.write(unique_fragments_dict[fragment][0])
        outbam.write(unique_fragments_dict[fragment][1])

    # for key in unique_fragments_dict:
    #     print(key, unique_fragments_dict[key][0].reference_start, unique_fragments_dict[key][1].reference_start)

inbam.close()
outbam.close()
print("All done!")

# print out summary statistics:
print("Soft_clipped (reads): ", soft_clip_num)
print("No corrected barcode (of all paired fragments): ", no_corrected_barcode_num)
print("Not properly mapped (fragments): ", not_properly_mapped_num)
print("Total unique fragments: ", total_unique_fragment_num)
print("Total duplicate fragments: ", total_duplicate_fragment_num, "\n")
summary = "Summary (remove_duplicate.py): total unique fragments: " + str(total_unique_fragment_num) + ", total duplicate fragments: " + str(total_duplicate_fragment_num) + ", other fragments(unproper mapped or with no corrected barcodes): " + str(no_corrected_barcode_num + not_properly_mapped_num) + "."
print(summary)
with open("summary_rm_dup_" + os.path.basename(args.inbam.name) + ".txt", "w") as f:
    f.write(summary)
