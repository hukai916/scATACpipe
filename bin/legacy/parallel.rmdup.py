#!/usr/bin/env python

"""
Created on Fri Sep 28 12:08:23 2021
@authors: liuh, huk
Usage:
python remove_duplicate.py -h

Known issues/bugs: remove_duplicate.py has fixed them.
1. Out-of-boundary issue still exist: start can't be 0;
2. outbam better has a default;
3. If contig contains no reads, the pysam.merge() can be issue: remove_duplicate.py doesn't have such problem;
4. @ nprocs are hard-coded.
5. Need to get rid of intermediate files.

"""

import pysam
import sys, argparse, re, os, copy
from multiprocessing import Pool
from statistics import mean
import functools
import os.path
import tempfile

def rm_dup(contig, inbam, header_len_dict,
           soft_clip_5, soft_clip_3,
           shift_forward = 4, shift_reverse = -5,
           barcode_regex = "[^:]+",
           outdir = ".",
           barcode_tag = "N/A"):
    print("Processing ", contig, " ...")
    prefix = re.sub(".bam", "", os.path.basename(inbam))
    outname = os.path.join(outdir, prefix + "_" + contig + ".bam")

    read_dict = {} # dict to store reads by name
    unique_fragments_dict = {} # dict to store unique fragments

    ## summary statistics
    soft_clip_num = 0
    no_corrected_barcode_num     = 0
    not_properly_mapped_num      = 0
    total_unique_fragment_num    = 0
    total_duplicate_fragment_num = 0

    inbam   = pysam.AlignmentFile(inbam, "rb")
    outbam  = pysam.AlignmentFile(outname, "wb", template = inbam)

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
                        left_read.reference_start = left_read.reference_start - int(match_5.group(1)) + shift_forward
                    else:
                        left_read.reference_start = left_read.reference_start + shift_forward
                    if match_3:
                        left_read_reference_end = left_read.reference_end + int(match_3.group(1))
                    if match_3 or match_5:
                        soft_clip_num += 1
                    left_read.query_sequence = left_read.query_sequence[shift_forward:]
                    left_read.query_qualities = read.query_qualities[shift_forward:] # use read since left_read.query_qualities will be None as long as query_sequence is modified.
                    left_read.cigar =((0, left_read.infer_query_length() - shift_forward),)  # return cigar contains only 'M'.
                elif read.flag in [147, 83]: # reads mapped to reverse strand
                    right_read = copy.deepcopy(read)
                    right_read_reference_end = right_read.reference_end
                    right_read_check = 1
                    match_3    = re.search(soft_clip_3, cigar)
                    match_5    = re.search(soft_clip_5, cigar)
                    if match_3:
                        # frag_end_pos = right_read.reference_end + int(match_3.group(1)) + shift_reverse
                        right_read_reference_end = right_read.reference_end + int(match_3.group(1)) + shift_reverse
                    else:
                        # frag_end_pos = right_read.reference_end + shift_reverse
                        right_read_reference_end = right_read.reference_end + shift_reverse
                    if match_5:
                        right_read.reference_start = right_read.reference_start - int(match_5.group(1))
                    if match_3 or match_5:
                        soft_clip_num += 1
                    if not shift_reverse == 0:
                        right_read.query_sequence = right_read.query_sequence[:shift_reverse]
                        right_read.query_qualities = read.query_qualities[:shift_reverse]
                    # right_read.cigar = ((0, frag_end_pos - right_read.reference_start),)
                    right_read.cigar = ((0, right_read.infer_query_length() + shift_reverse),)
                    # infer_query_length should always equal to len(query_sequence)
                else:
                    not_properly_mapped_check = 1
                    if not not_properly_mapped_check:
                        not_properly_mapped_num += 1 # count the number of fragments that are not properly mapped
                    continue
                if left_read_check and right_read_check:
                    left_read.next_reference_start  = right_read.reference_start
                    left_read.template_length       = max(left_read_reference_end,
                                                          right_read_reference_end) - min(left_read.reference_start,
                                                                                          right_read.reference_start)
                    right_read.next_reference_start = left_read.reference_start
                    right_read.template_length      = -left_read.template_length
                    # move barcode to tag if in name:
                    cell_barcode = ''
                    try:
                        if barcode_tag == "N/A":
                            cell_barcode = re.search(barcode_regex, left_read.query_name).group() # this always matchs something
                            left_read.set_tag(tag = "CB", value = cell_barcode, value_type = "Z")
                            right_read.set_tag(tag = "CB", value = cell_barcode, value_type = "Z")
                        else:
                            cell_barcode = left_read.get_tag(barcode_tag)
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
                            unique_fragments_dict[fragment_id] = [left_read, right_read] # last element is the count
                            total_unique_fragment_num += 1
                        else:
                            total_duplicate_fragment_num += 1
                    else:
                        no_corrected_barcode_num += 1
        else:
            not_properly_mapped_num += 1 # unpaired fragments
            continue

    for fragment in unique_fragments_dict:
        #remove out-of-bound fragments
        left_read = unique_fragments_dict[fragment][0]
        right_read = unique_fragments_dict[fragment][1]
        if left_read.reference_start < 0 or right_read.reference_end > header_len_dict[right_read.reference_name]:
            continue
        outbam.write(left_read)
        outbam.write(right_read)

    inbam.close()
    outbam.close()
    return  [outname, dict(soft_clip_num = soft_clip_num,
                          no_corrected_barcode_num = no_corrected_barcode_num,
                          not_properly_mapped_num = not_properly_mapped_num,
                          total_unique_fragment_num = total_unique_fragment_num,
                          total_duplicate_fragment_num = total_duplicate_fragment_num)
            ]

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Remove PCR duplicates with mapping coordinates and cell barcode.')
    parser.add_argument('--inbam', type = str,
                        help = 'Path to paired bam file.')
    parser.add_argument('--outbam', type = str, default = "out.bam",
                        help = 'Filename for output bam file.')
    parser.add_argument('--outdir', type = str,
                        help = 'directory for output intermediate bam file.')
    parser.add_argument('--nproc', type = int, default = 1,
                        help = "Number of cores for parallel computing. (default: %(default)i)")
    parser.add_argument('--shift_forward', type = int, default = 4,
                        help = "Number of bases to shift for the reads mapped to forward strand (+4 for Tn5). (default: %(default)i)")
    parser.add_argument('--shift_reverse', type = int, default = -5,
                        help = "Number of bases to shift for the reads mapped to reverse strand (-5 for Tn5). (default: %(default)i)")
    parser.add_argument('--barcode_regex', type = str, default = '[^:]*',
                        help = 'Regular expression (must be double quoted) that matches the barcode in name. (default: "%(default)s")')
    parser.add_argument('--barcode_tag', type = str, default = 'N/A',
                        help = "Bam tag that stores the cell barcode if available (default: %(default)s)")

    args = parser.parse_args()

    # read in bam file, close them when done:
    if (os.path.exists(args.inbam)):
        # for inbam.fetch() to work, must be put before reading in bam.
        try:
            inbam   = pysam.AlignmentFile(args.inbam, "rb")
        except:
            pysam.index(args.inbam)
            inbam   = pysam.AlignmentFile(args.inbam, "rb")
    else:
        raise ValueError(args.inbam + " doesn't exist!")

    # process on a per chromosome basis:
    header_chr = list(inbam.references)
    header_length =  list(inbam.lengths)
    header_len_dict = dict(zip(header_chr , header_length))

    soft_clip_5 = re.compile(r'^(\d+)S')
    soft_clip_3 = re.compile(r'(\d+)S$')
    inbam.close()

    if (not os.path.exists(args.outdir)):
        os.makedirs(args.outdir)

    nproc = args.nproc
    with Pool(nproc) as p:
        contig_bam_lists = p.map_async(functools.partial(rm_dup, inbam = args.inbam,
                                                    header_len_dict = header_len_dict,
                                                    soft_clip_5 = soft_clip_5,
                                                    soft_clip_3 = soft_clip_3,
                                                    shift_forward = args.shift_forward,
                                                    shift_reverse = args.shift_reverse,
                                                    outdir = args.outdir,
                                                    barcode_regex = args.barcode_regex,
                                                    barcode_tag = args.barcode_tag),
                                    header_chr,
                                    chunksize = len(header_chr) // nproc).get()

    soft_clip_num = 0
    no_corrected_barcode_num     = 0
    not_properly_mapped_num      = 0
    total_unique_fragment_num    = 0
    total_duplicate_fragment_num = 0

    contig_bam_files = [res[0] for res in contig_bam_lists] # list of BAM filenames
    contig_bam_stats = [res[1] for res in contig_bam_lists] # list of dict
    # cat files and write to output

    out_bams = []
    for bam in contig_bam_files:
        prefix = os.path.basename(bam)
        prefix = re.sub(".bam", "", prefix)
        outname =  os.path.join(args.outdir, prefix + ".srt.bam")
        pysam.sort("-o", outname,"-m", "4G", "-@", "4", bam)
        pysam.index(outname)
        out_bams.append(outname)
        #os.remove(i)

    merge_param = ["-f", "-@", "8", args.outbam] +  out_bams
    pysam.merge(*merge_param)
    pysam.index(args.outbam)

    for stat in contig_bam_stats:
        soft_clip_num += stat[soft_clip_num]
        no_corrected_barcode_num += stat["no_corrected_barcode_num"]
        not_properly_mapped_num += stat["not_properly_mapped_num"]
        total_unique_fragment_num += stat["total_unique_fragment_num"]
        total_duplicate_fragment_num += stat["total_duplicate_fragment_num"]

    # print out summary statistics:
    print("Soft_clipped (reads): ", soft_clip_num)
    print("No corrected barcode (of all paired fragments): ", no_corrected_barcode_num)
    print("Not properly mapped (fragments): ", not_properly_mapped_num)
    print("Total unique fragments: ", total_unique_fragment_num)
    print("Total duplicate fragments: ", total_duplicate_fragment_num, "\n")
    summary = ("Summary (remove_duplicate.py): total unique fragments: " +
               str(total_unique_fragment_num) +
               ", total duplicate fragments: " + str(total_duplicate_fragment_num) +
               ", other fragments(unproper mapped or with no corrected barcodes): " +
               str(no_corrected_barcode_num + not_properly_mapped_num) + ".")
    print(summary)
    with open("summary_rm_dup_" + os.path.basename(args.inbam) + ".txt", "w") as f:
        f.write(summary)
