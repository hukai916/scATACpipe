#!/usr/bin/env python

"""
Created on Fri Sep 28 12:08:23 2021
@authors: liuh, huk

Usage:
python remove_duplicate.py -h

Dev notes:
1. After extending for soft-clips and Tn5 shifts, need to ensure reads don't exceed chr boundaries, otherwise samtools sort/index can be problematic.
  1.1 Instead of discarding those reads, fixing them by not 'extending' too much.
  1.2 The pos field must not be 0 (1-based SAM, not everything you can see is 1-based, but in BAM, it is 0-based) in resulting BAM file, otherwise samtools sort/index can be problematic.
2. When shifting for Tn5, must shift positive value for forward reads and negative for reverse reads: this is to comply with boudary check.
  2.1 When shifting, for forward reads: start - shift_forward
  2.2 When shifting, for reverse reads: end + shift_reverse (shift_reverse is negative)
3. For paired-read check, since we split BAM into smaller chunks, "len(read_dict[query_name]) == 2" criteria is no longer valid (R1/R2 may end up to different chunks).
  3.1 To fix this, 'extend' both ends (only extend right-end should be fine) of the chunk by max_frag_len in order to include paired reads.

"""

import pysam
import sys, argparse, re, os, copy
from multiprocessing import Pool
from statistics import mean
import functools
import os.path
import tempfile
from sinto import utils
import math

def rm_dup(intervals, inbam, header_len_dict,
           soft_clip_5, soft_clip_3,
           shift_forward = 4, shift_reverse = -5,
           barcode_regex = "[^:]*",
           outdir = ".",
           barcode_tag = "N/A",
           max_frag_len = 5000):

    temp_name_start = "_".join([str(intervals[0][i]) for i in [0, 1]])
    temp_name_end   = "_".join([str(intervals[-1][i]) for i in [0, 2]])
    temp_name       = temp_name_start + "_to_" + temp_name_end

    print("Processing ", "chunk_" + temp_name, " ...")
    prefix = re.sub(".bam$", "", os.path.basename(inbam))
    outname = os.path.join(outdir, "tmp_" + prefix + "_chunk_" + temp_name + ".bam")

    ## summary statistics:
    soft_clip_num = 0
    no_corrected_barcode_num     = 0
    not_properly_mapped_num      = 0
    total_unique_fragment_num    = 0
    total_duplicate_fragment_num = 0

    inbam   = pysam.AlignmentFile(inbam, "rb")
    outbam  = pysam.AlignmentFile(outname, "wb", template = inbam)

    # fill in the read_dict dictionary:
    for i in intervals:
        read_dict = {} # dict to store reads by name
        unique_fragments_dict = {} # dict to store unique fragments

        s_extend = max(0, i[1] - max_frag_len) # pysam regin is 0-based.
        e_extend = min(header_len_dict[i[0]], i[2] + max_frag_len)

        for read in inbam.fetch(i[0], s_extend, e_extend):
            if not read.query_name in read_dict:
                read_dict[read.query_name] = [read]
            else:
                read_dict[read.query_name].append(read)

        # keep only those reads that the left_read_start is within chunk range:
        tokeep = []
        for query_name in read_dict.keys():
            for read in read_dict[query_name]:
                if read.flag in [99, 163]: # reads mapped to forward strand
                    if read.reference_start in range(math.floor(i[1]), math.ceil(i[2])):
                        tokeep.append(query_name)
                        break

        read_dict = { query_name: read_dict[query_name] for query_name in tokeep }

        # accounting for soft-clipping and Tn5 shifts:
        for query_name in read_dict.keys():
            #if len(read_dict[query_name]) == 2: since we count reads by bam chunks, this criteria is no longer valid for determining paired reads.
            if len(read_dict[query_name]) == 2:
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

                        # ensure reads are within boundaries:
                        left_read.reference_start = max(1, left_read.reference_start) # reference start must not be 0.
                        left_read_reference_end   = min(header_len_dict[left_read.reference_name], left_read_reference_end) # reference end must not greater than chr end.

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

                        # ensure reads are within boundaries:
                        right_read.reference_start = max(1, right_read.reference_start) # reference start must not be 0.
                        right_read_reference_end   = min(header_len_dict[right_read.reference_name], right_read_reference_end) # reference

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
                not_properly_mapped_num += 1

        for fragment in unique_fragments_dict:
            left_read = unique_fragments_dict[fragment][0]
            right_read = unique_fragments_dict[fragment][1]
            outbam.write(left_read)
            outbam.write(right_read)

    inbam.close()
    outbam.close()

    outname_sorted = os.path.join(outdir, "tmp_" + prefix + "_chunk_" + temp_name + ".srt.bam")
    pysam.sort("-o", outname_sorted, "-m", "4G", "-@", "1", outname)
    pysam.index(outname_sorted)
    try:
        os.remove(outname)
        os.remove(outname + ".bai")
    except:
        pass

    return  [outname_sorted, dict(soft_clip_num = soft_clip_num,
                          no_corrected_barcode_num = no_corrected_barcode_num,
                          not_properly_mapped_num = not_properly_mapped_num,
                          total_unique_fragment_num = total_unique_fragment_num,
                          total_duplicate_fragment_num = total_duplicate_fragment_num)
            ]

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Remove PCR duplicates with mapping coordinates and cell barcode.')
    parser.add_argument('--inbam', type = str,
                        help = 'Path to paired bam file.')
    parser.add_argument('--outbam', type = str, default = "dedup.bam",
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
    parser.add_argument('--max_frag_len', type = int, default = 5000,
                        help = "When chunking BAM for parallel computing, extend MAX_FRAG_LEN for each chunk at both ends to include both reads in the fragment.")


    args = parser.parse_args()

    assert args.shift_forward >= 0, "Can not shift negative value for forward reads!"
    assert args.shift_reverse <= 0, "Can not shift positive value for reverse reads!"
    # Above is to ensure the out-of-boundary issue doesn't exist, otherwise samtools sort/index might be problematic.

    # read in bam file, close them when done:
    if os.path.exists(args.inbam):
        if not os.path.exists(args.inbam + ".bai"):
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
    bam_temp = pysam.AlignmentFile(args.inbam, "rb")
    # intervals = utils.chunk_bam(bam_temp, nproc)
    intervals = utils.chunk_bam_by_chr(bam_temp, nproc)
    bam_temp.close()

    with Pool(nproc) as p:
        chunk_bam_lists = p.map_async(
            functools.partial(rm_dup, inbam = args.inbam,
                            header_len_dict = header_len_dict,
                            soft_clip_5 = soft_clip_5,
                            soft_clip_3 = soft_clip_3,
                            shift_forward = args.shift_forward,
                            shift_reverse = args.shift_reverse,
                            outdir = args.outdir,
                            barcode_regex = args.barcode_regex,
                            barcode_tag = args.barcode_tag,
                            max_frag_len = args.max_frag_len),
                            intervals.values()).get()

    soft_clip_num = 0
    no_corrected_barcode_num     = 0
    not_properly_mapped_num      = 0
    total_unique_fragment_num    = 0
    total_duplicate_fragment_num = 0

    chunk_bam_files  = [res[0] for res in chunk_bam_lists] # list of BAM filenames
    chunk_bam_stats = [res[1] for res in chunk_bam_lists] # list of dict
    # cat files and write to output

    out_bams = [bam for bam in chunk_bam_files]
    # for bam in chunk_bam_files:
    #     prefix = os.path.basename(bam)
    #     prefix = re.sub(".bam$", "", prefix)
    #     outname =  os.path.join(args.outdir, prefix + ".srt.bam")
    #     pysam.sort("-o", outname, "-m", "4G", "-@", str(nproc), bam)
    #     pysam.index(outname)
    #     out_bams.append(outname)

    merge_param = ["-f", "-@", str(nproc), args.outbam] +  out_bams
    pysam.merge(*merge_param)

    # clean intermediate files:
    for bam in chunk_bam_files:
        try:
            os.remove(bam)
            os.remove(bam + ".bai")
        except:
            pass

    for stat in chunk_bam_stats:
        soft_clip_num += stat["soft_clip_num"]
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
