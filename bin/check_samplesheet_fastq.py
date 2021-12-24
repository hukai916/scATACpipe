#!/usr/bin/env python

# TODO nf-core: Update the script to check the samplesheet
# This script is based on the example at: https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv

import os
import sys
import errno
import argparse
from itertools import combinations


def parse_args(args=None):
    Description = "Reformat nf-core/scatacseqflow samplesheet file and check its contents."
    Epilog = "Example usage: python check_samplesheet.py <FILE_IN> <FILE_OUT>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("FILE_IN", help="Input samplesheet file.")
    parser.add_argument("FILE_OUT", help="Output file.")
    return parser.parse_args(args)


def make_dir(path):
    if len(path) > 0:
        try:
            os.makedirs(path)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise exception


def print_error(error, context="Line", context_str=""):
    error_str = "ERROR: Please check samplesheet -> {}".format(error)
    if context != "" and context_str != "":
        error_str = "ERROR: Please check samplesheet -> {}\n{}: '{}'".format(
            error, context.strip(), context_str.strip()
        )
    print(error_str)
    sys.exit(1)


# TODO nf-core: Update the check_samplesheet function
def check_samplesheet(file_in, file_out):
    """
    This function checks that the samplesheet follows the following structure:

    sample_name,path_fastq_1,path_fastq_2,path_barcode
    SAMPLE_1,/Full_path/xxx.R1.fq.gz,/Full_path/xxx.R2.fq.gz,/Full_path/xxx.barcode.fq.gz
    SAMPLE_2,/Full_path/xxx.R1.fq.gz,/Full_path/xxx.R2.fq.gz,/Full_path/xxx.barcode.fq.gz

    Note that must use 'full path'. If multiple sequence files exist for a single column, seperate them by ';'.
    """

    sample_mapping_dict = {}
    with open(file_in, "r") as fin:

        ## Check header
        MIN_COLS = 3
        # TODO nf-core: Update the column names for the input samplesheet
        HEADER = ["sample_name", "path_fastq_1", "path_fastq_2", "path_barcode"]
        header = [x.strip('"') for x in fin.readline().strip().split(",")]
        if header[: len(HEADER)] != HEADER:
            print("ERROR: Please check samplesheet header -> {} != {}".format(",".join(header), ",".join(HEADER)))
            sys.exit(1)

        ## Check sample entries
        for line in fin:
            lspl = [x.strip().strip('"') for x in line.strip().split(",")]

            # Check valid number of columns per row
            if len(lspl) < len(HEADER):
                print_error(
                    "Invalid number of columns (minimum = {})!".format(len(HEADER)),
                    "Line",
                    line,
                )
            num_cols = len([x for x in lspl if x])
            if num_cols < MIN_COLS:
                print_error(
                    "Invalid number of populated columns (minimum = {})!".format(MIN_COLS),
                    "Line",
                    line,
                )

            ## Check sample name entries
            sample, fastq_1, fastq_2, barcode = lspl[: len(HEADER)]
            if sample:
                if sample.find(" ") != -1:
                    print_error("Sample entry contains spaces!", "Line", line)
            else:
                print_error("Sample entry has not been specified!", "Line", line)

            ## Check FastQ file extension
            for fastq in [fastq_1, fastq_2, barcode]:
                if fastq:
                    if fastq.find(" ") != -1:
                        print_error("FastQ file contains spaces!", "Line", line)
                    single_files = [x.strip() for x in fastq.split(";")]
                    for single_file in single_files:
                        if not single_file.endswith(".fastq.gz") and not single_file.endswith(".fq.gz"):
                            print_error(
                                "FastQ file does not have extension '.fastq.gz' or '.fq.gz'!",
                                "Line",
                                line,
                            )

            ## Auto-detect paired-end/single-end
            sample_info = []  ## [single_end, fastq_1, fastq_2]
            sample_info = [sample, fastq_1, fastq_2, barcode]

            ## Create sample mapping dictionary = { sample: [ single_end, fastq_1, fastq_2, barcode ] }
            if sample not in sample_mapping_dict:
                sample_mapping_dict[sample] = [sample_info]
            else:
                if sample_info in sample_mapping_dict[sample]:
                    print_error("Samplesheet contains duplicate rows!", "Line", line)
                else:
                    sample_mapping_dict[sample].append(sample_info)

    ## Write validated samplesheet with appropriate columns
    if len(sample_mapping_dict) > 0:
        ## Check if sample_name "contained" in another sample_name:
        sample_names = set(sample_mapping_dict.keys())
        for pair in combinations(sample_names, 2):
            if pair[0] in pair[1] or pair[1] in pair[0]:
                print_error("certain sample name 'contained' in another sample name!\nSupplied sample names: " + ", ".join(sample_names))
                # sample_1, sample_11 is not okay
                # sample_01, sample_11 is fine

        out_dir = os.path.dirname(file_out)
        make_dir(out_dir)
        with open(file_out, "w") as fout:
            fout.write(",".join(["sample_name", "path_fastq_1", "path_fastq_2", "path_barcode"]) + "\n")
            for library in sorted(sample_mapping_dict.keys()):
                for sample in sample_mapping_dict[library]:
                    fout.write(",".join(sample) + "\n")
    else:
        print_error("No entries to process!", "Samplesheet: {}".format(file_in))


def main(args=None):
    args = parse_args(args)
    check_samplesheet(args.FILE_IN, args.FILE_OUT)


if __name__ == "__main__":
    sys.exit(main())
