#!/usr/bin/env python

# This script is based on the example at: https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv

import os
import sys
import errno
import argparse


def parse_args(args=None):
    Description = "Reformat nf-core/smrnaseq samplesheet file and check its contents."
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


def check_samplesheet(file_in, file_out):
    """
    This function checks that the samplesheet follows the following structure:
    sample,fastq_1
    SAMPLE_PE,SAMPLE_PE_RUN1_1.fastq.gz
    SAMPLE_PE,SAMPLE_PE_RUN2_1.fastq.gz
    SAMPLE_SE,SAMPLE_SE_RUN1_1.fastq.gz
    For an example see:
    https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv
    """

    sample_mapping_dict = {}
    with open(file_in, "r") as fin:

        ## Check header
        MIN_COLS = 2
        HEADER = ["sample", "fastq_1"]
        header = [x.strip('"') for x in fin.readline().strip().split(",")]
        if any([item not in header for item in HEADER]):
            missing = [item for item in HEADER if item not in header]
            eprint("ERROR: Please check samplesheet header. Missing columns: '{}'".format(",".join(missing)))
            sys.exit(1)

        ## Check sample entries
        for line_number, line in enumerate(fin):
            lspl = [x.strip().strip('"') for x in line.strip().split(",")]
            row = {k: v for k, v in zip(header, lspl)}

            # Check valid number of columns per row
            if len(lspl) != len(header):
                print_error(
                    "Invalid number of columns: found {} columns (header has {})".format(
                        line_number, len(lspl), len(header)
                    ),
                    f"Line #{line_number+2}",
                    line,
                )

            ## Check sample name entries
            sample = row.get("sample", "").replace(" ", "_")
            if not sample:
                print_error("Sample entry has not been specified!", f"Line #{line_number+2}", line)

            ## Check FastQ file extension
            fastq = row.get("fastq_1", None)
            if fastq:
                if fastq.find(" ") != -1:
                    print_error("FastQ file contains spaces!", f"Line #{line_number+2}", line)
                if not fastq.endswith(".fastq.gz") and not fastq.endswith(".fq.gz"):
                    print_error(
                        "FastQ file does not have extension '.fastq.gz' or '.fq.gz'!",
                        "Line",
                        line,
                    )

            ## Create sample mapping dictionary
            sample_info = {"single_end": "1", "fastq_1": fastq}
            if sample not in sample_mapping_dict:
                sample_mapping_dict[sample] = [sample_info]
            else:
                if sample_info in sample_mapping_dict[sample]:
                    print_error("Samplesheet contains duplicate rows!", "Line", line)
                else:
                    sample_mapping_dict[sample].append(sample_info)

    ## Write validated samplesheet with appropriate columns
    output_cols = ["id", "intrasample_id", "single_end", "fastq_1"]
    if len(sample_mapping_dict) > 0:
        out_dir = os.path.dirname(file_out)
        make_dir(out_dir)
        with open(file_out, "w") as fout:
            fout.write(",".join(output_cols) + "\n")
            for sample in sorted(sample_mapping_dict.keys()):
                for intrasample_id, val in enumerate(sample_mapping_dict[sample]):
                    sample_info = {**{"id": sample, "intrasample_id": str(intrasample_id)}, **val}
                    outrow = [sample_info.get(colname, None) for colname in output_cols]
                    fout.write(",".join(outrow) + "\n")
    else:
        print_error("No entries to process!", "Samplesheet: {}".format(file_in))


def main(args=None):
    args = parse_args(args)
    check_samplesheet(args.FILE_IN, args.FILE_OUT)


if __name__ == "__main__":
    sys.exit(main())
