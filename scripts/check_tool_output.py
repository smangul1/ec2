#!/usr/bin/python

import os 
import argparse
import gzip
from Bio import SeqIO
import pickle

#kmer_sizes = [27, 31]
kmer_sizes = [15]
#kmer_sizes = [15, 19, 24, 27, 31]

missing_files_log = "missing_files.log"

def count_gzipped_reads(f, f_format="fastq"):
    read_count = 0
    with gzip.open(f, 'r') as handle:
        for record in SeqIO.parse(handle, f_format):
            read_count += 1
        return read_count


def count_reads(f, f_format="fastq"):
    read_count = 0
    
    f2, ext = os.path.splitext(f)
    # if file is gzipped, need to open differently 
    if ext == ".gz":
        f_format = os.path.splitext(f2)[1].split(".")[1]
        with gzip.open(f, 'r') as handle:
            for record in SeqIO.parse(handle, f_format):
                read_count += 1
    
    # else assume it's not
    else:
        f_format = ext.split(".")[1]
        with open(f, 'r') as handle:
            for record in SeqIO.parse(handle, f_format):
                read_count += 1
    return read_count




def parse_filename(filename):
    name, ext = os.path.splitext(filename)
    split_name = name.split('_')
    dataset = split_name[0]
    read_len = split_name[3]
    coverage = split_name[5]
    read_end = split_name[6]
    return (dataset, read_len, coverage, read_end)

def main():
    parser = argparse.ArgumentParser(description="Checks output of Error Correction tools")
    parser.add_argument('-i', '--input-dir', dest='input_dir', required=True, help='Input RAW DATA directory')
    parser.add_argument('-t', '--tool', dest='tool', required=True, help='Name of tool to check')
    parser.add_argument('-o', '--output-dir', dest='output_dir', required=True, help='Output directory containing corrected reads')
    args = parser.parse_args()

    input_dir = args.input_dir
    tool_name = args.tool
    output_dir = args.output_dir

    if not os.path.exists(input_dir):
        print "Error", input_dir, "does not exist"
        exit()

    if not os.path.exists(output_dir):
        print "Error", output_dir, "does not exists"
        exit()

    # example: TCRA_sim_rl_100_cov_4_R2.fastq
    output_string_format = '{}_sim_rl_{}_cov_{}_{}'

    read_lengths = set()
    coverages = set()
    datasets = set()
    read_ends = set()
    other_files = set()
    input_files = []
    input_file_readcount = {}

    # for each file in input directory
    for root, dirs, files in os.walk(input_dir):
        for f in files: 
            print f
            if os.path.splitext(f)[1] not in ['.fastq', '.gz']:
                continue
            # skip R2 files, since we only want filenames with R1 suffix
            elif os.path.basename(f).split(".")[0].split("_")[-1] == "R2":
                continue
            else: 
                #derive from filename
                input_files.append(os.path.basename(f).split(".")[0])
                if not os.path.exists("input_readcount.pkl"):
                    input_file_readcount[f] = count_reads(os.path.join(root, f))

    print "Files:", set(input_files)
    if os.path.exists("input_readcount.pkl"):
        input_file_readcount = pickle.load("input_readcount.pkl")
    else:
        pickle.dump(input_file_readcount, "input_readcount.pkl")

    missing_files = []
    bad_read_count_files = []
    # check if files are in output directory
    for root, dirs, files in os.walk(output_dir):
        # for each input file, make sure corresponding output file is present
        for f in input_files:
            for k in kmer_sizes:
                ec_file = "_".join([tool_name, f, str(k)])
                ec_file = ".".join([ec_file, "corrected", "fastq", "gz"])
                ec_file = os.path.join(output_dir, ec_file)
                print ec_file
                # if file does not exist, log and output warning
                if not os.path.exists( ec_file ):
                    missing_files.append(ec_file)
                #else:
                    #raw_count = count_reads(
    with open(missing_files_log, 'w') as missing_log:
        missing_log.write("======== MISSING FILES: ========\n")
        for missing_f in missing_files:
            missing_log.write(missing_f + "\n")


if __name__ == "__main__":
    main()
