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
    parser.add_argument('-i', '--input-file', dest='input_file', required=True, help='Input (original) fastq file')
    parser.add_argument('-o', '--output-file', dest='output_file', required=True, help='Output (corrected) fastq file')
    args = parser.parse_args()

    input_file = args.input_file
    output_file = args.output_file

    # make sure files exist
    if not os.path.exists(input_file):
        print "Error", input_file, "does not exist"
        exit()

    if not os.path.exists(output_file):
        print "Error", output_file, "does not exists"
        exit()

    # count number of reads and compare
    num_input_reads = count_reads(input_file)
    num_output_reads = count_reads(output_file)

    # 2*num_input_reads because num_output_reads should 
    # contain both paired end inputs
    if not (2*num_input_reads) == num_output_reads:
        print "ERROR: number of input reads no not match output reads!", 2*num_input_reads, " vs. ", num_output_reads
        exit(1)
    else:
        print "Success"
    exit()  


if __name__ == "__main__":
    main()
