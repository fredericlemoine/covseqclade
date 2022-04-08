#! /usr/bin/env python3
# coding: utf-8

import sys, getopt
import csv
import io
import tarfile
import logging
import os
import re
import subprocess

# This script read the output of bam-readcount, and a
# mutation/clade definition file, tab separated: 
# 1 - position (0-based, on the mapped reference basis)
# 2 - the expected nucleotide in the clade, at this position
# the header line contains:
# 1 - "Position"
# 2 - Name of the clade
# The output is a file with a header:
# 1 - "Position"
# 2 - Name of the clade
# 3 - "Total"
# A one line per genomic position that differs between the 2 clades and the following columns:
# 1 - the position (0-based, ont the mapped reference basis)
# 2 - the number of read mapping on the nucleotide specific to the clade
# 3 - the total number of reads passing quality thresholds
def main(argv):
    cpus=1
    offset=0
    try:
        opts, args = getopt.getopt(argv,"hi:c:o:d:")
    except getopt.GetoptError:
        print('computefreq.py -i <bam-readcount output> -c <co-infection desc file> -o <output file> -d <position offset>')
        sys.exit(1)
    for opt, arg in opts:
        if opt == '-h':
            print('computefreq.py -i <bam-readcount output> -c <co-infection desc file> -o <output file> -d <position offset')
            sys.exit(0)
        elif opt in ("-i"):
            countfile = arg
        elif opt in ("-c"):
            coinffile = arg
        elif opt in ("-o"):
            outfile = arg
        elif opt in ("-d"):
            offset = int(arg)

    logging.info(f'Input readcount file : {countfile}')
    logging.info(f'Input co-inf file file : {coinffile}')
    logging.info(f'Output file : {outfile}')

    # List context seq ID
    positions = dict()
    lin=""
    with open(coinffile) as file:
        tsv_file = csv.reader(file, delimiter="\t")
        header = next(tsv_file)
        if header != None:
            lin=header[1]
            for line in tsv_file:
                pos=int(line[0])
                n=line[1]
                positions[pos] = n

    with open(outfile, 'w') as ofile:
        print(f'Position\tCount\tLineage\tTotal',file=ofile)            
        covered = dict()
        with open(countfile) as file:
            tsv_file = csv.reader(file, delimiter="\t")
            for line in tsv_file:
                pos=int(line[1])
                pos+=offset
                counts=0
                total=int(line[3])
                for col in line[4:]:
                    vals=col.split(":")
                    nuc=vals[0]
                    count=int(vals[1])
                    if (pos) in positions:
                        if nuc==positions[pos]:
                            counts=count
                if (pos) in positions:
                    covered[pos]=1
                    print(f'{pos}\t{counts}\t{lin}\t{total}',file=ofile)
            for pos in positions:
                if pos not in covered:
                    print(f'{pos}\t0\t{lin}\t0',file=ofile)

        
if __name__ == "__main__":
    logging.basicConfig(level=os.environ.get("LOGLEVEL", "INFO"))
    main(sys.argv[1:])
