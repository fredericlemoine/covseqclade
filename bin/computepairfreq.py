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

# This script read the output of bam-readcount, and a coinfection file with one line per genomic position that differs between the two strains and the following columns:
# 1 - position (0-based, on the mapped reference basis)
# 2 - the expected nucleotide in the first strain
# 3 - the expected nucleotide in the second strain
# the header line contains:
# 1 - "Position"
# 2 - Name of the first co-infecting strain
# 3 - Name of the second co-infecting strain
# The output is a file with a header:
# 1 - "Position"
# 2 - Name of the first co-infecting strain
# 3 - Name of the second co-infecting strain
# 4 - "Total"
# A one line per genomic position that differs between the 2 strains and the following columns:
# 1 - the position (0-based, ont the mapped reference basis)
# 2 - the number of read mapping on first strain nucleotide
# 3 - the number of read mapping on the second strain nucleotide
# 4 - the total number of reads passing quality thresholds
def main(argv):
    cpus=1
    try:
        opts, args = getopt.getopt(argv,"hi:c:o:")
    except getopt.GetoptError:
        print('computefreq.py -i <bam-readcount output> -c <co-infection desc file> -o <output file>')
        sys.exit(1)
    for opt, arg in opts:
        if opt == '-h':
            print('computefreq.py -i <bam-readcount output> -c <co-infection desc file> -o <output file>')
            sys.exit(0)
        elif opt in ("-i"):
            countfile = arg
        elif opt in ("-c"):
            coinffile = arg
        elif opt in ("-o"):
            outfile = arg

    logging.info(f'Input readcount file : {countfile}')
    logging.info(f'Input co-inf file file : {coinffile}')
    logging.info(f'Output file : {outfile}')

    # List context seq ID
    positions = dict()
    s1=""
    s2=""
    with open(coinffile) as file:
        tsv_file = csv.reader(file, delimiter="\t")
        header = next(tsv_file)
        if header != None:
            s1=header[1]
            s2=header[2]
            for line in tsv_file:
                pos=int(line[0])
                n1=line[1]
                n2=line[2]
                positions[pos] = {'s1' : n1, 's2': n2 }

    with open(outfile, 'w') as ofile:
        print(f'Position\t{s1}\t{s2}\tTotal',file=ofile)
        # MN908947_Ampliseq
        # 2
        #A
        #202
        #=:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00
        # A:202:59.69:33.64:10.40:106:96:0.00:0.00:0.00:0:0.00:0.00:0.00
        #C:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00
        #G:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00
        #T:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00
        #N:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00
        gaps=[0]*40000
        with open(countfile) as file:
            tsv_file = csv.reader(file, delimiter="\t")
            for line in tsv_file:
                pos=int(line[1])
                counts1=0
                counts2=0
                total=int(line[3])+gaps[pos-1]
                for col in line[4:]:
                    vals=col.split(":")
                    nuc=vals[0]
                    count=int(vals[1])
                    if nuc.startswith("-"):
                        for i in range(0,(len(nuc) - 1)):
                            gaps[pos - 1 + i]+=int(count)
                    if (pos) in positions:
                        if nuc==positions[pos]["s1"]:
                            counts1=count
                        elif nuc==positions[pos]["s2"]:
                            counts2=count
                if (pos) in positions:
                    if gaps[pos-1]>0 and "-"==positions[pos]["s1"]:
                        counts1=gaps[pos-1]
                    if gaps[pos-1]>0 and "-"==positions[pos]["s2"]:
                        counts2=gaps[pos-1]
                    print(f'{pos}\t{counts1}\t{counts2}\t{total}',file=ofile)

if __name__ == "__main__":
    logging.basicConfig(level=os.environ.get("LOGLEVEL", "INFO"))
    main(sys.argv[1:])
