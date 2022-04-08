#! /usr/bin/env python3
# coding: utf-8

import sys, getopt
import csv
import io
import logging
import math
import os
import json
import re
from Bio import SeqIO

def main(argv):
    cpus=1
    try:
        opts, args = getopt.getopt(argv,"hv:c:d:r:o:")
    except getopt.GetoptError:
        print('specific_mutations.py -v <virus properties file from nextclade data> -c <clade1> -d <clade2> -r <reference> -o <output file>')
        sys.exit(1)
    for opt, arg in opts:
        if opt == '-h':
            print('specific_mutations.py -v <virus properties file from nextclade data> -c <clade1> -d <clade2> -r <reference> -o <output file>')
            sys.exit(0)
        elif opt in ("-v"):
            virusprop = arg
        elif opt in ("-c"):
            clade1 = arg
        elif opt in ("-d"):
            clade2 = arg
        elif opt in ("-r"):
            reffile = arg
        elif opt in ("-o"):
            outfile = arg

    offset = 54
    logging.info(f'Virus prop file : {virusprop}')
    logging.info(f'Clade 1  : {clade1}')
    logging.info(f'Clade 2  : {clade2}')
    logging.info(f'Output file : {outfile}')

    # Import virus properties from nextclade
    lineages = dict()
    with open(virusprop) as file:
        lineages=json.load(file)

    with open(reffile,'rt') as file:
        for record in SeqIO.parse(file, "fasta"):
            refseq=record.seq


    lcount1 = lineages["nucMutLabelMapReverse"].get(clade1,None)
    lcount2 = lineages["nucMutLabelMapReverse"].get(clade2,None)

    pos1= dict()
    pos2= dict()

    for p in lcount1:
        ptrim = p[:-1]
        pos1[ptrim] = p[-1]
    for p in lcount2:
        ptrim = p[:-1]
        pos2[ptrim] = p[-1]


    with open(outfile, 'w') as ofile:
        print(f'Pos\t{clade1}\t{clade2}',file=ofile)
        for m in lcount2:
            ptrim = m[:-1]
            alt=m[-1]
            if(m not in lcount1):
                if ptrim in pos1:
                    ref=pos1[ptrim]
                else:
                    ref=refseq[int(ptrim)-1]
                print(f'{ptrim}\t{ref}\t{alt}',file=ofile)

        for m in lcount1:
            ptrim = m[:-1]
            alt=m[-1]
            if(m not in lcount2):
                if ptrim in pos2:
                    ref=pos2[ptrim]
                else:
                    ref=refseq[int(ptrim)-1]
                print(f'{ptrim}\t{alt}\t{ref}',file=ofile)

if __name__ == "__main__":
    logging.basicConfig(level=os.environ.get("LOGLEVEL", "INFO"))
    main(sys.argv[1:])
