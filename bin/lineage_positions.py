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
        opts, args = getopt.getopt(argv,"hv:o:")
    except getopt.GetoptError:
        print('lineage_positions.py -v <virus properties file from nextclade data> -o <suffix of output file>')
        sys.exit(1)
    for opt, arg in opts:
        if opt == '-h':
            print('lineage_positions.py -v <virus properties file from nextclade data> -o <suffix of output file>')
            sys.exit(0)
        elif opt in ("-v"):
            virusprop = arg
        elif opt in ("-o"):
            outfile = arg

    logging.info(f'Virus prop file : {virusprop}')
    logging.info(f'Output file : {outfile}')

    # Import virus properties from nextclade
    lineages = dict()
    with open(virusprop) as file:
        lineages=json.load(file)

    for lin in lineages["nucMutLabelMapReverse"]:
        with open(f'{lin}_{outfile}', 'w') as ofile:
            print(f'Pos\t{lin}',file=ofile)
            for p in lineages["nucMutLabelMapReverse"][lin]:
                pos = p[:-1]
                alt = p[-1]
                print(f'{pos}\t{alt}',file=ofile)

if __name__ == "__main__":
    logging.basicConfig(level=os.environ.get("LOGLEVEL", "INFO"))
    main(sys.argv[1:])
