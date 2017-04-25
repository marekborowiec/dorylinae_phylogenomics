#! /usr/bin/env python3

# This script takes tab-delimited AMAS summary and sorts loci by
# average branch length. It then slices the list of sorted alignments
# to get desired number of equally distributed loci.
# This code was written to subsample loci for time-calibrated analyses.

from math import ceil
from sys import argv
from operator import itemgetter
import csv

# input is the argument to script
in_fn = argv[1]

with open(in_fn, 'r') as f:
    # skip header
    next(f)

    reader = csv.reader(f,delimiter='\t')
    # define columns of interest and return as tuple
    def get_line_data(line):
        aln_name = line[0]
        aln_length = int(line[2])
        #prop_pars_inf = float(line[9])
        avg_br = float(line[37])
        tpl = (aln_name, aln_length, avg_br)
        return tuple(tpl)
    # create a list of tuples 
    tuple_list = [get_line_data(line) for line in reader]
    # sort by tuple index number
    avg_br_sorted = (sorted(tuple_list,key=itemgetter(2)))
    # get total number of sites
    all_sites = sum(length for fn, length, avg_br in avg_br_sorted)
    selected = avg_br_sorted[0::int(ceil(len(avg_br_sorted) / 100))]
    #print(len(selected))
    #print(selected)
    for uce, length, avg_br in selected:
        print('cp {}_alignment_masked.fasta TimeTree100'.format(uce))
