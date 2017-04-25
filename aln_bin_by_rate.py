#! /usr/bin/env python3

# This script takes tab-delimited AMAS summary and sorts loci by average 
# branch length (or other statistic of choice)
# It then classifies them into bins of approximately same number of sites
# using the information on each locus' length.

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
        avg_br = float(line[36])
        tpl = (aln_name, aln_length, avg_br)
        return tuple(tpl)

    # create a list of tuples 
    tuple_list = [get_line_data(line) for line in reader]
    # sort by tuple index number
    avg_br_sorted = (sorted(tuple_list,key=itemgetter(2)))
    # get total number of sites
    all_sites = sum(length for fn, length, avg_br in avg_br_sorted)
    #print(avg_br_sorted)
    # define number of bins
    bins = 5
    # get bin length
    bin_length = all_sites / bins

    # function to increment numbers
    def increment(abin, num):
        abin += abin * num
        return abin

    # create a list with upper bounds for each bin
    all_bins = [increment(bin_length, n) for n in range(bins)]
    # create a counter to track no of sites added
    total = 0
    # create empty list of lists for each bin
    alns = [[] for i in range(bins)]

    for tpl in avg_br_sorted:
        # for each tuple check what bin it fits in
        for index, abin_len in enumerate(all_bins):
            if index == 0:
                if total < abin_len:
                    alns[0].append(tpl)
                    total += tpl[1] # alignment length
            else:
                # starting with bin no 2 check if an alignment has already been added
                if tpl not in alns[index-1] and total > all_bins[index-1] and total < abin_len:
                    alns[index].append(tpl)
                    total += tpl[1] # alignment length

    for index, l in enumerate(alns):
        for aln in l:
            print('{}\t{}'.format(aln[0], index+1))
