#!/home/mlborowiec/anaconda/bin/python
# -*- coding: utf-8 -*-
# modified from Phyluce's phyluce_assembly_get_fasta_lengths

# This script prints custom UPP (Nguyen et al. 2015 Genome Biol 16:124)
# command based on the distribution of sequence lengths in unaligned fasta

from __future__ import division
import argparse
import heapq
import os
import gzip
import numpy
import re
from itertools import groupby
from phyluce.helpers import FullPaths, is_file

def get_args():
    """Get arguments from CLI"""
    parser = argparse.ArgumentParser(
            description="""Print UPP aligner commands for a directory of FASTA files""")
    parser.add_argument(
            "--input",
            required=True,
            type=is_file,
            action=FullPaths,
            help="""The input fasta file"""
        )
    parser.add_argument(
            "--backb_prop",
            type=float,
            default=0.3,
            help="""Proprtion of taxa to be used in backbone UPP alignment"""
        )
    return parser.parse_args()

def fasta_iter(fasta):
    """modified from @brent_p on stackoverflow.  yield tuple of header, sequence"""
    if fasta.endswith('.gz'):
        with gzip.open(fasta) as f:
            faiter = (x[1] for x in groupby(f, lambda line: line[0] == ">"))
            for header in faiter:
                yield sum(len(s.strip()) for s in faiter.next())
    else:
        with open(fasta) as f:
            faiter = (x[1] for x in groupby(f, lambda line: line[0] == ">"))
            for header in faiter:
                yield sum(len(s.strip()) for s in faiter.next())

def main():
    args = get_args()
    lengths = numpy.array([int(record) for record in fasta_iter(args.input)])

    #get the no. of taxa as a fraction of all seqs
    prop_no_taxa = int(round((len(lengths) * args.backb_prop), 0)) ## for uce-10004 this is 32 taxa
    
    # use that proportion to get longest sequences
    n_longest = heapq.nlargest(prop_no_taxa, lengths)
    
    # determine the threshold from the longest that
    # would capture the desired no of sequences
    thresh = round(((max(n_longest) - min(n_longest)) / max(n_longest)), 5) ## for uce-10004 this is 0.6301

    # helper stuff
    upp_loc = 'run_upp.py'
    in_fn = args.input
    m = re.search('uce-[0-9]+', in_fn)
    locus = m.group()
    len_str = ','.join(map(str, lengths))
    # print the command
    print('{0} -s ./{1}/{1}.unaligned.fasta -M {3} -T {2} -d ./{1} -o output_{1} -p /datadrive/sepp_tmp -l 5'.format(upp_loc, locus, thresh, max(n_longest)))

if __name__ == '__main__':
    main()
