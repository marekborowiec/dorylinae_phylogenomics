#! /usr/bin/env python

# compo_tests_wrapper.py by Marek Borowiec

# This is just a wrapper for compo_tests_corrected.py

import glob, subprocess as sp

in_files = glob.glob('*fasta')

no_simulations = 200

for file_name in in_files:
    print(file_name)
    command = 'python /datadrive/p4-phylogenetics/compo_tests_corrected.py {} {}'.format(file_name, no_simulations)
    sp.call(command, shell=True)

