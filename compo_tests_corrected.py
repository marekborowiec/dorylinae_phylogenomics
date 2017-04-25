#! /usr/bin/env python

# compo_tests_corrected.py by Marek Borowiec

# This script takes the input alignment file name 
# and number of simulations as its arguments
# to perform Chi-squared and corrected tests
# of compositional heterogeneity using p4 (Foster 2004)

from sys import argv
import p4

in_file = argv[1]
n = argv[2]

# don't check for empty sequences or sites since
# p4 does not consider those in the test anyways
p4.var.doCheckForAllGapColumns = False

print(
'''

========== calculating test stats for {} ========== 

'''.format(in_file)
)

p4.read(in_file)

a = p4.var.alignments[0]
dm = a.pDistances()
t = dm.bionj()
d = p4.Data()
t.data = d
t.newComp(free=1, spec='empirical')
t.newRMatrix(free=1, spec='ones')
t.setNGammaCat(nGammaCat=4)
t.newGdasrv(free=1, val=0.5)
t.setPInvar(free=0, val=0.0)
t.optLogLike()
t.name = 'homogOpt'
t.tPickle()


# Then, do the test ...
p4.read('homogOpt.p4_tPickle')
t = p4.var.trees[0]
t.data = d
t.compoTestUsingSimulations(nSims=int(n), doChiSquare=True)

