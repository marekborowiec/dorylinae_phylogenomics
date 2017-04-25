This is code used in a phylogenomic study of the ant subfamily Dorylinae. Please consider citing it if you are using the scripts:

@@@

Structure:

```
.
├── aln_bin_by_rate.py 	# concatenate loci by rate into n bins
├── aln_choose.py 	# sample n loci across rate range
├── chronos_dating.R 	# Chronos divergence dating
├── compo_tests_corrected.py 	# compositional heterogeneity tests in p4
├── compo_tests_wrapper.py 	# wrapper for tests on multiple loci
├── compute_rcfv.py 	# compute RCFV (compositional heterogeneity)
├── print_upp_command.py 	# print UPP alignment commands
├── significance_tests.R 	# testing locus properties in matrices
├── tree_props.R 	# avg. bootstrap, branch length, saturation, plotting gene trees
└── upp_script.sh 	# UPP commands used for alignment of each locus
```