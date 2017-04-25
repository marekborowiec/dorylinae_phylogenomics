#! /usr/bin/env python3

# compute_rcfv.py by Marek Borowiec
# RCFV calculation on parsed alignment

from amas import AMAS
from collections import Counter
from sys import argv

# in file name
in_fn = argv[1]
# print per taxon RCFVs: yes/no
per_taxon = argv[2]

# use AMAS to get alignment parsed as a dictionary { taxon : sequence }
meta_aln = AMAS.MetaAlignment(in_files=[in_fn], data_type="dna",in_format="phylip", cores=1)
parsed_aln = meta_aln.get_parsed_alignments()
aln = parsed_aln[0]

# determine what letter you will be looking at; here just ACGT
alphabet = ["A","C","G","T","K","M","R","Y","S","W","B","V","H","D","X", \
     "N", "O", "-","?"]
missing_ambiguous_chars = ["K","M","R","Y","S","W","B","V","H","D","X", \
     "N", "O", "-","?"]
used_alphabet = list(Counter(alphabet) - Counter(missing_ambiguous_chars))

# For each character of interest (ACGT etc.) determine the overall frequency across taxa
def get_ov_count(char, aln_dict):
    char_count = sum(seq.count(char) for seq in aln_dict.values())
    return char_count

overall_count = sum(get_ov_count(char, aln) for char in used_alphabet)

def get_ov_freq(char, aln_dict):
    char_freq = (get_ov_count(char, aln_dict)) / overall_count
    return char_freq

ov_freqs = {char : get_ov_freq(char, aln) for char in used_alphabet}

# For each taxon determine the frequency of each character of interest
differences = []
per_taxon_differences = []
for taxon, seq in aln.items():
    taxon_diffrs = []
    for char, freq in ov_freqs.items():
        taxon_freq = seq.count(char) / sum(seq.count(char) for char in used_alphabet)
        # For each character subtract taxon-specific frequency from overall frequencies
        difference = abs(freq - taxon_freq)
        record = (taxon, char, difference)
        differences.append(record)
        taxon_diffrs.append(difference)
    per_taxon_differences.append({taxon : (sum(taxon_diffrs) / len(aln.keys()))})

# Sum all the differences and divide by no. of taxa
RCFV = sum(diffr for taxon, char, diffr in differences) / len(aln.keys())

print('{}\t{:.12f}'.format(in_fn, RCFV))

# print per taxon rcfvs if second argument was 'yes'
if per_taxon == 'yes':
    for dff in per_taxon_differences:
        for tax, rcfv in sorted(dff.items()):
            print('{}\t{}\t{:.12f}'.format(in_fn, tax, rcfv))
