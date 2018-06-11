#!/usr/bin/python

'''

'''


#-----------------------------------------------------
# Step 1
# Import variables & load input files
#-----------------------------------------------------

import sys
import argparse
import re
from sets import Set
from collections import defaultdict
from operator import itemgetter

ap = argparse.ArgumentParser()

ap.add_argument('--subset',required=True,type=str,help='list of a subset of genes to extract reciprical hits from')
ap.add_argument('--hits_a_vs_b',required=True,type=str,help='Blast hits in format 6 of organism A vs organism B')
ap.add_argument('--hits_b_vs_a',required=True,type=str,help='Blast hits in format 6 of organism B vs organism A')

conf = ap.parse_args()

with open(conf.subset) as f:
    subset_lines = f.readlines()

with open(conf.hits_a_vs_b) as f:
    hits_a_lines = f.readlines()

with open(conf.hits_b_vs_a) as f:
    hits_b_lines = f.readlines()

#----
#
#----

hits_dict_A = defaultdict(str)
for line in hits_a_lines:
    line = line.rstrip()
    split_line = line.split("\t")
    hits_dict_A[split_line[0]] = split_line[1]

hits_dict_B = defaultdict(str)
for line in hits_b_lines:
    line = line.rstrip()
    split_line = line.split("\t")
    hits_dict_B[split_line[0]] = split_line[1]

#----
#
#----

for line in subset_lines:
    gene = line.rstrip()
    hit = hits_dict_A[gene]
    recipricol_hit = hits_dict_B[hit]
    # print("\t".join([gene, hit, recipricol_hit]))
    if gene == recipricol_hit:
        print("\t".join([gene, hit]))
    else:
        print("\t".join([gene, '']))
