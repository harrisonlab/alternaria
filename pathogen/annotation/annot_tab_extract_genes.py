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

ap = argparse.ArgumentParser()

ap.add_argument('--annot',required=True,type=str,help='Annotation table')
ap.add_argument('--genes',required=True,type=str,help='list of genes to extract')

conf = ap.parse_args()


with open(conf.annot) as f:
    annot_lines = f.readlines()

with open(conf.genes) as f:
    gene_lines = f.readlines()

print annot_lines[0].rstrip()
for gene in gene_lines:
    gene = gene.rstrip()
    for line in annot_lines[1:]:
        line = line.rstrip()
        # split_line = line.split("\t")
        ID = line.split(".")[0]
        if gene == ID:
            print line
