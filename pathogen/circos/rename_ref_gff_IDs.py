#!/usr/bin/python

'''
'''

#-----------------------------------------------------
# Step 1
# Import variables & load input files
#-----------------------------------------------------

import sys
import argparse

ap = argparse.ArgumentParser()
ap.add_argument('--gff',required=True,type=str,help='Reference Gff file')

conf = ap.parse_args()

with open(conf.gff) as f:
    gff_lines = f.readlines()

for line in gff_lines:
    line = line.rstrip()
    if line.startswith('#'):
        print line
        continue
    split_line = line.split()
    if any(x in split_line[2] for x in ['mRNA', 'transcript']):
        # print line
        col9 = split_line[8]
        split_col9 = col9.split(';')
        ID = split_col9[0]
        prot_id = split_col9[4].split('|')[-1]
        prot_id = prot_id.replace('mRNA.','')
        line = line.replace(ID, 'ID=' + prot_id)
    print line
