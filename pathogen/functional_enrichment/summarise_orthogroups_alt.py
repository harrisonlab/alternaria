#!/usr/bin/python

'''
Summarise ortholog groups by organism
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

ap.add_argument('--orthogroups',required=True,type=str,help='A fasta file of the assembled contigs')
ap.add_argument('--orthomclIDs',required=True,type=str,nargs='+',help='All isolates used in orthology analysis')
ap.add_argument('--species',required=True,type=str,nargs='+',help='Species designations for each organism in the orthomclID list')
ap.add_argument('--prefix',required=True,type=str,help='Outfile prefix')


conf = ap.parse_args()
# grp1_isolates = conf.set1
# grp2_isolates = conf.set2
all_isolates = conf.orthomclIDs
species_list = conf.species
o = conf.prefix

with open(conf.orthogroups) as f:
    orthogroup_lines = f.readlines()

orthogroup_content_dict = defaultdict(list)
species_counts_dict = defaultdict(list)
orthogroup_content_dict["OrthogroupID"].extend(["OrthomclID", "Species"])
for isolate, species in zip(all_isolates, species_list):
    orthogroup_content_dict[isolate].extend([isolate, species])

enrichment_lines = []
for line in orthogroup_lines:
    line = line.rstrip("\n")
    split_line = line.split(" ")
    orthogroup_id = split_line[0].replace(":", "")
    orthogroup_content_dict["OrthogroupID"].append(orthogroup_id)
    species_counts_dict.clear()
    counts_list = []
    for isolate, species in zip(all_isolates, species_list):
        num_genes = line.count((isolate + "|"))
        orthogroup_content_dict[isolate].append(str(num_genes))
        species_counts_dict[species].append(num_genes)
        counts_list.append(str(num_genes))
    enrichment = ''
    if min(species_counts_dict['A.alternata_ssp._arborescens']) > max(species_counts_dict['A.alternata_ssp._tenuissima']):
        enrichment = 'arborescens expanded'
    elif min(species_counts_dict['A.alternata_ssp._tenuissima']) > max(species_counts_dict['A.alternata_ssp._arborescens']):
        enrichment = 'tenuissima expanded'
    enrichment_lines.append("\t".join([orthogroup_id, "\t".join(counts_list), enrichment]))

outfile = open(o + "_expansion_status.tsv", 'w')
outfile.write("\t".join(["orthogroupID", "\t".join(all_isolates), "expansion_status"]) + "\n")
outfile.write("\n".join(enrichment_lines))

outfile = open(o + ".tsv", 'w')
outfile.write("\t".join(orthogroup_content_dict["OrthogroupID"])+ "\n")
for isolate in all_isolates:
    outfile.write("\t".join(orthogroup_content_dict[isolate]) + "\n")
