#!/usr/bin/python

'''

'''


#-----------------------------------------------------
# Step 1
# Import variables
#-----------------------------------------------------

import sys
import argparse
import re
from sets import Set
from collections import defaultdict
from operator import itemgetter

ap = argparse.ArgumentParser()

ap.add_argument('--isolate_list_A',required=True,nargs='+',type=str,help='List of isolates in species 1')
ap.add_argument('--interpro_list_A',required=True,nargs='+',type=str,help='List of IPR files for species 1')

ap.add_argument('--isolate_list_B',required=True,nargs='+',type=str,help='List of isolates in species 2')
ap.add_argument('--interpro_list_B',required=True,nargs='+',type=str,help='List of IPR files for species 2')

conf = ap.parse_args()


list_A = conf.isolate_list_A
list_B = conf.isolate_list_B

# print "Isolates listed"
# print "group A:"
# print "\t".join(list_A)
# print "group B:"
# print "\t".join(list_B)


#-----------------------------------------------------
# Step 2
# Define classes
#-----------------------------------------------------

class Annot_obj(object):
    """A gene identified as differentially expressed in one isolate.
    Attributes:
        gene_name: A string representing the gene name.
        conditions_tested: List of conditions that the gene was tested for DE.
        conditions_positive: A binary list of integers representing whether the
            gene tested positive for each condition listed in conditions tested.
    """

    def __init__(self, ID, desc):
        """Return a Annot_obj whose name is *gene_name*"""
        self.ipr_ID = ID
        self.ipr_desc = desc
        self.presence_dict = defaultdict(set)
        self.genes_by_isolate = []

    def add_interpro(self, isolate, gene_id):
        """Add InterPro information"""
        self.presence_dict[isolate].add(gene_id)

    def summarise(self, isolate_list):
        """"""
        for key in isolate_list:
            gene_count = len(self.presence_dict[key])
            self.genes_by_isolate.append(str(gene_count))
        return "\t".join(self.genes_by_isolate)


#-----------------------------------------------------
# Step 3
# Build interproscan dictionary
#-----------------------------------------------------

IPR_obj_dict = defaultdict()

for isolate, ipr_file in zip(list_A + list_B, conf.interpro_list_A + conf.interpro_list_B):
    # print ("\t".join([isolate, ipr_file]))
    with open(ipr_file) as f:
        ipr_lines = f.readlines()
    for line in ipr_lines:
        if not 'IPR' in line:
            continue
        line = line.rstrip()
        split_line = line.split('\t')
        gene_id = split_line[0]
        IPR_ID = split_line[11]
        IPR_desc = split_line[12]
        # print "\t".join([gene_id, IPR_ID])
        if IPR_ID in IPR_obj_dict:
            IPR_obj_dict[IPR_ID].add_interpro(isolate, gene_id)
        else:
            IPR_obj = Annot_obj(IPR_ID, IPR_desc)
            IPR_obj.add_interpro(isolate, gene_id)
            IPR_obj_dict[IPR_ID] = IPR_obj

print "\t".join(['IPR_ID', 'IPR_description', "\t".join(list_A + list_B)])

for key in IPR_obj_dict.keys():
    IPR_obj = IPR_obj_dict[key]
    genes_by_isolate = IPR_obj.summarise(list_A + list_B)
    print "\t".join([IPR_obj.ipr_ID, IPR_obj.ipr_desc, genes_by_isolate])
