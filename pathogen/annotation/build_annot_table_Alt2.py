#!/usr/bin/python

'''
This program is used to build information on all the genes predicted in
an annotated genome. These commands take information on location of genes
& suppliment this information with information on interproscan domains
and swissprot annotations.
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

# ap.add_argument('--genome',required=True,type=str,help='A fasta file of the assembled contigs')
ap.add_argument('--genes_gff',required=True,type=str,help='A gff file of the genes')
ap.add_argument('--SigP',required=True,type=str,help='A file containing a list of signal-peptide containing genes')
ap.add_argument('--TM_list',required=True,type=str,help='txt file of headers from gene testing positive for tranmembrane proteins by TMHMM')
# # ap.add_argument('--GPI_list',required=True,type=str,help='txt file of headers from gene testing positive for GPI anchors as identified by GPI-SOM')
ap.add_argument('--EffP_list',required=True,type=str,help='A file containing results of effectorP')
ap.add_argument('--CAZY_list',required=True,type=str,help='A file containing results of CAZY')
ap.add_argument('--PhiHits',required=True,type=str,help='BLAST results of Phibase proteins vs predicted gene models')
ap.add_argument('--ToxinHits',required=True,type=str,help='BLAST results of CDC toxins vs predicted gene models')
ap.add_argument('--InterPro',required=True,type=str,help='The Interproscan functional annotation .tsv file')
ap.add_argument('--Swissprot',required=True,type=str,help='A parsed table of BLAST results against the Swissprot database. Note - must have been parsed with swissprot_parser.py')
ap.add_argument('--Antismash',required=True,type=str,help='Output of Antismash parsed into a tsv file of gene names, contig, secmet function and cluster ID')
# ap.add_argument('--Smurf',required=True,type=str,help='Output of Smurf parsed into a tsv file of gene names, contig, secmet function and cluster ID')
ap.add_argument('--TFs',required=True,type=str,help='Tab seperated of putative transcription factors and their domains as identified by interpro2TFs.py')
ap.add_argument('--orthogroups', required=True,type=str,help='A file containing results of orthology analysis')
ap.add_argument('--strain_id',required=True,type=str,help='The identifier of this strain as used in the orthology analysis')
ap.add_argument('--OrthoMCL_all',required=True,type=str,nargs='+',help='The identifiers of all strains used in the orthology analysis')

conf = ap.parse_args()



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

    def __init__(self):
        """Return a Annot_obj whose name is *gene_name*"""
        self.gene_name = ''
        self.contig = ''
        self.start = ''
        self.stop = ''
        self.strand = ''
        self.sigp = ''
        self.transmem = ''
        self.secreted = ''
        self.effp = ''
        self.cazy = ''
        self.interpro = set()
        self.swissprot = ''
        self.phi = ''
        self.toxin = ''
        self.antismash = ''
        self.smurf = ''
        self.cluster_id = ''
        self.orthogroup_id = ''
        self.content_counts = ''
        self.content_str = ''
        self.transcriptionfactor = set()

    def set_conditions(self, gff_elements):
        """Reset conditions_tested to a list of conditions"""
        self.contig = gff_elements[0]
        self.start = gff_elements[3]
        self.stop = gff_elements[4]
        self.strand = gff_elements[6]
        gene_features = gff_elements[8].split(';')
        gene_id = gene_features[0]
        gene_id = gene_id.replace('ID=', '')
        self.gene_name = gene_id

    def add_sigp(self):
        """Add SignalP information"""
        self.sigp = 'Yes'
        self.secreted = 'Yes'
    def add_transmem(self):
        """Add transmembrane protein information"""
        # print self.gene_name
        if self.secreted == 'Yes':
            self.transmem = 'Yes'
            self.secreted = ''
    def add_effp(self):
        """Add EffectorP information"""
        self.effp = 'Yes'
    def add_cazy(self, line):
        """Add CAZY information"""
        split_line = line.split()
        self.cazy = 'CAZY:' + split_line[0].replace('.hmm', '')
    def add_antismash(self, line, num):
        split_line = line.split()
        function = split_line[2]
        cluster_name = "AS_" + str(num)
        self.antismash = ":".join([cluster_name, function])
    def add_smurf(self, line):
        split_line = line.split()
        function = split_line[2]
        cluster_name = "SM" + split_line[3].replace('Cluster_', '')
        self.smurf = ":".join([cluster_name, function])
    def add_cluster_id(self, i):
        cluster_id = "Cluster_" + str(i)
        evidence = []
        if '' != self.antismash:
            evidence.append(self.antismash)
        if '' != self.smurf:
            evidence.append(self.smurf)
        self.cluster_id = ":".join([cluster_id, "|".join(evidence)])

    def add_interpro(self, line):
        """Add InterPro information"""
        split_line = line.split("\t")
        interpro_columns = []
        index_list = [4, 5, 11, 12]
        for x in index_list:
            if len(split_line) > x:
                interpro_columns.append(split_line[x])
        ipr_line = ";".join(interpro_columns)
        self.interpro.add(ipr_line)
    def add_swissprot(self, line):
        """Add swissprot information"""
        split_line = line.split("\t")
        self.swissprot = "|".join(itemgetter(14, 12, 13)(split_line))
        # self.cazy = 'CAZY:' + split_line[0].replace('.hmm', '')
    def add_phi(self, line):
        line = line.rstrip()
        split_line = line.split()
        self.phi = split_line[1]
    def add_toxin(self, line):
        line = line.rstrip()
        split_line = line.split()
        self.toxin = split_line[1]
    def add_orthogroup(self, orthogroup_id, content_counts, content_str):
        """Add swissprot information"""
        self.orthogroup_id = orthogroup_id
        self.content_counts = content_counts
        self.content_str = content_str
    def add_TF(self, line):
        """Add swissprot information"""
        split_line = line.split("\t")
        TF_function = split_line[2]
        self.transcriptionfactor.add(TF_function)

#-----------------------------------------------------
# Step X
# Read input files
#-----------------------------------------------------

# with open(conf.genome) as f:
#     contig_lines = f.readlines()
with open(conf.genes_gff) as f:
    gff_lines = f.readlines()
with open(conf.SigP) as f:
    sigP_lines = f.readlines()
with open(conf.TM_list) as f:
    TM_lines = f.readlines()
with open(conf.EffP_list) as f:
    effP_lines = f.readlines()
with open(conf.TFs) as f:
    TF_lines = f.readlines()
with open(conf.CAZY_list) as f:
    cazy_lines = f.readlines()
with open(conf.Antismash) as f:
    antismash_lines = f.readlines()
# with open(conf.Smurf) as f:
#     smurf_lines = f.readlines()
with open(conf.PhiHits) as f:
    phibase_lines = f.readlines()
with open(conf.ToxinHits) as f:
    toxin_lines = f.readlines()
with open(conf.InterPro) as f:
    InterPro_lines = f.readlines()
with open(conf.Swissprot) as f:
    swissprot_lines = f.readlines()
with open(conf.orthogroups) as f:
    orthogroup_lines = f.readlines()
gene_dict = defaultdict(list)


#-----------------------------------------------------
# Step X
# Create an annotation object for each gene
#-----------------------------------------------------


for line in gff_lines:
    if "gff-version" in line:
        continue
    if line.startswith('#'):
        continue
    line = line.rstrip()
    split_line = line.split("\t")
    if 'mRNA' in split_line[2]:
        gene_features = split_line[8].split(';')
        gene_id = gene_features[0]
        gene_id = gene_id.replace('ID=', '')
        gene_obj = Annot_obj()
        gene_obj.set_conditions(split_line)
        gene_dict[gene_id] = gene_obj


#-----------------------------------------------------
# Step X
# And annotation information to gene objects
#-----------------------------------------------------


for line in sigP_lines:
    line = line.rstrip()
    if line.startswith('>'):
        split_line = line.split()
        gene_id = split_line[0].replace('>', '')
        gene_dict[gene_id].add_sigp()

for line in TM_lines:
    line = line.rstrip()
    gene_id = line.split("\t")[0]
    gene_dict[gene_id].add_transmem()

for line in effP_lines:
    line = line.rstrip()
    gene_id = line
    gene_dict[gene_id].add_effp()

for line in TF_lines:
    line = line.rstrip("\n")
    split_line = line.split("\t")
    gene_id = split_line[0]
    gene_dict[gene_id].add_TF(line)

for line in phibase_lines:
    line = line.rstrip()
    gene_id = line.split("\t")[0]
    gene_dict[gene_id].add_phi(line)

for line in toxin_lines:
    line = line.rstrip()
    gene_id = line.split("\t")[0]
    gene_dict[gene_id].add_toxin(line)

cazy_dict = defaultdict(list)
for line in cazy_lines:
    if line.startswith('#'):
        continue
    line = line.rstrip()
    split_line = line.split()
    gene_id = split_line[2]
    gene_dict[gene_id].add_cazy(line)

cluster_id_set = set()
for line in antismash_lines:
    line = line.rstrip("\n")
    split_line = line.split("\t")
    gene_id = split_line[0]
    cluster_id = split_line[3]
    cluster_id_set.add(cluster_id)
    cluster_num = len(cluster_id_set)
    gene_dict[gene_id].add_antismash(line, cluster_num)

# for line in smurf_lines:
#     line = line.rstrip("\n")
#     split_line = line.split("\t")
#     gene_id = split_line[0]
#     gene_dict[gene_id].add_smurf(line)


for line in swissprot_lines:
    line = line.rstrip("\n")
    split_line = line.split("\t")
    gene_id = split_line[0]
    gene_dict[gene_id].add_swissprot(line)


interpro_set =  Set([])
interpro_dict = defaultdict(list)

for line in InterPro_lines:
    line = line.rstrip()
    split_line = line.split()
    gene_id = split_line[0]
    gene_dict[gene_id].add_interpro(line)

strain_id = conf.strain_id + "|"
all_isolates = conf.OrthoMCL_all
orthogroup_dict = defaultdict(list)
orthogroup_content_dict = defaultdict(list)
for line in orthogroup_lines:
    line = line.rstrip("\n")
    split_line = line.split(" ")
    orthogroup_id = split_line[0].replace(":", "")
    orthogroup_contents = []
    orthogroup_content_dict.clear()
    for isolate in all_isolates:
        num_genes = line.count((isolate + "|"))
        orthogroup_contents.append(str(isolate) + "(" + str(num_genes) + ")")
        content_counts = ":".join(orthogroup_contents)
        orthogroup_content_dict[isolate] = num_genes
    for gene_id in split_line[1:]:
        content_str = ",".join(split_line[1:])
        if strain_id in gene_id:
            gene_id = gene_id.replace(strain_id, '')
            gene_dict[gene_id].add_orthogroup(orthogroup_id, content_counts, content_str)

gene_list = sorted(gene_dict.keys(), key=lambda x:int(x[1:].split('.')[0]))

i = 0
prev_gene = False
for gene_id in gene_list:
    gene_obj = gene_dict[gene_id]
    if any('' != x for x in ([gene_obj.antismash, gene_obj.antismash])):
        secmet_cluster = True
        if prev_gene == False:
            i += 1
        gene_obj.add_cluster_id(i)
        prev_gene = True
    else:
        prev_gene = False


print "\t".join([
    'gene_name',
    'contig',
    'start',
    'stop',
    'strand',
    'sigp',
    'transmem',
    'secreted',
    'effp',
    'cazy',
    'transcriptionfactor',
    'cluster_id',
    'orthogroup_id',
    'content_counts',
    'content_str',
    'phi',
    'toxin',
    'swissprot',
    'interpro'
    ])

for gene_id in gene_list:
    gene_obj = gene_dict[gene_id]
    print "\t".join([
        gene_obj.gene_name,
        gene_obj.contig,
        gene_obj.start,
        gene_obj.stop,
        gene_obj.strand,
        gene_obj.sigp,
        gene_obj.transmem,
        gene_obj.secreted,
        gene_obj.effp,
        gene_obj.cazy,
        ";".join(gene_obj.transcriptionfactor),
        gene_obj.cluster_id,
        gene_obj.orthogroup_id,
        gene_obj.content_counts,
        gene_obj.content_str,
        gene_obj.phi,
        gene_obj.toxin,
        gene_obj.swissprot,
        "|".join(gene_obj.interpro)
        ])



#-----------------------------------------------------
# Step 3
# Build DEG gene table for each condition
#-----------------------------------------------------
#
# gene_dict = defaultdict(list)
# for file, condition in zip(input_list, conditions_list):
#     with open(file) as f:
#         lines = f.readlines()
#         for x in lines:
#             gene_name = x.rstrip()
#             if gene_dict[gene_name]:
#                 gene_dict[gene_name].add_positive(condition)
#             else:
#                 gene_obj = DEG_obj(gene_name)
#                 gene_obj.set_conditions(conditions_list)
#                 gene_obj.set_positive(conditions_list)
#                 gene_obj.add_positive(condition)
#                 gene_dict[gene_name] = gene_obj
