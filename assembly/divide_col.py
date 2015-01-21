#!/usr/bin/env python
import os
import sys

infile = sys.argv[1]
# divide this column:
div_column = int(sys.argv[2])
# by this column:
by_column = int(sys.argv[3])

usage = "div_col.py <infile.csv> [div_this_col] [by_this_col]"

with open (infile, 'r') as infile_fh:
		for read_line in infile_fh:
			cur_line = read_line.strip()
			contig_lgth = cur_line.split("\t")[div_column]
			no_reads = cur_line.split("\t")[by_column]
			if no_reads == '0':
				coverage = 0
			else: 
				coverage = int(contig_lgth) / int(no_reads)
			print cur_line + "\t" + str(coverage)
