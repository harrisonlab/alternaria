#!/bin/bash

QUERY=/home/groups/harrisonlab/project_files/idris/analysis/blast_homology/PHI_36_accessions

PROGRAM=/home/armita/git_repos/emr_repos/tools/pathogen/blast/blast_pipe.sh 

qsub $PROGRAM $QUERY assembly/velvet/A.alternata_ssp._tenuissima/1082/1082_assembly.41/sorted_contigs.fa

qsub $PROGRAM $QUERY assembly/velvet/A.alternata_ssp._tenuissima/1164/1164_assembly.41/sorted_contigs.fa

qsub $PROGRAM $QUERY assembly/velvet/A.alternata_ssp._tenuissima/1177/1177_assembly.41/sorted_contigs.fa

qsub $PROGRAM $QUERY assembly/velvet/A.alternata_ssp._tenuissima/24350/24350_assembly.41/sorted_contigs.fa

qsub $PROGRAM $QUERY assembly/velvet/A.alternata_ssp._tenuissima/635/635_assembly.41/sorted_contigs.fa

qsub $PROGRAM $QUERY assembly/velvet/A.alternata_ssp._tenuissima/648/648_assembly.41/sorted_contigs.fa

qsub $PROGRAM $QUERY assembly/velvet/A.alternata_ssp._tenuissima/743/743_assembly.41/sorted_contigs.fa

qsub $PROGRAM $QUERY assembly/velvet/A.alternata_ssp._arborescens/675/675_assembly.41/sorted_contigs.fa

qsub $PROGRAM $QUERY assembly/velvet/A.alternata_ssp._arborescens/97.0013/97.0013_assembly.41/sorted_contigs.fa

qsub $PROGRAM $QUERY assembly/velvet/A.alternata_ssp._arborescens/97.0016/97.0016_assembly.41/sorted_contigs.fa

qsub $PROGRAM $QUERY assembly/velvet/A.alternata_ssp._gaisen/650/650_assembly.41/sorted_contigs.fa

