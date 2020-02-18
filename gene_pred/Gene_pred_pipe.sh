#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -pe smp 16
#$ -l virtual_free=1G

# Pipeline used for Gene prediction used in Alternaria using rnaSeq data assembled with trinity.
# Training was performed using RNASeq data.


# USAGE="Gene_pred_pipe.sh <transcriptome_assembly.fa> <genomic_contigs.fa>

#------------------------------------------------------
# 		Step 0.		Initialise values
#------------------------------------------------------

AUG_DIR=$(which augustus)
PASA_DIR=$(which pasa | sed s%/pasa%%)

TRANSCRIPTOME=$1
#TRANSCRIPTOME=assembly/trinity/A.alternata_ssp._gaisen/650/650_rna_contigs/Trinity.fasta
#TRANSCRIPTOME=assembly/trinity/genbank/P.cactorum/P.cactorum_rna_contigs/Trinity.fasta 

GENOMIC_CONTIGS=$2
#GENOMIC_CONTIGS=assembly/velvet/A.alternata_ssp._gaisen/650/650_assembly.41/sorted_contigs.fa
#GENOMIC_CONTIGS=repeat_masked/P.cactorum/10300/version1_repmask/10300_contigs_unmasked.fa
ORGANISM=$(echo $GENOMIC_CONTIGS | rev | cut -d "/" -f4 | rev)
STRAIN=$(echo $GENOMIC_CONTIGS | rev | cut -d "/" -f3 | rev)

CUR_PATH=$PWD
#WORK_DIR=$TMPDIR/$STRAIN
WORK_DIR=/tmp/$STRAIN

PASA_DB="$STRAIN"_new_db
# When Pasa makes a MYSQL database it wont work if there is already a database of this name.
# For this reason the PASA_DB value needs to be changed everytime this script is run.
# IMPORTANT - find a way to remove MYSQL databases.
# 			- MYSQL is not being called correctly from the slave nodes. 
#				Forcing the script to be run from the head.

#------------------------------------------------------
# 		Step 1.		Extend reads using FLASH
#------------------------------------------------------
#flash.sh <F_FILE.fastq.gz> <R_FILE.fastq.gz>
#------------------------------------------------------
# 		Step 2.		Assemble transcriptome using Trinity
#------------------------------------------------------
#transcriptome_assembly_trinity.sh <F_FILE.fastq.gz> <R_FILE.fastq.gz>
#------------------------------------------------------
# 		Step 3.		Generate spliced alignents to genome using PASA
#------------------------------------------------------
#run seqclean to identify and discard polyA trim vectors, and low quality sequences 

mkdir -p gene_pred/pasa/$ORGANISM/$STRAIN/
cd gene_pred/pasa/$ORGANISM/$STRAIN/
seqclean $CUR_PATH/$TRANSCRIPTOME

#outputs files <name>.cidx <name>.clean, <name>.cln, seql_<name>.log, err_seql_<name>.log, cleaning_1, outparts_cln.sort

cp $PASA_DIR/../pasa_conf/pasa.alignAssembly.Template.txt .
perl -pi -e "s%<__MYSQLDB__>%$PASA_DB%g" pasa.alignAssembly.Template.txt				
perl -pi -e "s%<__MIN_PERCENT_ALIGNED__>%90%g" pasa.alignAssembly.Template.txt			
perl -pi -e "s%<__MIN_AVG_PER_ID__>%95%g" pasa.alignAssembly.Template.txt				


Launch_PASA_pipeline.pl -c pasa.alignAssembly.Template.txt -C -R -g $CUR_PATH/$GENOMIC_CONTIGS -t $CUR_PATH/$TRANSCRIPTOME --ALIGNERS blat,gmap


pasa_asmbls_to_training_set.dbi --pasa_transcripts_fasta $PASA_DB.assemblies.fasta  --pasa_transcripts_gff $PASA_DB.pasa_assemblies.gff3

grep "complete" $PASA_DB.assemblies.fasta.transdecoder.pep | cut -d " " -f1 | cut -c 2- > complete_genes.txt
#cd $CUR_PATH

#------------------------------------------------------
# 		Step 4.		Convert .gff output to .gb format
#------------------------------------------------------
# gff3_2_auggff.pl <infile.gff3> > <augustus_format_outfile.gff3>

/home/armita/git_repos/emr_repos/tools/gene_prediction/augustus/gff3_2_auggff.pl $PASA_DB.assemblies.fasta.transdecoder.genome.gff3 > $PASA_DB.assemblies.transdecoder.genome.aug.gff3

# As part of this the flanking region length between genes needs to be set. This is recomended
# to be set to half the average gene size (as used on the Augustus web server).
# In autoAUG.pl it is set to 4000bp.
MAX_FLANK_DNA=4000

# #gff2gbSmallDNA.pl gene_pred/pasa/$ORGANISM/$STRAIN/temp_db1.pasa_assemblies.gff3 assembly/velvet/A.alternata_ssp._gaisen/650/650_assembly.41/sorted_contigs.fa $MAX_FLANK_DNA "$ORGANISM"_"$STRAIN"_evidence.gb
# gff2gbSmallDNA.pl gene_pred/pasa/temp_db1.pasa_assemblies.gff3 assembly/velvet/A.alternata_ssp._gaisen/650/650_assembly.41/sorted_contigs.fa 4000 A.alternata_ssp.gaisen_650_evidence.gb
# gene_size.pl A.alternata_ssp.gaisen_650_evidence.gb > tmp.txt
# MAX_FLANK_DNA = $(tail -n 1 tmp.txt)
# rm tmp.txt

gff2gbSmallDNA.pl $PASA_DB.assemblies.fasta.transdecoder.genome.gff3 ../../../../$GENOMIC_CONTIGS 4000 "$ORGANISM"_"$STRAIN"_evidence.gb --good=complete_genes.txt


#------------------------------------------------------
# 		Step 5.		Extract a subset of the aligned reads to use as a test set
#------------------------------------------------------
# The genes.gb must be in the directory this command is being run from
randomSplit.pl "$ORGANISM"_"$STRAIN"_evidence.gb 100



#------------------------------------------------------
# 		Step 2.		Create a metafile for the new species
#------------------------------------------------------

new_species.pl --species="$ORGANISM"_"$STRAIN"

perl -pi -e "s%codingseq           off%codingseq           on%g" ~/prog/augustus/config/species/"$ORGANISM"_"$STRAIN"/"$ORGANISM"_"$STRAIN"_parameters.cfg
perl -pi -e "s%stopCodonExcludedFromCDS false%stopCodonExcludedFromCDS true%g" ~/prog/augustus/config/species/"$ORGANISM"_"$STRAIN"/"$ORGANISM"_"$STRAIN"_parameters.cfg
perl -pi -e "s%alternatives-from-evidence  false%alternatives-from-evidence  true%g" ~/prog/augustus/config/species/"$ORGANISM"_"$STRAIN"/"$ORGANISM"_"$STRAIN"_parameters.cfg



#------------------------------------------------------
# 		Step 2.		Train Augustus using aligned reads
#------------------------------------------------------
#cp $AUG_DIR/config/species/generic/generic_metapars.cfg .

etraining --species="$ORGANISM"_"$STRAIN" "$ORGANISM"_"$STRAIN"_evidence.gb.train

augustus --species="$ORGANISM"_"$STRAIN" "$ORGANISM"_"$STRAIN"_evidence.gb.train | tee "$ORGANISM"_"$STRAIN"_sum.txt 

# at the bottom of this output file is a table with assembly accuracy statistics. The 
# key values to note are: gene level sensitivity (proportion of genes predicted exactly)
# exon level sensitivity (proportion of exons predicted exactly); exon level specificity 
# (proportion of exons predicted exactly as in the test set).
#
# ---------------------------------------------|
# nucleotide level |       0.943 |       0.474 |
# ---------------------------------------------/
# 
# ----------------------------------------------------------------------------------------------------------\
#            |  #pred |  #anno |      |    FP = false pos. |    FN = false neg. |             |             |
#            | total/ | total/ |   TP |--------------------|--------------------| sensitivity | specificity |
#            | unique | unique |      | part | ovlp | wrng | part | ovlp | wrng |             |             |
# ----------------------------------------------------------------------------------------------------------|
#            |        |        |      |                340 |                 50 |             |             |
# exon level |    537 |    247 |  197 | ------------------ | ------------------ |       0.798 |       0.367 |
#            |    537 |    247 |      |   31 |    2 |  307 |   32 |    2 |   16 |             |             |
# ----------------------------------------------------------------------------------------------------------/
# 
# ----------------------------------------------------------------------------\
# transcript | #pred | #anno |   TP |   FP |   FN | sensitivity | specificity |
# ----------------------------------------------------------------------------|
# gene level |   215 |   100 |   62 |  153 |   38 |        0.62 |       0.288 |
# ----------------------------------------------------------------------------/

optimize_augustus.pl --species="$ORGANISM"_"$STRAIN" "$ORGANISM"_"$STRAIN"_evidence.gb.train

#------------------------------------------------------
# 		Step 2.		Train Maker using aligned reads
#------------------------------------------------------

#------------------------------------------------------
# 		Step 2.		Assemble transcriptome using Trinity
#------------------------------------------------------