#$ -S /bin/bash
#$ -cwd
#$ -pe smp 24
#$ -l virtual_free=1.25G
#$ -l h=blacklace01.blacklace|blacklace02.blacklace|blacklace04.blacklace|blacklace05.blacklace|blacklace06.blacklace|blacklace07.blacklace|blacklace08.blacklace|blacklace09.blacklace|blacklace10.blacklace|blacklace12.blacklace

# Testing parallelisation of GATk HaplotypeCaller - may crash. (It did not! Resulted in 2x speedup)
# NOTE: this is a haploid organism. For diploid organism, change "ploidy" argument to 2.
# Changes required in the script:
# VARIABLES
# Reference - the genome reference used in read mapping.
# INSIDE THE GATK command:
# To specify which BAM mapping files (Out1 from pre_SNP_calling_cleanup.sh, RefName ending with "_rg" -> that is, with
# read group added) are to be used in SNP calling, use the -I argument with full path to each file following after that.
# Each new BAM file has to be specified after a separate -I

# Project=/home/groups/harrisonlab/project_files/idris
Project=/data/scratch/armita/alternaria
OutDir=analysis/popgen/SNP_calling
Reference=$(ls $Project/repeat_masked/A.alternata_ssp_tenuissima/1166/filtered_contigs/1166_contigs_unmasked.fa)
RefName=$(basename "$Reference")
Out1="${RefName%.*}_temp.vcf"
Out2="${RefName%.*}.vcf"

ProgDir=/home/sobczm/bin/GenomeAnalysisTK-3.6

java -jar $ProgDir/GenomeAnalysisTK.jar \
     -R $Reference \
     -T HaplotypeCaller \
     -ploidy 1 \
     -nct 24 \
     --allow_potentially_misencoded_quality_scores \
     -I $Project/analysis/popgen/A.alternata_ssp._arborescens/675/675_vs_1166_aligned_sorted_nomulti_proper_sorted_nodup_rg.bam \
     -I $Project/analysis/popgen/A.alternata_ssp._arborescens/97.0013/97.0013_vs_1166_aligned_sorted_nomulti_proper_sorted_nodup_rg.bam \
     -I $Project/analysis/popgen/A.alternata_ssp._arborescens/97.0016/97.0016_vs_1166_aligned_sorted_nomulti_proper_sorted_nodup_rg.bam \
     -I $Project/analysis/popgen/A.alternata_ssp._gaisen/650/650_vs_1166_aligned_sorted_nomulti_proper_sorted_nodup_rg.bam \
     -I $Project/analysis/popgen/A.alternata_ssp._tenuissima/1082/1082_vs_1166_aligned_sorted_nomulti_proper_sorted_nodup_rg.bam \
     -I $Project/analysis/popgen/A.alternata_ssp._tenuissima/1164/1164_vs_1166_aligned_sorted_nomulti_proper_sorted_nodup_rg.bam \
     -I $Project/analysis/popgen/A.alternata_ssp._tenuissima/1166/1166_vs_1166_aligned_sorted_nomulti_proper_sorted_nodup_rg.bam \
     -I $Project/analysis/popgen/A.alternata_ssp._tenuissima/1177/1177_vs_1166_aligned_sorted_nomulti_proper_sorted_nodup_rg.bam \
     -I $Project/analysis/popgen/A.alternata_ssp._tenuissima/24350/24350_vs_1166_aligned_sorted_nomulti_proper_sorted_nodup_rg.bam \
     -I $Project/analysis/popgen/A.alternata_ssp._tenuissima/635/635_vs_1166_aligned_sorted_nomulti_proper_sorted_nodup_rg.bam \
     -I $Project/analysis/popgen/A.alternata_ssp._tenuissima/648/648_vs_1166_aligned_sorted_nomulti_proper_sorted_nodup_rg.bam \
     -I $Project/analysis/popgen/A.alternata_ssp._tenuissima/743/743_vs_1166_aligned_sorted_nomulti_proper_sorted_nodup_rg.bam \
     -o $Out1

#Break down complex SNPs into primitive ones with VariantsToAllelicPrimitives
#This tool will take an MNP (e.g. ACCCA -> TCCCG) and break it up into separate records for each component part (A-T and A->G).
#This tool modifies only bi-allelic variants.

java -jar $ProgDir/GenomeAnalysisTK.jar \
   -T VariantsToAllelicPrimitives \
   -R $Reference \
   -V $Out1 \
   -o $Out2 \


#####################################
# Notes on GATK parallelisation
#####################################
# http://gatkforums.broadinstitute.org/gatk/discussion/1975/how-can-i-use-parallelism-to-make-gatk-tools-run-faster
