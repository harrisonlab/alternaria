#!/usr/bin/bash
# 
# bowtie 1082_bowtie_index 
# 
# fastq-mcf $ILLUMINA_ADAPTERS $F_FILE $R_FILE -o F/"$QC_OUTFILE"_F.fastq -o R/"$QC_OUTFILE"_R.fastq -C 1000000 -u -k 20 -t 0.01 -q 30 -p 5

gunzip 1177_appended_1.fastq.gz
gunzip 1177_appended_2.fastq.gz

# bowtie 1082_bowtie_index 1177_appended_1.fastq.gz 1177_appended_2.fastq.gz
# 
# bowtie 1082_bowtie_index 1177_appended_1.fastq 1177_appended_1_aligned

#——
# 58% reads aligned
#——

#—————————————————

fastq-mcf /home/groups/harrisonlab/project_files/alternaria/raw_dna/paired/A.alternata_ssp._tenuissima/1177/F/1177_appended_1.fastq.gz /home/groups/harrisonlab/project_files/alternaria/raw_dna/paired/A.alternata_ssp._tenuissima/1177/R/1177_appended_2.fastq.gz -o 1177_F_trim.fq -o 1177_R_trim.fq -C 1000000 -u -k 20 -t 0.01 -q 30

bowtie 1082_bowtie_index 1177_F_trim.fq 1177_F_trim_out.txt --al 1177_F_trim_aln.fq --un 1177_trim_F_unaln.fq -p 16

#——
# 58% reads aligned
#——
# 
bowtie 1082_bowtie_index 1177_R_trim.fq 1177_R_trim_out.txt --al 1177_R_trim_aln.fq --un 1177_trim_R_unaln.fq -p 16

#———
# 56.38% reads aligned
#———
# 
# bowtie 1082_bowtie_index -1 1177_F_trim.fq -2 1177_R_trim.fq 1177_R_trim_out.txt --al 1177_R_trim_aln.fq --un 1177_trim_R_unaln.fq -p 16

#——
# Store the reads from 1177 with both F&R reads that didnt align to non-pathotype
#-----

cat 1177_trim_F_unaln.fq 1177_trim_R_unaln.fq | grep '@M00712' |  cut -d ' ' -f 1 | sort | uniq -d > duplicate_headers.txt

# cat duplicate_headers.txt | xargs -I{} grep -w {} 1177_F_trim.fq -A3 > 1177_F_LS.fq

# cat duplicate_headers.txt | xargs -I{} grep -w {} 1177_R_trim.fq -A3 >> 1177_R_LS.fq
count_nucl.pl -i 1177_F_trim.fq -i 1177_R_trim.fq -g 35

# while read line; do 
# 	grep -A3 -w "$line" 1177_trim_F_unaln.fq >>1177_F_LS2.fq 
# 	grep -A3 -w "$line" 1177_trim_R_unaln.fq >>1177_R_LS2.fq
# done<duplicate_headers.txt

SCRIPT_DIR=/home/armita/git_repos/emr_repos/tools/pathogen/lineage_specific_regions
$SCRIPT_DIR/find_novel_reads.py duplicate_headers.txt 1177_trim_F_unaln.fq > 1177_F_LS.fq
$SCRIPT_DIR/find_novel_reads.py duplicate_headers.txt 1177_trim_R_unaln.fq > 1177_R_LS.fq

velveth tenuissima_1177 41,81,10, -fastq -shortPaired -separate 1177_F_LS.fq 1177_R_LS.fq

for DIRECTORY in $(ls -d tenuissima_1177_*); do
	echo $DIRECTORY;
	cd $DIRECTORY;
	velvetg . -exp_cov 30 -cov_cutoff 15 -ins_length 600 -min_contig_lgth 500;
	cd ../;
	process_contigs.pl -i $DIRECTORY/contigs.fa -o $DIRECTORY;
done

--------------

cat 1177_contigs_unmasked.fa 743_contigs_unmasked.fa > contigs_appended_pathotypes.fa
cat 1082_contigs_unmasked.fa 1164_contigs_unmasked.fa 648_contigs_unmasked.fa 24350_contigs_unmasked.fa > contigs_appended_non-pathotypes.fa
bowtie-build contigs_appended_non-pathotypes.fa contigs_appended_non-pathotypes.index
bowtie contigs_appended_non-pathotypes.index ../1177_F_trim.fq 1177_F_trim_out.txt --al 1177_F_trim_aln.fq --un 1177_trim_F_unaln.fq -p 16
# reads processed: 7794525
# reads with at least one reported alignment: 5886321 (75.52%)
# reads that failed to align: 1908204 (24.48%)
bowtie contigs_appended_non-pathotypes.index ../1177_R_trim.fq 1177_R_trim_out.txt --al 1177_R_trim_aln.fq --un 1177_trim_R_unaln.fq -p 16
# reads processed: 7794525
# reads with at least one reported alignment: 5694535 (73.06%)
# reads that failed to align: 2099990 (26.94%)
bowtie-build contigs_appended_pathotypes.fa contigs_appended_pathotypes.index
bowtie contigs_appended_pathotypes.index 1177_trim_F_unaln.fq 1177_F_vs_pathotype --al 1177_F_vs_pathotype_aln.fq --un 1177_F_vs_pathotype_unaln.fq -p 16
# reads processed: 1908204
# reads with at least one reported alignment: 1587996 (83.22%)
# reads that failed to align: 320208 (16.78%)
bowtie contigs_appended_pathotypes.index 1177_trim_R_unaln.fq 1177_R_vs_pathotype --al 1177_R_vs_pathotype_aln.fq --un 1177_R_vs_pathotype_unaln.fq -p 16
# reads processed: 2099990
# reads with at least one reported alignment: 1534399 (73.07%)
# reads that failed to align: 565591 (26.93%)
cat 1177_F_vs_pathotype_aln.fq 1177_R_vs_pathotype_aln.fq | grep '@M00712' |  cut -d ' ' -f 1 | sort | uniq -d > 1177_pathotype_read_headers.txt
# 582468 lines in 1177_pathotype_read_headers.txt
SCRIPT_DIR=/home/armita/git_repos/emr_repos/tools/pathogen/lineage_specific_regions
$SCRIPT_DIR/find_novel_reads.py 1177_pathotype_read_headers.txt 1177_F_vs_pathotype_aln.fq > 1177_F_LS.fq
$SCRIPT_DIR/find_novel_reads.py 1177_pathotype_read_headers.txt 1177_R_vs_pathotype_aln.fq > 1177_R_LS.fq
# 
# velveth tenuissima_1177 41,81,10, -fastq -shortPaired -separate 1177_F_LS.fq 1177_R_LS.fq
# 
# for DIRECTORY in $(ls -d tenuissima_1177_*); do
# 	echo $DIRECTORY;
# 	cd $DIRECTORY;
# 	velvetg . -exp_cov 30 -cov_cutoff 10 -ins_length 600 -min_contig_lgth 500;
# 	cd ../;
# 	process_contigs.pl -i $DIRECTORY/contigs.fa -o $DIRECTORY;
# done
# 


bowtie contigs_appended_pathotypes.index --fr -1 1177_trim_F_unaln.fq -2 1177_trim_R_unaln.fq 1177_vs_pathotype -S -p 16
cp 1177_vs_pathotype 1177_vs_pathotype.sam
samtools sort 1177_vs_pathotype.bam 1177_vs_pathotype_sorted
samtools index 1177_vs_pathotype_sorted.bam 
samtools view -bS 1177_vs_pathotype.sam > 1177_vs_pathotype.bam
samtools view -F4 1177_vs_pathotype_sorted.bam
#fold contigs_appended_pathotypes.fa > contigs_appended_pathotypes.fold.fa
#samtools faidx contigs_appended_pathotypes.fold.fa
#samtools tview 1177_vs_pathotype_sorted.bam contigs_appended_pathotypes.fa

-----
1177 vs 1082

bowtie2-build ../1082_contigs_unmasked.fa 1082_bowtie_index
bowtie2 -x 1082_bowtie_index -1 ../1177_F_trim.fq -2 ../1177_R_trim.fq -S 1177_vs_1082.sam -p 16
# 7794525 reads; of these:
#   7794525 (100.00%) were paired; of these:
#     2385562 (30.61%) aligned concordantly 0 times
#     5372794 (68.93%) aligned concordantly exactly 1 time
#     36169 (0.46%) aligned concordantly >1 times
#     ----
#     2385562 pairs aligned concordantly 0 times; of these:
#       722774 (30.30%) aligned discordantly 1 time
#     ----
#     1662788 pairs aligned 0 times concordantly or discordantly; of these:
#       3325576 mates make up the pairs; of these:
#         3191109 (95.96%) aligned 0 times
#         122623 (3.69%) aligned exactly 1 time
#         11844 (0.36%) aligned >1 times
# 79.53% overall alignment rate
bowtie2 -x 1082_bowtie_index -1 ../1177_F_trim.fq -2 ../1177_R_trim.fq -S 1177_vs_1082.sam -p 16 --sensitive-local
# 7794525 reads; of these:
#   7794525 (100.00%) were paired; of these:
#     2074711 (26.62%) aligned concordantly 0 times
#     5546383 (71.16%) aligned concordantly exactly 1 time
#     173431 (2.23%) aligned concordantly >1 times
#     ----
#     2074711 pairs aligned concordantly 0 times; of these:
#       760342 (36.65%) aligned discordantly 1 time
#     ----
#     1314369 pairs aligned 0 times concordantly or discordantly; of these:
#       2628738 mates make up the pairs; of these:
#         2549614 (96.99%) aligned 0 times
#         42135 (1.60%) aligned exactly 1 time
#         36989 (1.41%) aligned >1 times
# 83.64% overall alignment rate
samtools view -bS 1177_vs_1082.sam > 1177_vs_1082.bam
samtools sort 1177_vs_1082.bam 1177_vs_1082_sorted
samtools index 1177_vs_1082_sorted.bam
samtools faidx 1082_contigs_unmasked.fa
samtools tview 1177_vs_1082_sorted.bam 1082_contigs_unmasked.fa
samtools view -F4 1177_vs_1082_sorted.bam

#----------------------------
#	Identify 1177 LS regions
#----------------------------
#	The methodology of alignment
#	Has been established.
#	Now assemble non-pathotype reads
#	Against a pathotype and
#	Extract LS regions

#	1082 reads vs 1177
fastq-mcf /home/armita/git_repos/emr_repos/tools/seq_tools/illumina_full_adapters.fa /home/groups/harrisonlab/project_files/alternaria/raw_dna/paired/A.alternata_ssp._tenuissima/1082/F/1082_appended_1.fastq.gz /home/groups/harrisonlab/project_files/alternaria/raw_dna/paired/A.alternata_ssp._tenuissima/1082/R/1082_appended_2.fastq.gz -o 1082_F_trim.fq -o 1082_R_trim.fq -C 1000000 -u -k 20 -t 0.01 -q 30
bowtie2-build ../appended/1177_contigs_unmasked.fa 1177_bowtie_index
bowtie2 -x 1177_bowtie_index -1 1082_F_trim.fq -2 1082_R_trim.fq -S 1082_vs_1177.sam -p 16
# 2544523 reads; of these:
#   2544523 (100.00%) were paired; of these:
#     1522357 (59.83%) aligned concordantly 0 times
#     979641 (38.50%) aligned concordantly exactly 1 time
#     42525 (1.67%) aligned concordantly >1 times
#     ----
#     1522357 pairs aligned concordantly 0 times; of these:
#       1268224 (83.31%) aligned discordantly 1 time
#     ----
#     254133 pairs aligned 0 times concordantly or discordantly; of these:
#       508266 mates make up the pairs; of these:
#         429386 (84.48%) aligned 0 times
#         38470 (7.57%) aligned exactly 1 time
#         40410 (7.95%) aligned >1 times
samtools view -bS 1082_vs_1177.sam > 1082_vs_1177.bam
samtools sort 1082_vs_1177.bam 1082_vs_1177_sorted
samtools index 1082_vs_1177_sorted.bam
samtools faidx 1177_contigs_unmasked.fa
samtools tview 1082_vs_1177_sorted.bam 1177_contigs_unmasked.fa
samtools idxstats 1082_vs_1177_sorted.bam > 1082_vs_1177_sorted_indexstats.csv
#samtools view -F4 1082_vs_1177_sorted.bam
/home/armita/git_repos/emr_repos/scripts/alternaria/assembly/divide_col.py 1082_vs_1177_sorted_indexstats.csv 1 2 > 1082_vs_1177_sorted_indexstats_coverage.csv
printf "occurence\taligned_reads_per_base\n" > aligned_reads_per_base.csv
cat 1082_vs_1177_sorted_indexstats_coverage.csv | cut -f5 | sort -n | uniq -c | sed 's/ *//' | sed 's/ /\t/g' >> aligned_reads_per_base.csv | cut -f1 | python -c"import sys; print(sum(map(int, sys.stdin)))"
# These contigs represent 294200bp of the assembly. There were 521530 bp of unique contigs in 1082.

cat 1082_vs_1177_sorted_indexstats_coverage.csv | cut -f1,5 | grep -w '0' | grep -v '*' | cut -f1 > 1177_ls_contigs.txt

