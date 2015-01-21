#!/usr/bin/bash
bowtie 1082_bowtie_index 
1177_appended_1.fastq.gz 1177_appended_2.fastq.gz 1177_1177 117g 11gz 1qgz q.gztq.gstq.astqfast.fasowtie_index 1177_appended_1.fastq.gz 1177_appended_2.fa

fastq-mcf $ILLUMINA_ADAPTERS $F_FILE $R_FILE -o F/"$QC_OUTFILE"_F.fastq -o R/"$QC_OUTFILE"_R.fastq -C 1000000 -u -k 20 -t 0.01 -q 30 -p 5

gunzip 1177_appended_1.fastq.gz
gunzip 1177_appended_2.fastq.gz

bowtie 1082_bowtie_index 1177_appended_1.fastq.gz 1177_appended_2.fastq.gz

bowtie 1082_bowtie_index 1177_appended_1.fastq 1177_appended_1_aligned

#——
# 58% reads aligned
#——

#—————————————————

fastq-mcf /home/groups/harrisonlab/project_files/alternaria/raw_dna/paired/A.alternata_ssp._tenuissima/1177/F/1177_appended_1.fastq.gz /home/groups/harrisonlab/project_files/alternaria/raw_dna/paired/A.alternata_ssp._tenuissima/1177/R/1177_appended_2.fastq.gz -o 1177_F_trim.fq -o 1177_R_trim.fq -C 1000000 -u -k 20 -t 0.01 -q 30

bowtie 1082_bowtie_index 1177_F_trim.fq 1177_F_trim_out.txt --al 1177_F_trim_aln.fq --un 1177_trim_F_unaln.fq 

#——
# 58% reads aligned
#——

bowtie 1082_bowtie_index 1177_R_trim.fq 1177_R_trim_out.txt --al 1177_R_trim_aln.fq --un 1177_trim_R_unaln.fq -p 16

#———
# 56.38% reads aligned
#———

bowtie 1082_bowtie_index -1 1177_F_trim.fq -2 1177_R_trim.fq 1177_R_trim_out.txt --al 1177_R_trim_aln.fq --un 1177_trim_R_unaln.fq -p 16

#——

#Store the reads from 1177 with both F&R reads that didnt align to non-pathotype

cat 1177_trim_F_unaln.fq 1177_trim_R_unaln.fq | grep '@M00712' |  cut -d ' ' -f 1 | sort | uniq -d > duplicate_headers.txt

cat duplicate_headers.txt | xargs -I{} grep -w {} 1177_F_trim.fq -A3 > 1177_F_LS.fq

cat duplicate_headers.txt | xargs -I{} grep -w {} 1177_R_trim.fq -A3 >> 1177_R_LS.fq
count_nucl.pl -i 1177_F_trim.fq -i 1177_R_trim.fq -g 35

/home/groups/harrisonlab/project_files/alternaria/bowtie_alternaria$ while read line; do
	grep -A4 -w "$line" 1177_trim_F_unaln.fq >>1177_F_LS2.fq
done<duplicate_headers.txt

/home/groups/harrisonlab/project_files/alternaria/bowtie_alternaria$ while read line; do
	grep -A3 -w "$line" 1177_trim_R_unaln.fq >>1177_R_LS2.fq
done<duplicate_headers.txt

