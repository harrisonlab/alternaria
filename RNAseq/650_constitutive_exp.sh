# Alternaria alignment and expression pipe.

# Are there genes in 1166 that are predicted but do not show constitutive
# expression in RNAseq data.

#		Later test
# # 	Are there genes in 1166 that show no coverage by 1164 sequence reads.
# # 	Do these genes show constitutive expression in 1166 culture?


cd /home/groups/harrisonlab/project_files/alternaria/assembly
mkdir tmp
cd tmp
cp ../../repeat_masked/A.alternata_ssp._gaisen/650/650_assembly.41_repmask/650_contigs_hardmasked.fa .
bowtie2-build 650_contigs_hardmasked.fa 650_contigs_hardmasked_index
bowtie2 -x 650_contigs_hardmasked_index -1 ../../qc_rna/paired/A.alternata_ssp._gaisen/650/F/650_qc_F.fastq -2 ../../qc_rna/paired/A.alternata_ssp._gaisen/650/R/650_qc_R.fastq -S 650_RNA_alignment.sam -p 4
samtools view -bS 650_RNA_alignment.sam > 650_RNA_alignment.bam
samtools sort 650_RNA_alignment.bam 650_RNA_alignment_sorted
samtools index 650_RNA_alignment_sorted.bam
samtools faidx 650_contigs_hardmasked.fa
samtools tview 650_RNA_alignment_sorted.bam 650_contigs_hardmasked.fa
samtools idxstats 650_RNA_alignment_sorted.bam > 650_RNA_alignment_sorted_indexstats.csv
cat ../../gene_pred/augustus/A.alternata_ssp._gaisen/650/650_augustus_preds.gtf | grep -v '#' | grep -v -E "^Error" > 650_augustus.gff
bedtools intersect -a 650_RNA_alignment_sorted.bam -b 650_augustus.gff > 650_intersect1.txt
bedtools intersect -bed -u -a 650_RNA_alignment_sorted.bam -b 650_augustus.gff > 650_intersect2.txt
bedtools intersect -c -a 650_augustus.gff -b 650_RNA_alignment_sorted.bam > 650_intersect3.bed
cat 650_intersect3.bed | grep -E "\".0$" > 650_unexpressed_features.bed


