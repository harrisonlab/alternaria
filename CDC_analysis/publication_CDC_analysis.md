
# 1. Alignment of Alternaria raw reads vs the each genome

Alignment of reads from a single run:

```bash
for Reference in $(ls repeat_masked/*/*/*/*_contigs_unmasked.fa | grep -e '1166' -e '650' -e '97.0013'); do
for StrainPath in $(ls -d qc_dna/paired/*/*); do
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/spades
Strain=$(echo $StrainPath | rev | cut -f1 -d '/' | rev)
Organism=$(echo $StrainPath | rev | cut -f2 -d '/' | rev)
F_Read=$(ls $StrainPath/F/*_trim.fq.gz)
R_Read=$(ls $StrainPath/R/*_trim.fq.gz)
for FileF in $(ls qc_rna/paired/*/*/F/*.fastq.gz); do
Jobs=$(qstat | grep 'sub_bwa' | grep 'qw'| wc -l)
while [ $Jobs -gt 1 ]; do
sleep 1m
printf "."
Jobs=$(qstat | grep 'sub_bwa' | grep 'qw'| wc -l)
done
echo $F_Read
echo $R_Read
Prefix="${Organism}_${Strain}"
OutDir=/home/scratch/groups/harrisonlab/alternaria/analysis/genome_alignment/bwa/$Organism/$Strain/vs_${Reference}
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/genome_alignment
qsub $ProgDir/bwa/sub_bwa.sh $Prefix $Reference $F_Read $R_Read $OutDir
done
done
done
```

Identify read coverage over each bp

```bash
  for Bam in $(ls analysis/genome_alignment/bowtie/*/*/vs_*/*_aligned_sorted.bam); do
    Strain=$(echo $Bam | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Bam | rev | cut -f4 -d '/' | rev)
    echo "$Organism - $Strain"
    OutDir=$(dirname $Bam)
    # samtools depth -aa $Bam > $OutDir/${Organism}_${Strain}_vs_12008_depth.tsv
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/genome_alignment/coverage_analysis
    $ProgDir/cov_by_window.py --cov $OutDir/${Organism}_${Strain}_vs_12008_depth.tsv > $OutDir/${Organism}_${Strain}_vs_12008_depth_10kb.tsv
    sed -i "s/$/\t$Strain/g" $OutDir/${Organism}_${Strain}_vs_12008_depth_10kb.tsv
  done
  OutDir=analysis/genome_alignment/bowtie/grouped
  mkdir -p $OutDir
  cat analysis/genome_alignment/bowtie/*/*/*vs_12008unmasked/*_*_vs_12008_depth_10kb.tsv > analysis/genome_alignment/bowtie/grouped/vs_12008_grouped_depth.tsv
```


This allowed assessment of read depth against the reference genomes.

The analysis was repeated for each isolate vs its own assembly to give an control value for the isolate:

```bash
for Reference in $(ls repeat_masked/*/*/*/*_contigs_unmasked.fa | grep -e '1166' -e '650' -e '97.0013'); do
RefIsolate=$(echo $Reference | rev | cur -f3 -d '/'| rev)
StrainPath=$(ls -d qc_dna/paired/*/$RefIsolate)
Organism=$(echo $StrainPath | rev | cut -f2 -d '/' | rev)
F_Read=$(ls $StrainPath/F/*_trim.fq.gz)
R_Read=$(ls $StrainPath/R/*_trim.fq.gz)
for FileF in $(ls qc_rna/paired/*/*/F/*.fastq.gz); do
Jobs=$(qstat | grep 'sub_bwa' | grep 'qw'| wc -l)
while [ $Jobs -gt 1 ]; do
sleep 1m
printf "."
Jobs=$(qstat | grep 'sub_bwa' | grep 'qw'| wc -l)
done
echo $F_Read
echo $R_Read
Prefix="${Organism}_${Strain}"
OutDir=analysis/genome_alignment/bwa/$Organism/$Strain/vs_${Reference}
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/genome_alignment
qsub $ProgDir/bwa/sub_bwa.sh $Prefix $Reference $F_Read $R_Read $OutDir
done
done
done
```
