#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -pe smp 4
#$ -l virtual_free=2G

Usage="alt_tophat_align.sh Genome.fa ReadF1.fq ReadR1.fq ReadF2.fq ReadR2.fq <Output_directory>"
Genome=$1
ReadF1=$2
ReadR1=$3
ReadF2=$4
ReadR2=$5
OutDir=$6

Organism=$(echo $ReadF | rev | cut -d "/" -f4 | rev)
Strain=$(echo $ReadF | rev | cut -d "/" -f3 | rev)

OutName=$(basename "$Strain"_tophat)

CurPath=$PWD
WorkDir=$TMPDIR/tophat_align

mkdir -p $WorkDir
cd $WorkDir

cp $CurPath/$Genome "$OutName".fa

cp $CurPath/$ReadF1 ReadF1
gunzip -cf ReadF1 > ReadF1.fastq
cp $CurPath/$ReadR1 ReadR1
gunzip -cf ReadR1 > ReadR1.fastq
cp $CurPath/$ReadF2 ReadF2
gunzip -cf ReadF2 > ReadF2.fastq
cp $CurPath/$ReadR2 ReadR2
gunzip -cf ReadR2 > ReadR2.fastq

bowtie2-build "$OutName".fa "$OutName"
tophat -o "$OutName" -p 8 --library-type fr-unstranded "$OutName" ReadF1.fastq,ReadR1.fastq ReadF2.fastq,ReadR2.fastq 2>&1 | tee "$Strain"_report.txt
rm "$OutName".fa*
rm ReadF1
rm ReadR1
rm ReadF2
rm ReadR2
rm ReadF1.fastq
rm ReadR1.fastq
rm ReadF2.fastq
rm ReadR2.fastq

mkdir -p $CurPath/$OutDir
cp -r $WorkDir/$OutName/* $CurPath/$OutDir/.
rm -r $TMPDIR
