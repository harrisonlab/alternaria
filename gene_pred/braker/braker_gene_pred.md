
This document details commands used for an initial run of Braker to train and
predict gene models.

Braker is a pipeline based upon Augustus and GeneMark-est.

## 1) RNA QC

Data quality was visualised using fastqc:
```bash
	for RawData in raw_rna/paired/*/*/*/*.fastq.gz; do
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc
		echo $RawData;
		qsub $ProgDir/run_fastqc.sh $RawData
	done
```

Trimming was performed on data to trim adapters from
sequences and remove poor quality data. This was done with fastq-mcf

```bash
	for StrainPath in raw_rna/paired/*/*; do
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/rna_qc
		IlluminaAdapters=/home/armita/git_repos/emr_repos/tools/seq_tools/ncbi_adapters.fa
		ReadsF=$(ls $StrainPath/F/*.fastq.gz)
		ReadsR=$(ls $StrainPath/R/*.fastq.gz)
		echo $ReadsF
		echo $ReadsR
		qsub $ProgDir/rna_qc_fastq-mcf.sh $ReadsF $ReadsR $IlluminaAdapters RNA
	done
```

Data quality was visualised once again following trimming:
```bash
	for TrimData in qc_rna/paired/*/*/*/*.fastq.gz; do
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc
		echo $RawData;
		qsub $ProgDir/run_fastqc.sh $TrimData
	done
```


# 2) RNA Alignment

Insert sizes of the RNA seq library were unknown until a draft alignment could
be made. To do this tophat and cufflinks were run, aligning the reads against a
single genome. The fragment length and stdev were printed to stdout while
cufflinks was running.

```bash
for Strain in $(650); do
	Assembly=$(ls repeat_masked/*/$Strain/filtered_contigs_repmask/*_contigs_unmasked.fa)
	Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
	echo "$Organism"
	echo "$Strain"
	ProgDir=/home/armita/git_repos/emr_repos/scripts/alternaria/gene_pred/braker
	FileF1=qc_rna/paired/A.alternata_ssp._gaisen/650/F/Pear_S1_L001_R1_001_trim.fq.gz
	FileR1=qc_rna/paired/A.alternata_ssp._gaisen/650/R/Pear_S1_L001_R2_001_trim.fq.gz
	FileF2=qc_rna/paired/A.alternata_ssp._tenuissima/1166/F/Apple_S2_L001_R1_001_trim.fq.gz
	FileR2=qc_rna/paired/A.alternata_ssp._tenuissima/1166/R/Apple_S2_L001_R2_001_trim.fq.gz
	OutDir=aligment/$Organism/$Strain
	qsub $ProgDir/alt_tophat_alignment.sh $Assembly $FileF1 $FileR1 $FileF2 $FileR2 $OutDir
done
```

Cufflinks was run to produce the fragment length and stdev statistics:

```bash
	for Assembly in $(ls repeat_masked/*/650/*/*_contigs_unmasked.fa); do
		Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
		Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
		echo "$Organism - $Strain"
		for RNADir in $(ls -d qc_rna/paired/*/*); do
			Timepoint=$(echo $RNADir | rev | cut -f1 -d '/' | rev)
			echo "$Timepoint"
			FileF=$(ls $RNADir/F/*_trim.fq.gz)
			FileR=$(ls $RNADir/R/*_trim.fq.gz)
			OutDir=alignment/$Organism/$Strain/$Timepoint
			ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/RNAseq
			qsub $ProgDir/tophat_alignment.sh $Assembly $FileF $FileR $OutDir
		done
	done
```

Output from stdout included:
```
	Processed 22484 loci.                        [*************************] 100%
	Map Properties:
	     Normalized Map Mass: 50507412.55
	     Raw Map Mass: 50507412.55
	     Fragment Length Distribution: Empirical (learned)
	                   Estimated Mean: 181.98
	                Estimated Std Dev: 78.39
	[13:02:48] Assembling transcripts and estimating abundances.
	Processed 22506 loci.                        [*************************] 100%
```

The Estimated Mean: X allowed calculation of of the mean insert gap to be
Xbp X-(152*2) where 152 was the mean read length. This was provided to tophat
on a second run (as the -r option) along with the fragment length stdev to
increase the accuracy of mapping.


A second round of tophat RNAseq alignment was performed, this time using
parameters determined above.

```bash
	for Assembly in $(ls repeat_masked/*/*/*/*_contigs_unmasked.fa); do
		Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
		Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
		echo "$Organism - $Strain"
		for RNADir in $(ls -d qc_rna/paired/*/*); do
			Timepoint=$(echo $RNADir | rev | cut -f1 -d '/' | rev)
			echo "$Timepoint"
			FileF=$(ls $RNADir/F/*_trim.fq.gz)
			FileR=$(ls $RNADir/R/*_trim.fq.gz)
			OutDir=alignment/$Organism/$Strain/$Timepoint
			InsertGap=''
			InsertStdDev=''
			ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/RNAseq
			qsub $ProgDir/tophat_alignment.sh $Assembly $FileF $FileR $OutDir $InsertGap $InsertStdDev
		done
	done
```

These alignments were used to perform Braker gene prediction. As alternaria is
a fungus, the fungi specific gene prediction option was turned on.

```bash
	for Assembly in $(ls repeat_masked/*/Fus2/*/*_contigs_softmasked.fa); do
		Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
		Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
		echo "$Organism - $Strain"
		mkdir -p merge alignment/$Organism/$Strain/concatenated
		samtools merge -f alignment/$Organism/$Strain/concatenated/concatenated.bam \
			alignment/$Organism/$Strain/55_72hrs_rep1/accepted_hits.bam \
			alignment/$Organism/$Strain/55_72hrs_rep2/accepted_hits.bam \
			alignment/$Organism/$Strain/55_72hrs_rep3/accepted_hits.bam \
			alignment/$Organism/$Strain/FO47_72hrs_rep1/accepted_hits.bam \
			alignment/$Organism/$Strain/FO47_72hrs_rep2/accepted_hits.bam \
			alignment/$Organism/$Strain/FO47_72hrs_rep3/accepted_hits.bam \
			alignment/$Organism/$Strain/Fus2_0hrs_prelim/accepted_hits.bam \
			alignment/$Organism/$Strain/Fus2_16hrs_prelim/accepted_hits.bam \
			alignment/$Organism/$Strain/Fus2_24hrs_prelim_rep1/accepted_hits.bam \
			alignment/$Organism/$Strain/Fus2_36hrs_prelim/accepted_hits.bam \
			alignment/$Organism/$Strain/Fus2_48hrs_prelim/accepted_hits.bam \
			alignment/$Organism/$Strain/Fus2_4hrs_prelim/accepted_hits.bam \
			alignment/$Organism/$Strain/Fus2_72hrs_prelim/accepted_hits.bam \
			alignment/$Organism/$Strain/Fus2_72hrs_rep1/accepted_hits.bam \
			alignment/$Organism/$Strain/Fus2_72hrs_rep2/accepted_hits.bam \
			alignment/$Organism/$Strain/Fus2_72hrs_rep3/accepted_hits.bam \
			alignment/$Organism/$Strain/Fus2_8hrs_prelim/accepted_hits.bam \
			alignment/$Organism/$Strain/Fus2_96hrs_prelim/accepted_hits.bam \
			alignment/$Organism/$Strain/Fus2_CzapekDox/accepted_hits.bam \
			alignment/$Organism/$Strain/Fus2_GlucosePeptone/accepted_hits.bam \
			alignment/$Organism/$Strain/Fus2_PDA/accepted_hits.bam \
			alignment/$Organism/$Strain/Fus2_PDB/accepted_hits.bam
		OutDir=gene_pred/braker/$Organism/"$Strain"
		AcceptedHits=alignment/$Organism/$Strain/concatenated/concatenated.bam
		GeneModelName="$Organism"_"$Strain"_braker
		rm -r /home/armita/prog/augustus-3.1/config/species/"$Organism"_"$Strain"_braker
		ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/braker1
		qsub $ProgDir/sub_braker_fungi.sh $Assembly $OutDir $AcceptedHits $GeneModelName
	done
```
