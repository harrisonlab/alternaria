
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
```bash
  for Strain in 1166 650 97.0013; do
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
