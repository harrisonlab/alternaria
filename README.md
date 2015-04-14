#Alternaria
==========

Scripts used for the analysis of Alternaria genomes
Note - all this work was performed in the directory:
/home/groups/harrisonlab/project_files/alternaria

The following is a summary of the work presented in this Readme.

The following processes were applied to Alternaria genomes prior to analysis:
Data qc
Genome assembly
Repeatmasking
Gene prediction
Functional annotation

Analyses performed on these genomes involved BLAST searching for:
Known toxin genes

CDC contigs were identified using:
Alignment of raw reads to assembled genomes
Assembly of remaining reads


#Data qc

programs:
  fastqc
  fastq-mcf
  kmc

Data quality was visualised using fastqc:
```shell
	for RawData in raw_dna/paired/*/*/*/*.fastq*?; do 
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc
		echo $RawData; 
		qsub $ProgDir/run_fastqc.sh $RawData
	done
```

Trimming was performed on data to trim adapters from 
sequences and remove poor quality data. This was done with fastq-mcf

```shell
	for StrainPath in raw_dna/paired/*/*; do 
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/rna_qc
		IlluminaAdapters=/home/armita/git_repos/emr_repos/tools/seq_tools/illumina_full_adapters.fa
		ReadsF=$(ls $StrainPath/F/*.fastq*)
		ReadsR=$(ls $StrainPath/R/*.fastq*)
		echo $ReadsF
		echo $ReadsR
		qsub $ProgDir/rna_qc_fastq-mcf.sh $ReadsF $ReadsR $IlluminaAdapters DNA
	done
```

kmer counting was performed using kmc
This allowed estimation of sequencing depth and total genome size

```shell
	for TrimPath in qc_dna/paired/*/*; do 
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc
		TrimF=$(ls $TrimPath/F/*.fastq*)
		TrimR=$(ls $TrimPath/R/*.fastq*)
		echo $TrimF
		echo $TrimR
		qsub $ProgDir/kmc_kmer_counting.sh $TrimF $TrimR
	done
```

#Assembly

Assembly was performed using Velvet

A range of hash lengths were used and the best assembly selected for subsequent analysis

```shell
	for TrimPath in qc_dna/paired/*/*; do
	ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/velvet
	Strain=$(printf $TrimPath | rev | cut -f1 -d '/' | rev)
	
	MinHash=41
	MaxHash=81
	HashStep=2
	TrimF=$(ls $TrimPath/F/*.fastq*)
	TrimR=$(ls $TrimPath/R/*.fastq*)
	GenomeSz=36

	echo $Strain
	if [ "$Strain" == '675' ]; then 
		ExpCov=27
		MinCov=9
		InsLgth=600
		echo "$Strain set"
	elif [ "$Strain" == '97.0013' ]; then
		ExpCov=30
		MinCov=10
		InsLgth=600
		echo "$Strain set"
	elif [ "$Strain" == '97.0016' ]; then
		ExpCov=23
		MinCov=8
		InsLgth=600
		echo "$Strain set"
	elif [ "$Strain" == '650' ]; then
		ExpCov=27
		MinCov=9
		InsLgth=600
		echo "$Strain set"
	elif [ "$Strain" == '1082' ]; then
		ExpCov=16
		MinCov=6
		InsLgth=600
		echo "$Strain set"
	elif [ "$Strain" == '1164' ]; then
		ExpCov=24
		MinCov=8
		InsLgth=600
		echo "$Strain set"
	elif [ "$Strain" == '1166' ]; then
		ExpCov=34
		MinCov=11
		InsLgth=600
		echo "$Strain set"
	elif [ "$Strain" == '1177' ]; then
		ExpCov=50
		MinCov=17
		InsLgth=600
		echo "$Strain set"
	elif [ "$Strain" == '24350' ]; then
		ExpCov=30
		MinCov=10
		InsLgth=600
		echo "$Strain set"
	elif [ "$Strain" == '635' ]; then
		ExpCov=21
		MinCov=7
		InsLgth=600
		echo "$Strain set"
	elif [ "$Strain" == '648' ]; then
		ExpCov=43
		MinCov=14
		InsLgth=600
		echo "$Strain set"
	elif [ "$Strain" == '743' ]; then
		ExpCov=32
		MinCov=11
		InsLgth=600
		echo "$Strain set"
	fi
	
	qsub $ProgDir/submit_velvet_range.sh $MinHash $MaxHash $HashStep \
	$TrimF $TrimR $GenomeSz $ExpCov $MinCov $InsLgth
	done
	
```

#Repeatmasking

Repeat masking was performed and used the following programs:
	Repeatmasker
	Repeatmodeler
	
	

