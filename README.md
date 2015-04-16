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

Data quality was visualised once again following trimming:
```shell
	for RawData in qc_dna/paired/*/*/*/*.fastq*; do 
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc
		echo $RawData; 
		qsub $ProgDir/run_fastqc.sh $RawData
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


Assemblies were summarised to allow the best assembly to be determined by eye

```shell
	for StrainPath in $(ls -d assembly/velvet/A.alternata_ssp._*/*); do
		printf "N50\tMax_contig_size\tNumber of bases in contigs\tNumber of contigs\tNumber of contigs >=1kb\tNumber of contigs in N50\tNumber of bases in contigs >=1kb\tGC Content of contigs\n" > $StrainPath/assembly_stats.csv
		for StatsFile in $(ls $StrainPath/*/stats.txt); do 
			cat $StatsFile | rev | cut -f1 -d ' ' | rev | paste -d '\t' -s >> $StrainPath/assembly_stats.csv
		done
	done
	tail -n+1 assembly/velvet/A.alternata_ssp._*/*/assembly_stats.csv > assembly/velvet/alternaria_assembly_stats.csv
```


#Repeatmasking

Repeat masking was performed and used the following programs:
	Repeatmasker
	Repeatmodeler

The best assemblies were used to perform repeatmasking
	
```shell
	ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/repeat_masking
	BestAss675=assembly/velvet/A.alternata_ssp._arborescens/675/A.alternata_ssp._arborescens_675_69/sorted_contigs.fa
	BestAss970013=assembly/velvet/A.alternata_ssp._arborescens/97.0013/A.alternata_ssp._arborescens_97.0013_59/sorted_contigs.fa
	BestAss970016=assembly/velvet/A.alternata_ssp._arborescens/97.0016/A.alternata_ssp._arborescens_97.0016_77/sorted_contigs.fa
	BestAss650=assembly/velvet/A.alternata_ssp._gaisen/650/A.alternata_ssp._gaisen_650_67/sorted_contigs.fa
	BestAss1082=assembly/velvet/A.alternata_ssp._tenuissima/1082/A.alternata_ssp._tenuissima_1082_49/sorted_contigs.fa
	BestAss1164=assembly/velvet/A.alternata_ssp._tenuissima/1164/A.alternata_ssp._tenuissima_1164_67/sorted_contigs.fa
	BestAss1166=assembly/velvet/A.alternata_ssp._tenuissima/1166/A.alternata_ssp._tenuissima_1166_43/sorted_contigs.fa
	BestAss1177=assembly/velvet/A.alternata_ssp._tenuissima/1177/A.alternata_ssp._tenuissima_1177_63/sorted_contigs.fa
	BestAss24350=assembly/velvet/A.alternata_ssp._tenuissima/24350/A.alternata_ssp._tenuissima_24350_63/sorted_contigs.fa
	BestAss635=assembly/velvet/A.alternata_ssp._tenuissima/635/A.alternata_ssp._tenuissima_635_59/sorted_contigs.fa
	BestAss648=assembly/velvet/A.alternata_ssp._tenuissima/648/A.alternata_ssp._tenuissima_648_45/sorted_contigs.fa
	BestAss743=assembly/velvet/A.alternata_ssp._tenuissima/743/A.alternata_ssp._tenuissima_743_69/sorted_contigs.fa

	qsub $ProgDir/rep_modeling.sh $BestAss675
	qsub $ProgDir/rep_modeling.sh $BestAss970013
	qsub $ProgDir/rep_modeling.sh $BestAss970016
	qsub $ProgDir/rep_modeling.sh $BestAss650
	qsub $ProgDir/rep_modeling.sh $BestAss1082
	qsub $ProgDir/rep_modeling.sh $BestAss1164
	qsub $ProgDir/rep_modeling.sh $BestAss1166
	qsub $ProgDir/rep_modeling.sh $BestAss1177
	qsub $ProgDir/rep_modeling.sh $BestAss24350
	qsub $ProgDir/rep_modeling.sh $BestAss635
	qsub $ProgDir/rep_modeling.sh $BestAss648
	qsub $ProgDir/rep_modeling.sh $BestAss743

	qsub $ProgDir/transposonPSI.sh $BestAss675
	qsub $ProgDir/transposonPSI.sh $BestAss970013
	qsub $ProgDir/transposonPSI.sh $BestAss970016
	qsub $ProgDir/transposonPSI.sh $BestAss650
	qsub $ProgDir/transposonPSI.sh $BestAss1082
	qsub $ProgDir/transposonPSI.sh $BestAss1164
	qsub $ProgDir/transposonPSI.sh $BestAss1166
	qsub $ProgDir/transposonPSI.sh $BestAss1177
	qsub $ProgDir/transposonPSI.sh $BestAss24350
	qsub $ProgDir/transposonPSI.sh $BestAss635
	qsub $ProgDir/transposonPSI.sh $BestAss648
	qsub $ProgDir/transposonPSI.sh $BestAss743
```
	
#Gene Prediction

Gene prediction followed two steps:
	Gene model training
		- Gene models were trained for isolates 1166 and 650 using assembled RNAseq data
	Gene prediction
		- Gene models were used to predict genes in A. alternata genomes. This used RNAseq data as hints for gene models.

#Gene model training

Data quality was visualised using fastqc:
```shell
	for RawData in raw_rna/paired/*/*/*/*.fastq.gz; do 
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc
		echo $RawData; 
		qsub $ProgDir/run_fastqc.sh $RawData
	done
```

Trimming was performed on data to trim adapters from 
sequences and remove poor quality data. This was done with fastq-mcf

```shell
	for StrainPath in raw_rna/paired/*/*; do 
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/rna_qc
		IlluminaAdapters=/home/armita/git_repos/emr_repos/tools/seq_tools/illumina_full_adapters.fa
		ReadsF=$(ls $StrainPath/F/*.fastq.gz)
		ReadsR=$(ls $StrainPath/R/*.fastq.gz)
		echo $ReadsF
		echo $ReadsR
		qsub $ProgDir/rna_qc_fastq-mcf.sh $ReadsF $ReadsR $IlluminaAdapters RNA
	done
```

Data quality was visualised once again following trimming:
```shell
	for TrimData in qc_rna/paired/*/*/*/*.fastq.gz; do 
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc
		echo $RawData; 
		qsub $ProgDir/run_fastqc.sh $TrimData
	done
```	

RNAseq data was assembled into transcriptomes using Trinity
```shell
	for StrainPath in qc_rna/paired/*/*; do
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/transcriptome_assembly
		ReadsF=$(ls $StrainPath/F/*.fastq.gz)
		ReadsR=$(ls $StrainPath/R/*.fastq.gz)	
		echo $ReadsF
		echo $ReadsR
		qsub $ProgDir/transcriptome_assembly_trinity.sh $ReadsF $ReadsR
	done
```	
	
Gene prediction



/home/armita/git_repos/emr_repos/tools/gene_prediction/augustus/training_by_transcriptome.sh assembly/trinity/A.alternata_ssp._gaisen/650/650_rna_contigs/Trinity.fasta assembly/velvet/A.alternata_ssp._gaisen/650/A.alternata_ssp._gaisen_650_53/sorted_contigs.fa 


