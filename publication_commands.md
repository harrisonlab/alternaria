# fusarium_ex_strawberry
Commands used in the analysis of Alternaria spp.isolates ex. strawberry.

Alternaria alternata ssp.
====================

Commands used during analysis of the Alternaria alternata ssp. genome. Note - all this work was performed in the directory: /home/groups/harrisonlab/project_files/alternaria

The following is a summary of the work presented in this Readme:
Data organisation:
  * Preparing data  
Draft Genome assembly
  * Data qc
  * Genome assembly
  * Repeatmasking
  * Gene prediction
  * Functional annotation
Genome analysis
  * Homology between predicted genes & published effectors


#Data organisation

Data was copied from the raw_data repository to a local directory for assembly
and annotation.

```bash
    mkdir -p /home/groups/harrisonlab/project_files/alternaria
  	cd /home/groups/harrisonlab/project_files/alternaria
  	Species=F.oxysporum_fsp_fragariae
  	Strain=FeChina
    mkdir -p raw_dna/paired/fusarium_ex_strawberry/FeChina/F
    mkdir -p raw_dna/paired/fusarium_ex_strawberry/FeChina/R
    RawDat=/home/groups/harrisonlab/raw_data/raw_seq/raw_reads/150925_M01678_0029_AC669
    cp $RawDat/FeChina_S1_L001_R1_001.fastq.gz raw_dna/paired/fusarium_ex_strawberry/FeChina/F/.
    cp $RawDat/FeChina_S1_L001_R1_001.fastq.gz raw_dna/paired/fusarium_ex_strawberry/FeChina/R/.
```


#Data qc

programs: fastqc fastq-mcf kmc

Data quality was visualised using fastqc:

```bash
  	for RawData in $(ls raw_dna/paired/*/*/*/*.fastq.gz); do
	  echo $RawData;
	  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc;
	  qsub $ProgDir/run_fastqc.sh $RawData;
  	done
```

Trimming was performed on data to trim adapters from sequences and remove poor quality data.
This was done with fastq-mcf


```bash
	Read_F=raw_dna/paired/fusarium_ex_strawberry/FeChina/F/FeChina_S1_L001_R1_001.fastq.gz
	Read_R=raw_dna/paired/fusarium_ex_strawberry/FeChina/R/FeChina_S1_L001_R1_001.fastq.gz
	IluminaAdapters=/home/armita/git_repos/emr_repos/tools/seq_tools/illumina_full_adapters.fa
	ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/rna_qc
	qsub $ProgDir/rna_qc_fastq-mcf.sh $Read_F $Read_R $IluminaAdapters DNA
```


Data quality was visualised once again following trimming:

```bash
	for RawData in $(ls qc_dna/paired/*/*/*/*.fq.gz); do
  	echo $RawData;
  	ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc;
  	qsub $ProgDir/run_fastqc.sh $RawData;
  done
```


kmer counting was performed using kmc.
This allowed estimation of sequencing depth and total genome size:

```bash
	Trim_F=qc_dna/paired/fusarium_ex_strawberry/FeChina/F/FeChina_S1_L001_R1_001_trim.fq.gz
	Trim_R=qc_dna/paired/fusarium_ex_strawberry/FeChina/R/FeChina_S1_L001_R1_001_trim.fq.gz
	ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc
	qsub $ProgDir/kmc_kmer_counting.sh $Trim_F $Trim_R
```

** Estimated Genome Size is: 48490129

** Esimated Coverage is: 35

# Assembly
Assembly was performed using: Velvet / Abyss / Spades

## Velvet Assembly
A range of hash lengths were used and the best assembly selected for subsequent analysis

```bash

  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/velvet
  MinHash=41
  MaxHash=81
  HashStep=2
  Trim_F=qc_dna/paired/fusarium_ex_strawberry/FeChina/F/FeChina_S1_L001_R1_001_trim.fq.gz
	Trim_R=qc_dna/paired/fusarium_ex_strawberry/FeChina/R/FeChina_S1_L001_R1_001_trim.fq.gz
  GenomeSz=36
  echo $Strain
  ExpCov=35
  MinCov=10
  InsLgth=600
  qsub $ProgDir/submit_velvet_range.sh \
  $MinHash $MaxHash $HashStep $Trim_F $Trim_R $GenomeSz $ExpCov $MinCov $InsLgth
```

## Spades Assembly
<!--
First run error correction. (This job is CPU intensive rather than RAM intensive
and will run on any node of the cluster).

```bash
	F_Read=qc_dna/paired/fusarium_ex_strawberry/FeChina/F/FeChina_S1_L001_R1_001_trim.fq.gz
	R_Read=qc_dna/paired/fusarium_ex_strawberry/FeChina/R/FeChina_S1_L001_R1_001_trim.fq.gz
	ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/spades
	OutDir=qc_dna/paired/fusarium_ex_strawberry/FeChina/corrected
  qsub $ProgDir/sub_spades_correction.sh $F_Read $R_Read $OutDir
``` -->


```bash
  F_Read=qc_dna/paired/fusarium_ex_strawberry/FeChina/F/FeChina_S1_L001_R1_001_trim.fq.gz
  R_Read=qc_dna/paired/fusarium_ex_strawberry/FeChina/R/FeChina_S1_L001_R1_001_trim.fq.gz
	ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/spades
	OutDir=assembly/spades/usarium_ex_strawberry/FeChina
  qsub $ProgDir/submit_SPAdes.sh $F_Read $R_Read $OutDir correct
```
<!--
Quast

```bash
	ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
	Assembly=assembly/spades/N.ditissima/R0905_v2/filtered_contigs/contigs_min_500bp.fasta
	OutDir=assembly/spades/N.ditissima/R0905_v2/filtered_contigs
	qsub $ProgDir/sub_quast.sh $Assembly $OutDir
```

Assemblies were summarised to allow the best assembly to be determined by eye.

** Assembly stats are:
  * Assembly size:
  * N50:153669
  * N80:
  * N20:
  * Longest contig:687738
  **

As SPADes was run with the option to autodetect a minimum coverage the assembly was assessed to identify the coverage of assembled contigs. This was done using the following command:

	BestAss=assembly/spades/N.ditissima/R0905_v2/filtered_contigs/contigs_min_500bp.fasta
	cat $BestAss | grep '>' | cut -f6 -d'_' | sort -n | cut -f1 -d '.' | sort -n | uniq -c | less

From this it was determined that SPades could not be trusted to set its own minimum threshold for coverage.
In future an option will be be used to set a coverage for spades.
In the meantime contigs with a coverage lower than 10 were filtered out using the following commands:

	Headers=assembly/spades/N.ditissima/R0905_v2/filtered_contigs/contigs_min_500bp_10x_headers.txt
	cat $BestAss | grep '>' | grep -E -v 'cov_.\..*_' > $Headers
	FastaMinCov=assembly/spades/N.ditissima/R0905_v2/filtered_contigs/contigs_min_500bp_10x_headers.fasta
	cat $BestAss | sed -e 's/\(^>.*$\)/#\1#/' | tr -d "\r" | tr -d "\n" | sed -e 's/$/#/' | tr "#" "\n" | sed -e '/^$/d' | grep -A1 -f $Headers | grep -v -E '^\-\-' > $FastaMinCov

```bash
	~/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants/remove_contaminants.py --inp ../neonectria_ditissima/assembly/spades/N.ditissima/R0905_v2/filtered_contigs/contigs_min_500bp_10x_headers.fasta  --out assembly/spades/N.galligena/R0905_v2/filtered_contigs/contigs_min_500bp_10x_filtered_renamed.fasta  --coord_file editfile.tab

We run Quast again.

	ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
	Assembly=assembly/spades/N.ditissima/R0905_v2/filtered_contigs/contigs_min_500bp_10x_headers.fasta
	OutDir=assembly/spades/N.ditissima/R0905_v2/contigs_min_500bp_10x_headers
	qsub $ProgDir/sub_quast.sh $Assembly $OutDir

# Repeat masking
Repeat masking was performed and used the following programs: Repeatmasker Repeatmodeler

The best assembly was used to perform repeatmasking

```bash
	ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/repeat_masking
	BestAss=/assembly/spades/N.ditissima/R0905_v2/filtered_contigs/contigs_min_500bp_10x_headers.fasta
	qsub $ProgDir/rep_modeling.sh $BestAss
	qsub $ProgDir/transposonPSI.sh $BestAss
 ```

** % bases masked by repeatmasker: 11.73%

** % bases masked by transposon psi: **


# Gene Prediction
Gene prediction followed two steps:
Pre-gene prediction - Quality of genome assemblies were assessed using Cegma to see how many core eukaryotic genes can be identified.
Gene models were used to predict genes in the Neonectria genome. This used results from CEGMA as hints for gene models.

## Pre-gene prediction
Quality of genome assemblies was assessed by looking for the gene space in the assemblies.

```bash
  	ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/cegma
  	Assembly=/repeat_masked/spades/N.ditissima/R0905_v2/filtered_contigs_repmask/R0905_v2_contigs_unmasked.fa
  	qsub $ProgDir/sub_cegma.sh $Assembly dna
```

** Number of cegma genes present and complete: 95.16
** Number of cegma genes present and partial: 97.18

##Gene prediction

Gene prediction was performed for the neonectria genome.
CEGMA genes were used as Hints for the location of CDS.

```bash
	ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/augustus
  	Assembly=/repeat_masked/spades/N.ditissima/R0905_v2/filtered_contigs_repmask/R0905_v2_contigs_unmasked.fa
  	GeneModel=fusarium
  	qsub $ProgDir/submit_augustus.sh $GeneModel $Assembly false
```

** Number of genes predicted: 12712

#Functional annotation

Interproscan was used to give gene models functional annotations.

```bash
	ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/interproscan/
  	Genes=gene_pred/augustus/spades/N.ditissima/N.ditissima_aug_out.aa
  	$ProgDir/sub_interproscan.sh $Genes
```

```bash
	ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/interproscan
	Genes=gene_pred/augustus/spades/N.ditissima/N.ditissima_aug_out.aa
	InterProRaw=gene_pred/interproscan/spades/N.ditissima/raw
	ProgDir/append_interpro.sh $Genes $InterProRaw
```

#Genomic analysis
The first analysis was based upon BLAST searches for genes known to be involved in toxin production


##Genes with homology to PHIbase
Predicted gene models were searched against the PHIbase database using tBLASTx.

```bash
	ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
	Query=../../phibase/v3.8/PHI_accessions.fa
	Subject=repeat_masked/spades/N.ditissima/NG-R0905_repmask/N.ditissima_contigs_unmasked.fa
	qsub $ProgDir/blast_pipe.sh $Query protein $Subject
```

Top BLAST hits were used to annotate gene models.

```bash

```

** Blast results of note: **
  * 'Result A'
-->
