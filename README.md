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


Gene prediction followed three steps:
	Pre-gene prediction
		- Quality of genome assemblies were assessed using Cegma to see how many core eukaryotic genes can be identified. 
	Gene model training
		- Gene models were trained for isolates 1166 and 650 using assembled RNAseq data
	Gene prediction
		- Gene models were used to predict genes in A. alternata genomes. This used RNAseq data as hints for gene models.

#Pre-gene prediction

Quality of genome assemblies was assessed by looking for the gene space in the assemblies.
```shell
	ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/cegma
	cd /home/groups/harrisonlab/project_files/alternaria/
	for Genome in $(ls repeat_masked/A.*/*/*/*_contigs_unmasked.fa); do 
		echo $Genome; 
		qsub $ProgDir/sub_cegma.sh $Genome dna;
	done
```

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
Gene training was performed using RNAseq data. The cluster can not run this script using qlogin. As such it was run on the head node (-naughty) using screen.
Training for 650 and 1166 was performed in two instances of screen and occassionally viewed to check progress over time.
(screen is detached after opening using ctrl+a then ctrl+d. - if just ctrl+d is pressed the instance of screen is deleted. - be careful)
```shell
	screen -a
	ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/augustus
	Assembly650=assembly/trinity/A.alternata_ssp._gaisen/650/650_rna_contigs/Trinity.fasta 
	Genome650=repeat_masked/A.alternata_ssp._gaisen/650/A.alternata_ssp._gaisen_650_67_repmask/650_contigs_unmasked.fa
	$ProgDir/training_by_transcriptome.sh $Assembly650 $Genome650

	screen -a
	ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/augustus
	Assembly1166=assembly/trinity/A.alternata_ssp._tenuissima/1166/1166_rna_contigs/Trinity.fasta 
	Genome1166=repeat_masked/A.alternata_ssp._tenuissima/1166/A.alternata_ssp._tenuissima_1166_43_repmask/1166_contigs_unmasked.fa
	$ProgDir/training_by_transcriptome.sh $Assembly1166 $Genome1166
```

Quality of Trinity assemblies were assessed using Cegma to assess gene-space within the transcriptome
```shell
	ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/cegma
	for Transcriptome in $(ls assembly/trinity/A.*/*/*_rna_contigs/Trinity.fasta); do  
		echo $Transcriptome;  
		qsub $ProgDir/sub_cegma.sh $Transcriptome rna; 
	done
```

	
#Gene prediction

Gene prediction was performed for A. alternata isolates. 
RNAseq reads were used as Hints for the location of CDS. 
A concatenated dataset of both ssp. tenuissima and ssp. gaisen RNAseq reads were used as hints for all strains.
Genes were predicted for ssp. tenuissima using the gene model trained to ssp. tenunissima. 
Genes were predicted for ssp. gaisen using the gene model trained to ssp. gaisen. 
Genes were predicted for ssp. arborescens using the gene model trained to ssp. tenuisima.


```shell
	ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/augustus
	mkdir -p qc_rna/concatenated
	RnaFiles=$(ls qc_rna/paired/A.alternata_ssp._*/*/*/*.fastq.gz | paste -s -d ' ')
	ConcatRna=qc_rna/concatenated/A.alternata_RNA_650_1166.fa.gz
	cat $RnaFiles > $ConcatRna
	
	for Genome in repeat_masked/A.alternata_ssp._*/*/*/*_contigs_unmasked.fa; do
		Species=$(printf $Genome | rev | cut -f4 -d '/' | rev)
		if [ "$Species" == 'A.alternata_ssp._tenuissima' ]; then 
			GeneModel=A.alternata_ssp._tenuissima_1166
		elif [ "$Species" == 'A.alternata_ssp._gaisen' ]; then 
			GeneModel=A.alternata_ssp._gaisen_650
		elif [ "$Species" == 'A.alternata_ssp._arborescens' ]; then 
			GeneModel=A.alternata_ssp._tenuissima_1166
		fi
		qsub $ProgDir/augustus_pipe.sh $Genome $ConcatRna $GeneModel
	done
```

#Functional annotation

Interproscan was used to give gene models functional annotations.

```shell
	ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/interproscan/
	for Genes in $(ls gene_pred/augustus/A.alternata_ssp._*/*/*_augustus_preds.aa); do
		echo $Genes
		$ProgDir/sub_interproscan.sh $Genes
	done
```









#Genomic analysis


The first analysis was based upon BLAST searches for genes known to be involved in toxin production

#BLAST Searches
```shell
	ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast/
	Query=analysis/blast_homology/CDC_genes/A.alternata_CDC_genes.fa 
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
	qsub $ProgDir/blast_pipe.sh $Query dna $BestAss675
	qsub $ProgDir/blast_pipe.sh $Query dna $BestAss970013
	qsub $ProgDir/blast_pipe.sh $Query dna $BestAss970016
	qsub $ProgDir/blast_pipe.sh $Query dna $BestAss650
	qsub $ProgDir/blast_pipe.sh $Query dna $BestAss1082
	qsub $ProgDir/blast_pipe.sh $Query dna $BestAss1164
	qsub $ProgDir/blast_pipe.sh $Query dna $BestAss1166
	qsub $ProgDir/blast_pipe.sh $Query dna $BestAss1177
	qsub $ProgDir/blast_pipe.sh $Query dna $BestAss24350
	qsub $ProgDir/blast_pipe.sh $Query dna $BestAss635
	qsub $ProgDir/blast_pipe.sh $Query dna $BestAss648
	qsub $ProgDir/blast_pipe.sh $Query dna $BestAss743
```
BLAST search results were summarised into a presence/absence table

The presence /absence table determines presence if a hit is present and the alignment represents >50% of the query sequence.
This thresholding means that some hits have not been summarised including AMT11, AMT15 and ALT1.
```shell
	ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast/
	InFiles=$(ls analysis/blast_homology/A.alternata_ssp._*/*/*_A.alternata_CDC_genes.fa_homologs.csv | paste -s -d ' ')
	echo $InFiles
	$ProgDir/blast_differentials.pl $InFiles
	mv *.csv analysis/blast_homology/CDC_genes/.
``` 






#CDC Assembly

Raw reads were aligned against assembled genomes to identify contigs that were unique to a isolate or clade
```shell
	for Pathz in $(ls -d qc_dna/paired/A.alternata_ssp._*/*); do  
		Strain=$(echo $Pathz | cut -d '/' -f4)
		echo "using reads for $Strain"
		ProgPath=/home/armita/git_repos/emr_repos/tools/pathogen/lineage_specific_regions
		F_IN=$(ls $Pathz/F/*.fastq.gz)
		R_IN=$(ls $Pathz/R/*.fastq.gz)
		for Assemblyz in $(ls repeat_masked/A.alternata_ssp._*/*/*/*_contigs_unmasked.fa); do
			basename $Assemblyz
			qsub "$ProgPath"/bowtie2_alignment_pipe.sh $F_IN $R_IN $Assemblyz
		done
	done
```
A summary file was made from the alignment logs.
The percentage of reads aligning to each set of assembled contigs was determined.

```shell
	SummaryFile=analysis/ls_contigs/alignment_summaries.txt
	printf "" > "$SummaryFile"
	for OUTPUT in $(ls bowtie2_alignment_pipe.sh.e*); do 
		ID=$(echo $OUTPUT | rev | cut -d 'e' -f1 | rev | less); 
		cat bowtie2_alignment_pipe.sh.o"$ID" | grep -E "Trimmed .* reads .*/F/|Output files: " | sed -e 's/.*\/F\///g' | cut -f1 -d ')' | cut -f2 -d '"' >> "$SummaryFile"; 
		cat $OUTPUT >> "$SummaryFile"; 
		printf "\n" >> "$SummaryFile"; 
	done
	cat $SummaryFile | grep ' overall alignment rate' | cut -f1 -d ' ' | less
```
These are the reads percentages reported by bowtie but do not actually reflect the percentage unaligned reads where neither pair matched.

These were identified using SAM FLAGS to extract unaligned pairs.
```shell
	for Pathz in $(ls assembly/ls_contigs/A.alternata_ssp._*/*/vs_*/*_sorted.bam); do  
		OutFileF=$(echo $Pathz | sed 's/.bam/_unaligned_F.txt/g')
		OutFileR=$(echo $Pathz | sed 's/.bam/_unaligned_R.txt/g')
		OutFileSum=$(echo $Pathz | sed 's/.bam/_sum.txt/g')
		samtools view -f 77 "$Pathz" | cut -f1 > $OutFileF
		samtools view -f 141 "$Pathz" | cut -f1 > $OutFileR
		NoReads=$(samtools view -f 1 "$Pathz" | wc -l)
		NoReadsF=$(cat "$OutFileF" | wc -l)
		NoReadsR=$(cat "$OutFileR" | wc -l)
		printf "File\tNo. paired reads\tF reads unaligned\tR reads unaligned\n" > "$OutFileSum"
		printf "$Pathz\t$NoReads\t$NoReadsF\t$NoReadsR\n" >> "$OutFileSum"
	done

	for File in $(ls assembly/ls_contigs/A.alternata_ssp._*/*/vs_*/*_sum.txt); do 
		cat $File | tail -n+2; 
	done > analysis/ls_contigs/assembly_summaries2.txt 
```

The number of reads aligning per bp of assembly was determined. Typical alignment values were 0.20 reads per bp. Contigs were detemined as unique to that alignment if they contained an average of 0 reads per bp.
The number of bp unique to each assembly were identified.

```shell
	mkdir -p analysis/ls_contigs
	Outfile=analysis/ls_contigs/ls_contig_size.csv
	printf "Reads" > $Outfile
	for Genome in $(ls -d assembly/ls_contigs/A.alternata_ssp._arborescens/675/vs_A.alternata_ssp._*); do
		NameG=$(basename "$Genome" | sed 's/vs_//g' | sed 's/_repmask_contigs//g')
		printf "\t$NameG" >> $Outfile
	done
	printf "\n" >> $Outfile
	for Reads in $(ls -d assembly/ls_contigs/A.alternata_ssp._*/*); do
		NameR=$(basename "$Reads")
		printf "$NameR" >> $Outfile
		for Genome in $(ls -d $Reads/*); do
			printf "\t" >> $Outfile
			cat $Genome/*_sorted_indexstats_coverage.csv | head -n-1 | grep -E -w "0.00$" | cut -f2 | awk '{s+=$1} END {printf s}' >> $Outfile
		done
		printf "\n" >> $Outfile
	done
```
This did not give clear results.

Contigs that had no reads align to them were identified.
