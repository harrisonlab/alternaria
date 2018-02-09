# Minion_Assemblies
==========

This document details the commands used to assemble and annotate minion genomes.

Note - all this work was performed in the directory:
/home/groups/harrisonlab/project_files/alternaria

The following is a summary of the work presented in this Readme.

The following processes were applied to Fusarium genomes prior to analysis:
Data qc
Genome assembly
Repeatmasking
Gene prediction
Functional annotation


# 0. Building of directory structure

#Building of directory structure


Data was basecalled again using Albacore 2.02 on the minion server:

```bash
screen -a
ssh nanopore@nanopore

# upgrade albacore
wget https://mirror.oxfordnanoportal.com/software/analysis/ont_albacore-2.1.10-cp34-cp34m-manylinux1_x86_64.whl
~/.local/bin/read_fast5_basecaller.py --version
pip3 install --user ont_albacore-2.1.10-cp34-cp34m-manylinux1_x86_64.whl --upgrade
~/.local/bin/read_fast5_basecaller.py --version

mkdir Aalt_21-12-17
cd Aalt_21-12-17

# Oxford nanopore 07/03/17
Organism=A.alternata
Date=05-02-18
FlowCell="FLO-MIN106"
Kit="SQK-RBK001"
RawDatDir=/data/seq_data/minion/2017/20171221_Alternaria-AA/Alternaria-AA/GA30000/reads
OutDir=/data/scratch/nanopore_tmp_data/Alternaria/albacore_v2.1.10
mkdir -p $OutDir

mkdir -p ~/Aalt_21-12-17/$Date
cd ~/Aalt_21-12-17/$Date
~/.local/bin/read_fast5_basecaller.py \
  --flowcell $FlowCell \
  --kit $Kit \
  --input $RawDatDir \
  --recursive \
  --worker_threads 24 \
  --save_path Alt_albacore_v2.10_demultiplexed \
  --output_format fastq,fast5 \
  --reads_per_fastq_batch 4000 \
  --barcoding
  cat Alt_albacore_v2.10_demultiplexed/workspace/pass/barcode01/*.fastq | gzip -cf > Alt_albacore_v2.10_barcode01.fastq.gz
  cat Alt_albacore_v2.10_demultiplexed/workspace/pass/barcode02/*.fastq | gzip -cf > Alt_albacore_v2.10_barcode02.fastq.gz
  OutDir=/data/scratch/nanopore_tmp_data/Alternaria/albacore_v2.1.10
  mkdir -p $OutDir
  chmod +w $OutDir
  cp Alt_albacore_v2.10_barcode0*.fastq.gz $OutDir/.
  chmod +rw $OutDir/Alt_albacore_v2.10_barcode0*.fastq.gz
  # scp Alt_albacore_v2.10_barcode*.fastq.gz armita@192.168.1.200:$OutDir/.
  tar -cz -f Alt_albacore_v2.10_demultiplexed.tar.gz Alt_albacore_v2.10_demultiplexed
  OutDir=/data/scratch/nanopore_tmp_data/Alternaria/albacore_v2.1.10
  mv Alt_albacore_v2.10_demultiplexed.tar.gz $OutDir/.
  chmod +rw $OutDir/Alt_albacore_v2.10_demultiplexed.tar.gz
```

BUilding a director structure:

```bash
ProjDir=/home/groups/harrisonlab/project_files/alternaria
cd $ProjDir

Organism=A.alternata_ssp_tenuissima
Strain=1166
OutDir=raw_dna/minion/$Organism/$Strain
mkdir -p $OutDir
RawDat=$(ls /data/scratch/nanopore_tmp_data/Alternaria/albacore_v2.1.10/Alt_albacore_v2.10_barcode01.fastq.gz)
cd $OutDir
cp -s $RawDat .
cd $ProjDir

Organism=A.alternata_ssp_arborescens
Strain=650
OutDir=raw_dna/minion/$Organism/$Strain
mkdir -p $OutDir
RawDat=$(ls /data/scratch/nanopore_tmp_data/Alternaria/albacore_v2.1.10/Alt_albacore_v2.10_barcode02.fastq.gz)
cd $OutDir
cp -s $RawDat .
cd $ProjDir

```


<!--
### MiSeq data

```bash
	RawDatDir=/home/miseq_data/2017/RAW/170626_M04465_0043_000000000-B48RG/Data/Intensities/BaseCalls
	ProjectDir=/home/groups/harrisonlab/project_files/fusarium
	OutDir=$ProjectDir/raw_dna/paired/F.oxysporum_fsp_narcissi/FON_63
	mkdir -p $OutDir/F
	mkdir -p $OutDir/R
	cp $RawDatDir/FON63_S2_L001_R1_001.fastq.gz $OutDir/F/.
	cp $RawDatDir/FON63_S2_L001_R2_001.fastq.gz $OutDir/R/.
```


#### QC of MiSeq data

programs:
  fastqc
  fastq-mcf
  kmc

Data quality was visualised using fastqc:
```bash
	for RawData in $(ls raw_dna/paired/*/*/*/*.fastq.gz | grep 'FON_63'); do
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc
		echo $RawData;
		qsub $ProgDir/run_fastqc.sh $RawData
	done
```

```bash
	for StrainPath in $(ls -d raw_dna/paired/*/* | grep 'FON_63'); do
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/rna_qc
		IlluminaAdapters=/home/armita/git_repos/emr_repos/tools/seq_tools/ncbi_adapters.fa
		ReadsF=$(ls $StrainPath/F/*.fastq*)
		ReadsR=$(ls $StrainPath/R/*.fastq*)
		echo $ReadsF
		echo $ReadsR
		qsub $ProgDir/rna_qc_fastq-mcf.sh $ReadsF $ReadsR $IlluminaAdapters DNA
	done
``` -->

## Assembly

### Removal of adapters

Splitting reads and trimming adapters using porechop
```bash
	for RawReads in $(ls raw_dna/minion/F.oxysporum_fsp_narcissi/FON_63/all_reads_albacore1.1.1.fastq.gz); do
    Organism=$(echo $RawReads| rev | cut -f3 -d '/' | rev)
    Strain=$(echo $RawReads | rev | cut -f2 -d '/' | rev)
    echo "$Organism - $Strain"
  	OutDir=qc_dna/minion/$Organism/$Strain
  	ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc
  	qsub $ProgDir/sub_porechop.sh $RawReads $OutDir
  done
```


## Identify sequencing coverage

For Minion data:
```bash
	for RawData in $(ls qc_dna/minion/*/*/*q.gz | grep 'FON_63'); do
		echo $RawData;
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc;
		GenomeSz=60
		OutDir=$(dirname $RawData)
		qsub $ProgDir/sub_count_nuc.sh $GenomeSz $RawData $OutDir
	done
```

```bash
  for StrainDir in $(ls -d qc_dna/minion/*/* | grep 'FON_63'); do
    Strain=$(basename $StrainDir)
    printf "$Strain\t"
    for File in $(ls $StrainDir/*.txt); do
      echo $(basename $File);
      cat $File | tail -n1 | rev | cut -f2 -d ' ' | rev;
    done | grep -v '.txt' | awk '{ SUM += $1} END { print SUM }'
  done
```
MinION coverage
```
FON_63	122.23
```

For Miseq data:
```bash
	for RawData in $(ls qc_dna/paired/*/*/*/*q.gz | grep 'FON_63'); do
		echo $RawData;
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc;
		qsub $ProgDir/run_fastqc.sh $RawData;
		GenomeSz=60
		OutDir=$(dirname $RawData)
		qsub $ProgDir/sub_count_nuc.sh $GenomeSz $RawData $OutDir
	done
```

```bash
  for StrainDir in $(ls -d qc_dna/paired/*/* | grep 'FON_63'); do
    Strain=$(basename $StrainDir)
    printf "$Strain\t"
    for File in $(ls qc_dna/paired/*/"$Strain"/*/*.txt); do
      echo $(basename $File);
      cat $File | tail -n1 | rev | cut -f2 -d ' ' | rev;
    done | grep -v '.txt' | awk '{ SUM += $1} END { print SUM }'
  done
```
```
FON_63	52.39
```


### Read correction using Canu

```bash
  for TrimReads in $(ls qc_dna/minion/F.*/*/*_trim.fastq.gz | grep 'FON_63'); do
	# for TrimReads in $(ls assembly/canu-1.5/F.oxysporum_fsp_narcissi/FON_63/FON_63.correctedReads.fasta.gz); do
    Organism=$(echo $TrimReads | rev | cut -f3 -d '/' | rev)
    Strain=$(echo $TrimReads | rev | cut -f2 -d '/' | rev)
    # OutDir=assembly/canu-1.5/$Organism/"$Strain"
    OutDir=assembly/canu-1.6/$Organism/"$Strain"_test
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/canu
    qsub $ProgDir/sub_canu_correction.sh $TrimReads 60m $Strain $OutDir
  done
```

### Assembbly using SMARTdenovo

```bash
for CorrectedReads in $(ls assembly/canu-1.5/F.*/*/*.trimmedReads.fasta.gz | grep 'FON_63'); do
Organism=$(echo $CorrectedReads | rev | cut -f3 -d '/' | rev)
Strain=$(echo $CorrectedReads | rev | cut -f2 -d '/' | rev)
Prefix="$Strain"_smartdenovo
OutDir=assembly/SMARTdenovo/$Organism/"$Strain"
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/SMARTdenovo
qsub $ProgDir/sub_SMARTdenovo.sh $CorrectedReads $Prefix $OutDir
done
```


Quast and busco were run to assess the effects of racon on assembly quality:

```bash
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
for Assembly in $(ls assembly/SMARTdenovo/*/*/*.dmo.lay.utg | grep 'FON_63'); do
  Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
  Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)  
  OutDir=$(dirname $Assembly)
  qsub $ProgDir/sub_quast.sh $Assembly $OutDir
done
```


```bash
for Assembly in $(ls assembly/SMARTdenovo/*/*/*.dmo.lay.utg | grep 'FON_63'); do
Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
echo "$Organism - $Strain"
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
OutDir=gene_pred/busco/$Organism/$Strain/assembly
# OutDir=$(dirname $Assembly)
qsub $ProgDir/sub_busco2.sh $Assembly $BuscoDB $OutDir
done
```


Error correction using racon:

```bash
for Assembly in $(ls assembly/SMARTdenovo/*/*/*.dmo.lay.utg | grep 'FON_63'); do
Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
echo "$Organism - $Strain"
ReadsFq=$(ls qc_dna/minion/*/$Strain/all_reads_albacore1_trim.fastq.gz)
Iterations=10
OutDir=$(dirname $Assembly)"/racon_$Iterations"
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/racon
qsub $ProgDir/sub_racon.sh $Assembly $ReadsFq $Iterations $OutDir
done
```

```bash
	ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
	for Assembly in $(ls assembly/SMARTdenovo/F.*/*/racon_10/*.fasta | grep 'FON_63' |  grep 'round_10'); do
		OutDir=$(dirname $Assembly)
		echo "" > tmp.txt
		ProgDir=~/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
		$ProgDir/remove_contaminants.py --keep_mitochondria --inp $Assembly --out $OutDir/racon_min_500bp_renamed.fasta --coord_file tmp.txt > $OutDir/log.txt
	done
```

Quast and busco were run to assess the effects of racon on assembly quality:

```bash
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
for Assembly in $(ls assembly/SMARTdenovo/F.*/*/racon*/racon_min_500bp_renamed.fasta | grep 'FON_63'); do
  Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
  Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
  OutDir=$(dirname $Assembly)
  qsub $ProgDir/sub_quast.sh $Assembly $OutDir
done
```


```bash
# for Assembly in $(ls assembly/SMARTdenovo/F.*/*/racon*/*.fasta | grep 'FON_63' | grep 'racon_min_500bp_renamed'); do
for Assembly in $(ls assembly/SMARTdenovo/F.*/*/racon*/*.fasta | grep 'FON_63'); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
OutDir=gene_pred/busco/$Organism/$Strain/assembly
# OutDir=$(dirname $Assembly)
qsub $ProgDir/sub_busco2.sh $Assembly $BuscoDB $OutDir
done
```
```bash
printf "Filename\tComplete\tDuplicated\tFragmented\tMissing\tTotal\n"
for File in $(ls gene_pred/busco/F*/*/assembly/*/short_summary_*.txt | grep 'FON_63'); do
FileName=$(basename $File)
Complete=$(cat $File | grep "(C)" | cut -f2)
Duplicated=$(cat $File | grep "(D)" | cut -f2)
Fragmented=$(cat $File | grep "(F)" | cut -f2)
Missing=$(cat $File | grep "(M)" | cut -f2)
Total=$(cat $File | grep "Total" | cut -f2)
printf "$FileName\t$Complete\t$Duplicated\t$Fragmented\t$Missing\t$Total\n"
done
```


# Assembly correction using nanopolish

Fast5 files are very large and need to be stored as gzipped tarballs. These needed temporarily unpacking but must be deleted after nanpolish has finished running.

```bash
	TarBall=$(ls /home/miseq_data/minion/2017/FON-63_2017-05-22/FON-63_2017-05-22_albacore_output_1.1.1_fast5s.tar.gz)
	Fast5Dir=raw_dna/minion/F.oxysporum_fsp_narcissi/FON_63
	mkdir -p $Fast5Dir
	tar -zxvf $TarBall -C $Fast5Dir
```

```bash
for Assembly in $(ls assembly/SMARTdenovo/F.*/*/racon*/racon_min_500bp_renamed.fasta | grep 'FON_63'); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
# Step 1 extract reads as a .fq file which contain info on the location of the fast5 files
# Note - the full path from home must be used
Fast5Dir=$(ls -d /home/groups/harrisonlab/project_files/fusarium/raw_dna/minion/F.oxysporum_fsp_narcissi/FON_63)
ReadDir=raw_dna/nanopolish/$Organism/$Strain
if [ -d $ReadDir ]; then
	echo "reads already extracted"
else
	echo "extracting reads"
	mkdir -p $ReadDir
	CurDir=$PWD
	cd $ReadDir
	nanopolish extract -r $Fast5Dir | gzip -cf > "$Strain"_reads.fa.gz
	cd $CurDir
fi


RawReads=$(ls $ReadDir/"$Strain"_reads.fa.gz)
OutDir=$(dirname $Assembly)
mkdir -p $OutDir
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/nanopolish
# submit alignments for nanoppolish
qsub $ProgDir/sub_bwa_nanopolish.sh $Assembly $RawReads $OutDir/nanopolish
done
```

 Split the assembly into 50Kb fragments an submit each to the cluster for
 nanopolish correction

```bash
for Assembly in $(ls assembly/SMARTdenovo/F.*/*/racon*/racon_min_500bp_renamed.fasta | grep 'FON_63'); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
OutDir=$(dirname $Assembly)
RawReads=$(ls raw_dna/nanopolish/$Organism/$Strain/"$Strain"_reads.fa.gz)
AlignedReads=$(ls $OutDir/nanopolish/reads.sorted.bam)

NanoPolishDir=/home/armita/prog/nanopolish/nanopolish/scripts
python $NanoPolishDir/nanopolish_makerange.py $Assembly > $OutDir/nanopolish/nanopolish_range.txt

Ploidy=1
echo "nanopolish log:" > nanopolish_log.txt
for Region in $(cat $OutDir/nanopolish/nanopolish_range.txt | tail -n+21); do
Jobs=$(qstat | grep 'sub_nanopo' | grep 'qw' | wc -l)
while [ $Jobs -gt 1 ]; do
sleep 1m
printf "."
Jobs=$(qstat | grep 'sub_nanopo' | grep 'qw' | wc -l)
done		
printf "\n"
echo $Region
echo $Region >> nanopolish_log.txt
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/nanopolish
qsub $ProgDir/sub_nanopolish_variants.sh $Assembly $RawReads $AlignedReads $Ploidy $Region $OutDir/$Region
done
done
```

```bash
Assembly=$(ls assembly/SMARTdenovo/F.oxysporum_fsp_narcissi/FON_63/racon_10/racon_min_500bp_renamed.fasta)
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
OutDir=assembly/SMARTdenovo/$Organism/$Strain/nanopolish
mkdir -p $OutDir
# cat "" > $OutDir/"$Strain"_nanoplish.fa
NanoPolishDir=/home/armita/prog/nanopolish/nanopolish/scripts
python $NanoPolishDir/nanopolish_merge.py assembly/SMARTdenovo/$Organism/$Strain/racon_10/*/*.fa > $OutDir/"$Strain"_nanoplish.fa

echo "" > tmp.txt
ProgDir=~/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
$ProgDir/remove_contaminants.py --keep_mitochondria --inp $OutDir/"$Strain"_nanoplish.fa --out $OutDir/"$Strain"_nanoplish_min_500bp_renamed.fasta --coord_file tmp.txt > $OutDir/log.txt
```

Quast and busco were run to assess the effects of nanopolish on assembly quality:

```bash

for Assembly in $(ls assembly/SMARTdenovo/F.oxysporum_fsp_narcissi/FON_63/nanopolish/FON_63_nanoplish_min_500bp_renamed.fasta); do
  Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
  Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
	# Quast
  OutDir=$(dirname $Assembly)
	ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
  qsub $ProgDir/sub_quast.sh $Assembly $OutDir
	# Busco
	BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
	OutDir=gene_pred/busco/$Organism/$Strain/assembly
	ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
	qsub $ProgDir/sub_busco2.sh $Assembly $BuscoDB $OutDir
done
```


### Pilon assembly correction

Assemblies were polished using Pilon

```bash
	for Assembly in $(ls assembly/SMARTdenovo/F.oxysporum_fsp_narcissi/FON_63/nanopolish/FON_63_nanoplish_min_500bp_renamed.fasta); do
		Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
		Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
		echo "$Organism - $Strain"
		IlluminaDir=$(ls -d qc_dna/paired/*/$Strain)
		TrimF1_Read=$(ls $IlluminaDir/F/*_trim.fq.gz | head -n1)
		TrimR1_Read=$(ls $IlluminaDir/R/*_trim.fq.gz | head -n1)
		OutDir=$(dirname $Assembly)/../pilon
		Iterations=10
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/pilon
		qsub $ProgDir/sub_pilon.sh $Assembly $TrimF1_Read $TrimR1_Read $OutDir $Iterations
	done
```

Contigs were renamed
```bash
echo "" > tmp.txt
Assembly=$(ls assembly/SMARTdenovo/*/*/pilon/*.fasta | grep 'FON_63' | grep 'pilon_10')
OutDir=$(dirname $Assembly)
ProgDir=~/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
$ProgDir/remove_contaminants.py --keep_mitochondria --inp $Assembly --out $OutDir/pilon_min_500bp_renamed.fasta --coord_file tmp.txt > $OutDir/log.txt
```

Quast and busco were run to assess the effects of pilon on assembly quality:

```bash
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
for Assembly in $(ls assembly/SMARTdenovo/*/*/pilon/*.fasta | grep 'FON_63' | grep 'pilon_min_500bp_renamed.fasta'); do
  Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
  Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)  
  OutDir=$(dirname $Assembly)
  qsub $ProgDir/sub_quast.sh $Assembly $OutDir
done
```


```bash
for Assembly in $(ls assembly/SMARTdenovo/*/*/pilon/*.fasta | grep 'FON_63'); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
OutDir=gene_pred/busco/$Organism/$Strain/assembly
qsub $ProgDir/sub_busco2.sh $Assembly $BuscoDB $OutDir
done
```


```bash
printf "Filename\tComplete\tDuplicated\tFragmented\tMissing\tTotal\n"
for File in $(ls gene_pred/busco/F*/*/assembly/*/short_summary_*.txt | grep 'FON_63'); do  
FileName=$(basename $File)
Complete=$(cat $File | grep "(C)" | cut -f2)
Duplicated=$(cat $File | grep "(D)" | cut -f2)
Fragmented=$(cat $File | grep "(F)" | cut -f2)
Missing=$(cat $File | grep "(M)" | cut -f2)
Total=$(cat $File | grep "Total" | cut -f2)
printf "$FileName\t$Complete\t$Duplicated\t$Fragmented\t$Missing\t$Total\n"
done
```

```
short_summary_FON_63_nanoplish_min_500bp_renamed.txt	3569	34	74	82	3725
short_summary_FON_63_smartdenovo.dmo.lay.txt	636	2	582	2507	3725
short_summary_FON_63_smartdenovo_racon_round_10.txt	2577	18	605	543	3725
short_summary_FON_63_smartdenovo_racon_round_1.txt	2344	12	725	656	3725
short_summary_FON_63_smartdenovo_racon_round_2.txt	2528	16	620	577	3725
short_summary_FON_63_smartdenovo_racon_round_3.txt	2534	21	617	574	3725
short_summary_FON_63_smartdenovo_racon_round_4.txt	2520	19	636	569	3725
short_summary_FON_63_smartdenovo_racon_round_5.txt	2551	20	609	565	3725
short_summary_FON_63_smartdenovo_racon_round_6.txt	2543	13	635	547	3725
short_summary_FON_63_smartdenovo_racon_round_7.txt	2533	14	618	574	3725
short_summary_FON_63_smartdenovo_racon_round_8.txt	2558	19	616	551	3725
short_summary_FON_63_smartdenovo_racon_round_9.txt	2546	18	628	551	3725
short_summary_pilon_10.txt	3690	37	14	21	3725
short_summary_pilon_1.txt	3689	37	15	21	3725
short_summary_pilon_2.txt	3690	37	14	21	3725
short_summary_pilon_3.txt	3690	37	14	21	3725
short_summary_pilon_4.txt	3690	37	14	21	3725
short_summary_pilon_5.txt	3690	37	14	21	3725
short_summary_pilon_6.txt	3690	37	14	21	3725
short_summary_pilon_7.txt	3690	37	14	21	3725
short_summary_pilon_8.txt	3690	37	14	21	3725
short_summary_pilon_9.txt	3690	37	14	21	3725
short_summary_pilon_min_500bp_renamed.txt	3690	37	14	21	3725
short_summary_racon_min_500bp_renamed.txt	2577	18	605	543	3725
```

# Hybrid Assembly

<!--
## Spades Assembly

```bash
for TrimReads in $(ls qc_dna/minion/*/*/*_trim.fastq.gz | grep 'FON_63'); do
Organism=$(echo $TrimReads | rev | cut -f3 -d '/' | rev)
# Organism="F.oxysporum_fsp_mathioli"
Strain=$(echo $TrimReads | rev | cut -f2 -d '/' | rev)
IlluminaDir=$(ls -d qc_dna/paired/*/$Strain)
TrimF1_Read=$(ls $IlluminaDir/F/*_trim.fq.gz | head -n1)
TrimR1_Read=$(ls $IlluminaDir/R/*_trim.fq.gz | head -n1)
OutDir=assembly/spades_minion/$Organism/"$Strain"
echo $TrimR1_Read
echo $TrimR1_Read
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/spades
qsub $ProgDir/sub_spades_minion.sh $TrimReads $TrimF1_Read $TrimR1_Read $OutDir
done
```

Contigs shorter thaan 500bp were removed from the assembly

```bash
  for Contigs in $(ls assembly/spades_minion/*/*/contigs.fasta); do
    AssemblyDir=$(dirname $Contigs)
    mkdir $AssemblyDir/filtered_contigs
    FilterDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/abyss
    $FilterDir/filter_abyss_contigs.py $Contigs 500 > $AssemblyDir/filtered_contigs/contigs_min_500bp.fasta
  done
```
 -->

# Repeat Masking

Repeat masking was performed on the non-hybrid assembly.

```bash
	for Assembly in $(ls assembly/SMARTdenovo/*/*/pilon/pilon_min_500bp_renamed.fasta | grep 'FON_63'); do
		Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
		Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
		echo "$Organism - $Strain"
		OutDir=repeat_masked/$Organism/"$Strain"/filtered_contigs
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/repeat_masking
		qsub $ProgDir/rep_modeling.sh $Assembly $OutDir
		qsub $ProgDir/transposonPSI.sh $Assembly $OutDir
	done
```

The TransposonPSI masked bases were used to mask additional bases from the
repeatmasker / repeatmodeller softmasked and hardmasked files.

```bash

for File in $(ls repeat_masked/*/*/*/*_contigs_softmasked.fa | grep 'FON_63'); do
OutDir=$(dirname $File)
TPSI=$(ls $OutDir/*_contigs_unmasked.fa.TPSI.allHits.chains.gff3)
OutFile=$(echo $File | sed 's/_contigs_softmasked.fa/_contigs_softmasked_repeatmasker_TPSI_appended.fa/g')
echo "$OutFile"
bedtools maskfasta -soft -fi $File -bed $TPSI -fo $OutFile
echo "Number of masked bases:"
cat $OutFile | grep -v '>' | tr -d '\n' | awk '{print $0, gsub("[a-z]", ".")}' | cut -f2 -d ' '
done
# The number of N's in hardmasked sequence are not counted as some may be present within the assembly and were therefore not repeatmasked.
for File in $(ls repeat_masked/*/*/*/*_contigs_hardmasked.fa | grep 'FON_63'); do
OutDir=$(dirname $File)
TPSI=$(ls $OutDir/*_contigs_unmasked.fa.TPSI.allHits.chains.gff3)
OutFile=$(echo $File | sed 's/_contigs_hardmasked.fa/_contigs_hardmasked_repeatmasker_TPSI_appended.fa/g')
echo "$OutFile"
bedtools maskfasta -fi $File -bed $TPSI -fo $OutFile
done
```
```
Number of masked bases:
9075386
```
