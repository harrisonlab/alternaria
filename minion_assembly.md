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

Building a directory structure:

```bash
ProjDir=/home/groups/harrisonlab/project_files/alternaria
cd $ProjDir

Organism=A.alternata_ssp_tenuissima
Strain=1166
OutDir=raw_dna/minion/$Organism/$Strain
mkdir -p $OutDir
RawDat=$(ls /data/scratch/nanopore_tmp_data/Alternaria/albacore_v2.1.10/Alt_albacore_v2.10_barcode02.fastq.gz)
cd $OutDir
cp -s $RawDat .
cd $ProjDir

Organism=gaisen
Strain=650
OutDir=raw_dna/minion/$Organism/$Strain
mkdir -p $OutDir
RawDat=$(ls /data/scratch/nanopore_tmp_data/Alternaria/albacore_v2.1.10/Alt_albacore_v2.10_barcode01.fastq.gz)
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
	for RawReads in $(ls raw_dna/minion/*/*/*.fastq.gz); do
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
for RawData in $(ls ../../../../../home/groups/harrisonlab/project_files/alternaria/qc_dna/minion/*/*/*q.gz); do
echo $RawData;
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc;
GenomeSz=35
# OutDir=$(dirname $RawData)
OutDir=$(echo $RawData | cut -f11,12,13,14 -d '/')
mkdir -p $OutDir
qsub $ProgDir/sub_count_nuc.sh $GenomeSz $RawData $OutDir
done
```

```bash
  for StrainDir in $(ls -d qc_dna/minion/*/*); do
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
650	39.48
1166	44.07
```
<!--
For Miseq data:
```bash
	for RawData in $(ls qc_dna/paired/*/*/*/*q.gz | grep -e 'FON_63'); do
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
 -->

### Read correction using Canu

```bash
for TrimReads in $(ls ../../../../../home/groups/harrisonlab/project_files/alternaria/qc_dna/minion/*/*/*q.gz); do
Organism=$(echo $TrimReads | rev | cut -f3 -d '/' | rev)
Strain=$(echo $TrimReads | rev | cut -f2 -d '/' | rev)
OutDir=assembly/canu-1.6/$Organism/"$Strain"
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/canu
qsub $ProgDir/sub_canu_correction.sh $TrimReads 35m $Strain $OutDir
done
```

### Assembbly using SMARTdenovo

```bash
for CorrectedReads in $(ls assembly/canu-1.6/*/*/*.trimmedReads.fasta.gz); do
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
for Assembly in $(ls assembly/SMARTdenovo/*/*/*.dmo.lay.utg); do
  Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
  Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)  
  OutDir=$(dirname $Assembly)
  qsub $ProgDir/sub_quast.sh $Assembly $OutDir
done
```


```bash
for Assembly in $(ls assembly/SMARTdenovo/*/*/*.dmo.lay.utg); do
Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
echo "$Organism - $Strain"
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/ascomycota_odb9)
OutDir=gene_pred/busco/$Organism/$Strain/assembly
qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
done
```


Error correction using racon:

```bash
for Assembly in $(ls assembly/SMARTdenovo/*/*/*.dmo.lay.utg); do
Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
echo "$Organism - $Strain"
ReadsFq=$(ls ../../../../../home/groups/harrisonlab/project_files/alternaria/qc_dna/minion/*/$Strain/*q.gz)
# ReadsFq=$(ls qc_dna/minion/*/$Strain/*q.gz)
Iterations=10
OutDir=$(dirname $Assembly)"/racon2_$Iterations"
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/racon
qsub $ProgDir/sub_racon.sh $Assembly $ReadsFq $Iterations $OutDir
done
```

```bash
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
for Assembly in $(ls assembly/SMARTdenovo/A.*/*/racon2_10/*.fasta | grep 'round_10'); do
OutDir=$(dirname $Assembly)
echo "" > tmp.txt
ProgDir=~/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
$ProgDir/remove_contaminants.py --keep_mitochondria --inp $Assembly --out $OutDir/racon_min_500bp_renamed.fasta --coord_file tmp.txt > $OutDir/log.txt
done
```

Quast and busco were run to assess the effects of racon on assembly quality:

```bash
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
for Assembly in $(ls assembly/SMARTdenovo/A.*/*/racon2_10/racon_min_500bp_renamed.fasta); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
OutDir=$(dirname $Assembly)
qsub $ProgDir/sub_quast.sh $Assembly $OutDir
done
```


```bash
# for Assembly in $(ls assembly/SMARTdenovo/F.*/*/racon*/*.fasta | grep 'FON_63' | grep 'racon_min_500bp_renamed'); do
for Assembly in $(ls assembly/SMARTdenovo/A.*/*/racon2_10/*.fasta); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/ascomycota_odb9)
OutDir=gene_pred/busco/$Organism/$Strain/assembly
# OutDir=$(dirname $Assembly)
qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
done
```
```bash
printf "Filename\tComplete\tDuplicated\tFragmented\tMissing\tTotal\n"
for File in $(ls gene_pred/busco/A*/*/assembly/*/short_summary_*.txt); do
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

Raw reads were moved onto the cluster scratch space for this step and unpacked:

```bash
for Tar in $(ls /data/scratch/nanopore_tmp_data/Alternaria/albacore_v2.1.10/Alt_albacore_v2.10_demultiplexed.tar.gz); do
  ScratchDir=/data/scratch/armita/alternaria/albacore_v2.1.10
  mkdir -p $ScratchDir
  tar -zxvf $Tar -C $ScratchDir
done
```


```bash
for Assembly in $(ls assembly/SMARTdenovo/*/1166/racon2_10/racon_min_500bp_renamed.fasta); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
# Step 1 extract reads as a .fq file which contain info on the location of the fast5 files
# Note - the full path from home must be used
ReadDir=raw_dna/nanopolish/$Organism/$Strain
mkdir -p $ReadDir
ReadsFq=$(ls ../../../../../home/groups/harrisonlab/project_files/alternaria/raw_dna/minion/*/$Strain/*.fastq.gz)
ScratchDir=/data/scratch/armita/alternaria
Fast5Dir=$ScratchDir/albacore_v2.1.10/Alt_albacore_v2.10_demultiplexed/workspace/pass/barcode01
# nanopolish index -d $Fast5Dir $ReadsFq
OutDir=$(dirname $Assembly)/nanopolish
mkdir -p $OutDir
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/nanopolish
# submit alignments for nanoppolish
qsub $ProgDir/sub_minimap2_nanopolish.sh $Assembly $ReadsFq $OutDir/nanopolish
done

for Assembly in $(ls assembly/SMARTdenovo/*/650/racon2_10/racon_min_500bp_renamed.fasta); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
# Step 1 extract reads as a .fq file which contain info on the location of the fast5 files
# Note - the full path from home must be used
ReadDir=raw_dna/nanopolish/$Organism/$Strain
mkdir -p $ReadDir
ReadsFq=$(ls ../../../../../home/groups/harrisonlab/project_files/alternaria/raw_dna/minion/*/$Strain/*.fastq.gz)
ScratchDir=/data/scratch/armita/alternaria
Fast5Dir=$ScratchDir/albacore_v2.1.10/Alt_albacore_v2.10_demultiplexed/workspace/pass/barcode02
# nanopolish index -d $Fast5Dir $ReadsFq
OutDir=$(dirname $Assembly)
mkdir -p $OutDir
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/nanopolish
# submit alignments for nanoppolish
qsub $ProgDir/sub_minimap2_nanopolish.sh $Assembly $ReadsFq $OutDir/nanopolish
done
```

 Split the assembly into 50Kb fragments an submit each to the cluster for
 nanopolish correction

```bash
for Assembly in $(ls assembly/SMARTdenovo/*/*/racon2_10/racon_min_500bp_renamed.fasta); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
OutDir=$(dirname $Assembly)/nanopolish
RawReads=$(ls ../../../../../home/groups/harrisonlab/project_files/alternaria/raw_dna/minion/*/$Strain/*.fastq.gz)
AlignedReads=$(ls $OutDir/reads.sorted.bam)

NanoPolishDir=/home/armita/prog/nanopolish/nanopolish/scripts
python $NanoPolishDir/nanopolish_makerange.py $Assembly --segment-length 50000 > $OutDir/nanopolish_range.txt

Ploidy=1
echo "nanopolish log:" > $OutDir/nanopolish_log.txt
for Region in $(cat $OutDir/nanopolish_range.txt); do
Jobs=$(qstat | grep 'sub_nanopo' | grep 'qw' | wc -l)
while [ $Jobs -gt 1 ]; do
sleep 1m
printf "."
Jobs=$(qstat | grep 'sub_nanopo' | grep 'qw' | wc -l)
done		
printf "\n"
echo $Region
echo $Region >> $OutDir/nanopolish_log.txt
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/nanopolish
qsub $ProgDir/sub_nanopolish_variants.sh $Assembly $RawReads $AlignedReads $Ploidy $Region $OutDir/$Region
done
done
```

```bash
for Assembly in $(ls assembly/SMARTdenovo/*/*/racon2_10/racon_min_500bp_renamed.fasta | grep '1166'); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
OutDir=assembly/SMARTdenovo/$Organism/$Strain/nanopolish
mkdir -p $OutDir
# cat "" > $OutDir/"$Strain"_nanoplish.fa
NanoPolishDir=/home/armita/prog/nanopolish/nanopolish/scripts
InDir=$(dirname $Assembly)
python $NanoPolishDir/nanopolish_merge.py $InDir/nanopolish/*/*.fa > $OutDir/"$Strain"_nanoplish.fa

echo "" > tmp.txt
ProgDir=~/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
$ProgDir/remove_contaminants.py --keep_mitochondria --inp $OutDir/"$Strain"_nanoplish.fa --out $OutDir/"$Strain"_nanoplish_min_500bp_renamed.fasta --coord_file tmp.txt > $OutDir/log.txt
done
```

Quast and busco were run to assess the effects of nanopolish on assembly quality:

```bash

for Assembly in $(ls assembly/SMARTdenovo/*/*/nanopolish/*_nanoplish_min_500bp_renamed.fasta); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
# Quast
OutDir=$(dirname $Assembly)
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
qsub $ProgDir/sub_quast.sh $Assembly $OutDir
# Busco
BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/ascomycota_odb9)
OutDir=gene_pred/busco/$Organism/$Strain/assembly
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
done
```


### Pilon assembly correction

Assemblies were polished using Pilon

```bash
for Assembly in $(ls assembly/SMARTdenovo/*/*/nanopolish/*_nanoplish_min_500bp_renamed.fasta); do
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
echo "$Organism - $Strain"
IlluminaDir=$(ls -d ../../../../../home/groups/harrisonlab/project_files/alternaria/qc_dna/paired/*/$Strain)
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
for Assembly in $(ls assembly/SMARTdenovo/*/*/pilon/*.fasta | grep 'pilon_10'); do
echo $Assembly
echo "" > tmp.txt
OutDir=$(dirname $Assembly)
ProgDir=~/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
$ProgDir/remove_contaminants.py --keep_mitochondria --inp $Assembly --out $OutDir/pilon_min_500bp_renamed.fasta --coord_file tmp.txt > $OutDir/log.txt
done
```

Quast and busco were run to assess the effects of pilon on assembly quality:

```bash

for Assembly in $(ls assembly/SMARTdenovo/*/*/pilon/*.fasta); do
  Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
  Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
  echo "$Organism - $Strain"
  # OutDir=$(dirname $Assembly)
  # ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
  # qsub $ProgDir/sub_quast.sh $Assembly $OutDir
  Jobs=$(qstat | grep 'sub_busco' | grep 'qw'| wc -l)
  while [ $Jobs -gt 1 ]; do
  sleep 1m
  printf "."
  Jobs=$(qstat | grep 'sub_busco' | grep 'qw'| wc -l)
  done
  printf "\n"
  ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
  BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/ascomycota_odb9)
  # OutDir=gene_pred/busco/$Organism/$Strain/assembly
  OutDir=$(dirname $Assembly)
  qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
done
```

```bash
printf "Filename\tComplete\tDuplicated\tFragmented\tMissing\tTotal\n"
for File in $(ls assembly/SMARTdenovo/*/*/pilon/*/short_summary_*.txt); do  
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
short_summary_pilon_1.txt	1301	1	4	10	1315
short_summary_pilon_2.txt	1302	1	3	10	1315
short_summary_pilon_3.txt	1302	1	3	10	1315
short_summary_pilon_4.txt	1302	1	3	10	1315
short_summary_pilon_5.txt	1302	1	3	10	1315
short_summary_pilon_6.txt	1302	1	3	10	1315
short_summary_pilon_7.txt	1302	1	3	10	1315
short_summary_pilon_8.txt	1302	1	3	10	1315
short_summary_pilon_9.txt	1302	1	3	10	1315
short_summary_pilon_10.txt	1302	1	3	10	1315
short_summary_pilon_min_500bp_renamed.txt	1302	1	3	10	131
short_summary_pilon_1.txt	1290	2	8	17	1315
short_summary_pilon_2.txt	1294	2	5	16	1315
short_summary_pilon_3.txt	1298	2	5	12	1315
short_summary_pilon_4.txt	1298	2	5	12	1315
short_summary_pilon_5.txt	1298	2	5	12	1315
short_summary_pilon_6.txt	1298	2	5	12	1315
short_summary_pilon_7.txt	1298	2	5	12	1315
short_summary_pilon_8.txt	1298	2	5	12	1315
short_summary_pilon_9.txt	1298	2	5	12	1315
short_summary_pilon_10.txt	1298	2	5	12	1315
short_summary_pilon_min_500bp_renamed.txt	1298	2	5	12	1315
```

# Hybrid Assembly

## Spades Assembly

```bash
for TrimReads in $(ls ../../../../../home/groups/harrisonlab/project_files/alternaria/raw_dna/minion/*/*/*.fastq.gz); do
Organism=$(echo $TrimReads | rev | cut -f3 -d '/' | rev)
Strain=$(echo $TrimReads | rev | cut -f2 -d '/' | rev)
IlluminaDir=$(ls -d ../../../../../home/groups/harrisonlab/project_files/alternaria/qc_dna/paired/*/$Strain)
TrimF1_Read=$(ls $IlluminaDir/F/*_trim.fq.gz | head -n1)
TrimR1_Read=$(ls $IlluminaDir/R/*_trim.fq.gz | head -n1)
OutDir=assembly/spades_minion/$Organism/"$Strain"
echo $TrimF1_Read
echo $TrimR1_Read
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/spades
qsub $ProgDir/sub_spades_minion.sh $TrimReads $TrimF1_Read $TrimR1_Read $OutDir
done
```

Contigs shorter than 500bp were removed from the assembly

```bash
  for Contigs in $(ls assembly/spades_minion/*/*/contigs.fasta); do
    AssemblyDir=$(dirname $Contigs)
    mkdir $AssemblyDir/filtered_contigs
    FilterDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/abyss
    $FilterDir/filter_abyss_contigs.py $Contigs 500 > $AssemblyDir/filtered_contigs/contigs_min_500bp.fasta
  done
```


Quast and BUSCO

```bash
for Assembly in $(ls assembly/spades_minion/*/*/filtered_contigs/contigs_min_500bp.fasta); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
OutDir=$(dirname $Assembly)
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
qsub $ProgDir/sub_quast.sh $Assembly $OutDir
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
# BuscoDB="Fungal"
BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/ascomycota_odb9)
OutDir=$(dirname $Assembly)
qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
done
```


```bash
  for File in $(ls assembly/spades_minion/*/*/*/polished/*/short_summary_*.txt ); do
  Strain=$(echo $File| rev | cut -d '/' -f5 | rev)
  Organism=$(echo $File | rev | cut -d '/' -f6 | rev)
  Prefix=$(basename $File)
  Complete=$(cat $File | grep "(C)" | cut -f2)
  Fragmented=$(cat $File | grep "(F)" | cut -f2)
  Missing=$(cat $File | grep "(M)" | cut -f2)
  Total=$(cat $File | grep "Total" | cut -f2)
  echo -e "$Organism\t$Strain\t$Prefix\t$Complete\t$Fragmented\t$Missing\t$Total"
  done
```


## Merging pacbio and hybrid assemblies


Merging was noted to improve the 1166 apple pathotype assembly but not the 650 pear pathotype assembly

```bash
for MinionAssembly in $(ls assembly/SMARTdenovo/*/*/pilon/pilon_min_500bp_renamed.fasta | grep '1166'); do
Organism=$(echo $MinionAssembly | rev | cut -f4 -d '/' | rev)
Strain=$(echo $MinionAssembly | rev | cut -f3 -d '/' | rev)
HybridAssembly=$(ls assembly/spades_*/${Organism}/${Strain}*/filtered_contigs/contigs_min_500bp.fasta)
OutDir=assembly/merged_SMARTdenovo_spades/$Organism/${Strain}
AnchorLength=500000
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/quickmerge
qsub $ProgDir/sub_quickmerge.sh $MinionAssembly $HybridAssembly $OutDir $AnchorLength
done
```

Quast and BUSCO

```bash
for Assembly in $(ls assembly/merged_SMARTdenovo_spades/*/*/merged.fasta | grep '1166'); do
Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
OutDir=$(dirname $Assembly)
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
qsub $ProgDir/sub_quast.sh $Assembly $OutDir
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/ascomycota_odb9)
OutDir=$(dirname $Assembly)
qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
done
```

```bash
  for File in $(ls assembly/merged_SMARTdenovo_spades/*/*/*/short_summary_*.txt | grep '1166'); do
  Strain=$(echo $File| rev | cut -d '/' -f3 | rev)
  Organism=$(echo $File | rev | cut -d '/' -f4 | rev)
  Prefix=$(basename $File)
  Complete=$(cat $File | grep "(C)" | cut -f2)
  Fragmented=$(cat $File | grep "(F)" | cut -f2)
  Missing=$(cat $File | grep "(M)" | cut -f2)
  Total=$(cat $File | grep "Total" | cut -f2)
  echo -e "$Organism\t$Strain\t$Prefix\t$Complete\t$Fragmented\t$Missing\t$Total"
  done
```

This merged assembly was polished using Pilon

```bash
for Assembly in $(ls assembly/merged_SMARTdenovo_spades/*/*/merged.fasta | grep '1166' | grep -v -e '_100kb' -e 'hybrid_first'); do
Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
echo "$Organism - $Strain"
IlluminaDir=$(ls -d ../../../../../home/groups/harrisonlab/project_files/alternaria/qc_dna/paired/*/$Strain)
TrimF1_Read=$(ls $IlluminaDir/F/*_trim.fq.gz | head -n1)
TrimR1_Read=$(ls $IlluminaDir/R/*_trim.fq.gz | head -n1)
OutDir=$(dirname $Assembly)/pilon
Iterations=5
# OutDir=$(dirname $Assembly)/pilon2
# Iterations=10
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/pilon
qsub $ProgDir/sub_pilon.sh $Assembly $TrimF1_Read $TrimR1_Read $OutDir $Iterations
done
```

Contigs were renamed
```bash
for Assembly in $(ls assembly/merged_SMARTdenovo_spades/*/*/pilon/*.fasta | grep 'pilon_5' | grep '1166'); do
# for Assembly in $(ls assembly/merged_SMARTdenovo_spades/*/*/pilon2/*.fasta | grep 'pilon_10'); do
echo $Assembly
echo "" > tmp.txt
OutDir=$(dirname $Assembly)
ProgDir=~/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
$ProgDir/remove_contaminants.py --keep_mitochondria --inp $Assembly --out $OutDir/pilon_min_500bp_renamed.fasta --coord_file tmp.txt > $OutDir/log.txt
done
```


BUSCO

```bash

for Assembly in $(ls assembly/merged_SMARTdenovo_spades/*/*/pilon/*.fasta | grep '1166'); do
# for Assembly in $(ls assembly/merged_SMARTdenovo_spades/*/*/pilon2/*.fasta); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
Jobs=$(qstat | grep 'sub_busco' | grep 'qw'| wc -l)
while [ $Jobs -gt 1 ]; do
sleep 1m
printf "."
Jobs=$(qstat | grep 'sub_busco' | grep 'qw'| wc -l)
done
printf "\n"
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/ascomycota_odb9)
OutDir=$(dirname $Assembly)
qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
done
```

```bash
  for File in $(ls assembly/merged_SMARTdenovo_spades/*/*/pilon/*/short_summary_*.txt | grep '1166'); do
  Strain=$(echo $File| rev | cut -d '/' -f4 | rev)
  Organism=$(echo $File | rev | cut -d '/' -f5 | rev)
  Prefix=$(basename $File)
  Complete=$(cat $File | grep "(C)" | cut -f2)
  Fragmented=$(cat $File | grep "(F)" | cut -f2)
  Missing=$(cat $File | grep "(M)" | cut -f2)
  Total=$(cat $File | grep "Total" | cut -f2)
  echo -e "$Organism\t$Strain\t$Prefix\t$Complete\t$Fragmented\t$Missing\t$Total"
  done
```


The results of merging showed worse results than when the assembly was performed
using MinION only data.



## Identifying Mitochondrial genes in assemblies

Using a blast based approach of Mt genes:

```bash
for Assembly in $(ls assembly/merged_SMARTdenovo_spades/*/1166/pilon/pilon_min_500bp_renamed.fasta assembly/SMARTdenovo/A.gaisen/650/pilon/pilon_min_500bp_renamed.fasta); do
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
echo $Assembly
Query=$(ls ../../../../home/groups/harrisonlab/project_files/alternaria/analysis/blast_homology/Mt_dna/Mt_genes_A.alternata_Liao_2017.fasta)
OutDir=analysis/Mt_genes/$Organism/$Strain
mkdir -p $OutDir
ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
qsub $ProgDir/run_blast2csv.sh $Query protein $Assembly $OutDir
done
```

Using an exclusion database with deconseq:


```bash
  for Assembly in $(ls assembly/merged_SMARTdenovo_spades/*/1166/pilon/pilon_min_500bp_renamed.fasta assembly/SMARTdenovo/A.gaisen/650/pilon/pilon_min_500bp_renamed.fasta); do
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
    echo "$Organism - $Strain"
    for Exclude_db in "Aalt_mtDNA"; do
      AssemblyDir=$(dirname $Assembly)
      OutDir=$AssemblyDir/../deconseq_$Exclude_db
      ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
      qsub $ProgDir/sub_deconseq_no_retain.sh $Assembly $Exclude_db $OutDir
    done
  done
```

Results were summarised using the commands:

```bash
for Exclude_db in "Aalt_mtDNA"; do
echo $Exclude_db
for File in $(ls assembly/*/*/*/*/log.txt | grep "$Exclude_db"); do
Name=$(echo $File | rev | cut -f3 -d '/' | rev);
Good=$(cat $File |cut -f2 | head -n1 | tail -n1);
Bad=$(cat $File |cut -f2 | head -n3 | tail -n1);
printf "$Name\t$Good\t$Bad\n";
done
done
```

```
Aalt_mtDNA
1166	22	1
650	23	1
650	27	1
```

Quast was run on the removed mtDNA:

```bash
for Assembly in $(ls assembly/*/A*/*/deconseq_Aalt_mtDNA/*_cont.fa); do
Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
OutDir=$(dirname $Assembly)
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
qsub $ProgDir/sub_quast.sh $Assembly $OutDir
done
```

# Repeat Masking

Repeat masking was performed on the hybrid assembly.

```bash
  # for Assembly in $(ls assembly/merged_SMARTdenovo_spades/*/*/pilon2/pilon_min_500bp_renamed.fasta); do
  for Assembly in $(ls assembly/merged_SMARTdenovo_spades/A.alternata_ssp_tenuissima/1166/deconseq_Aalt_mtDNA/contigs_min_500bp_filtered_renamed.fasta  assembly/SMARTdenovo/A.gaisen/650/deconseq_Aalt_mtDNA/contigs_min_500bp_filtered_renamed.fasta | grep '650'); do
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

for File in $(ls repeat_masked/*/*/filtered_contigs/*_contigs_softmasked.fa); do
OutDir=$(dirname $File)
TPSI=$(ls $OutDir/*_contigs_unmasked.fa.TPSI.allHits.chains.gff3)
OutFile=$(echo $File | sed 's/_contigs_softmasked.fa/_contigs_softmasked_repeatmasker_TPSI_appended.fa/g')
echo "$OutFile"
bedtools maskfasta -soft -fi $File -bed $TPSI -fo $OutFile
echo "Number of masked bases:"
cat $OutFile | grep -v '>' | tr -d '\n' | awk '{print $0, gsub("[a-z]", ".")}' | cut -f2 -d ' '
done
# The number of N's in hardmasked sequence are not counted as some may be present within the assembly and were therefore not repeatmasked.
for File in $(ls repeat_masked/*/*/filtered_contigs/*_contigs_hardmasked.fa); do
OutDir=$(dirname $File)
TPSI=$(ls $OutDir/*_contigs_unmasked.fa.TPSI.allHits.chains.gff3)
OutFile=$(echo $File | sed 's/_contigs_hardmasked.fa/_contigs_hardmasked_repeatmasker_TPSI_appended.fa/g')
echo "$OutFile"
bedtools maskfasta -fi $File -bed $TPSI -fo $OutFile
done
```
```
repeat_masked/A.alternata_ssp_tenuissima/1166/filtered_contigs/1166_contigs_softmasked_repeatmasker_TPSI_appended.fa
Number of masked bases:
1268326
repeat_masked/A.gaisen/650/filtered_contigs/650_contigs_softmasked_repeatmasker_TPSI_appended.fa
Number of masked bases:
749954
```

Quast and BUSCO

```bash
for Assembly in $(ls repeat_masked/*/*/filtered_contigs/*_contigs_unmasked.fa); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
OutDir=$(dirname $Assembly)
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
qsub $ProgDir/sub_quast.sh $Assembly $OutDir
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/ascomycota_odb9)
OutDir=$(dirname $Assembly)
qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
done
```


```bash
  for File in $(ls repeat_masked/*/*/filtered_contigs/run_*_contigs_unmasked/short_summary_*.txt); do
  Strain=$(echo $File| rev | cut -d '/' -f4 | rev)
  Organism=$(echo $File | rev | cut -d '/' -f5 | rev)
  Complete=$(cat $File | grep "(C)" | cut -f2)
  Fragmented=$(cat $File | grep "(F)" | cut -f2)
  Missing=$(cat $File | grep "(M)" | cut -f2)
  Total=$(cat $File | grep "Total" | cut -f2)
  echo -e "$Organism\t$Strain\t$Complete\t$Fragmented\t$Missing\t$Total"
  done
```

```
A.alternata_ssp_tenuissima	1166	1302	3	10	1315
A.gaisen	650	1298	5	12	1315
```

```bash
  for File in $(ls repeat_masked/*/*/filtered_contigs/report.tsv); do
  Strain=$(echo $File| rev | cut -d '/' -f3 | rev)
  Organism=$(echo $File | rev | cut -d '/' -f4 | rev)
  Contigs=$(cat $File | grep "contigs (>= 0 bp)" | cut -f2)
  Length=$(cat $File | grep "Total length (>= 0 bp)" | cut -f2)
  Largest=$(cat $File | grep "Largest contig" | cut -f2)
  N50=$(cat $File | grep "N50" | cut -f2)
  echo -e "$Organism\t$Strain\t$Contigs\t$Length\t$Largest\t$N50"
  done
```

```
A.alternata_ssp_tenuissima	1166	22	35704180	3902980	2583941
A.gaisen	650	27	34346950	6257968	2110033
```

# Identify Telomere repeats:

Maria wrote a script to identify presence of Telomeric repeats in assemblies

```bash
for Assembly in $(ls repeat_masked/*/*/filtered_contigs/*_contigs_unmasked.fa | grep '650'); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
OutDir=analysis/telomere/$Organism/$Strain
mkdir -p $OutDir
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/telomeres
$ProgDir/annotate_telomeres.py --fasta $Assembly > $OutDir/telomere_hits.txt
# Motif="TTAGGG"
# for Strand in '+' '-'; do
# ProgDir=/home/armita/git_repos/emr_repos/scripts/popgen/codon
# python $ProgDir/identify_telomere_repeats.py $Assembly $Motif $Strand $OutDir/${Strain}_telomere_${Motif}_${Strand}.bed
# $ProgDir/how_many_repeats_regions.py $OutDir/${Strain}_telomere_${Motif}_${Strand}.bed
# done
done
cat $OutDir/telomere_hits.txt | sort -nr -k5 | less
```

# Gene Prediction

#### Aligning


```bash
for Assembly in $(ls repeat_masked/*/*/filtered_contigs/*_contigs_unmasked.fa | grep '650'); do
Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
for FileF in $(ls ../../../../home/groups/harrisonlab/project_files/alternaria/qc_rna/paired/*/*/F/*.fastq.gz); do
Jobs=$(qstat | grep 'sub_sta' | grep 'qw'| wc -l)
while [ $Jobs -gt 1 ]; do
sleep 1m
printf "."
Jobs=$(qstat | grep 'sub_sta' | grep 'qw'| wc -l)
done
printf "\n"
FileR=$(echo $FileF | sed 's&/F/&/R/&g' | sed 's/F.fastq.gz/R.fastq.gz/g')
echo $FileF
echo $FileR
Prefix=$(echo $FileF | rev | cut -f3 -d '/' | rev)
# Timepoint=$(echo $FileF | rev | cut -f2 -d '/' | rev)
Timepoint="treatment"
#echo "$Timepoint"
OutDir=alignment/star/$Organism/$Strain/$Timepoint/$Prefix
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/RNAseq
qsub $ProgDir/sub_star.sh $Assembly $FileF $FileR $OutDir
done
done
```

Accepted hits .bam file were concatenated and indexed for use for gene model training:

```bash
  for OutDir in $(ls -d alignment/star/*/* | grep '650'); do
    Strain=$(echo $OutDir | rev | cut -d '/' -f1 | rev)
    Organism=$(echo $OutDir | rev | cut -d '/' -f2 | rev)
    echo "$Organism - $Strain"
    # For all alignments
    BamFiles=$(ls $OutDir/treatment/*/*.sortedByCoord.out.bam | tr -d '\n' | sed 's/.bam/.bam /g')
    mkdir -p $OutDir/concatenated
    samtools merge -@ 12 -f $OutDir/concatenated/concatenated.bam $BamFiles
  done
```


#### Braker prediction

Before braker predictiction was performed, I double checked that I had the
genemark key in my user area and copied it over from the genemark install
directory:

```bash
	ls ~/.gm_key
	cp /home/armita/prog/genemark/gm_key_64 ~/.gm_key
```

```bash
  for Assembly in $(ls repeat_masked/*/*/filtered_contigs/*_contigs_softmasked_repeatmasker_TPSI_appended.fa | grep '650'); do
    Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
    Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
    echo "$Organism - $Strain"
    OutDir=gene_pred/braker/$Organism/"$Strain"_braker
    AcceptedHits=$(ls alignment/star/$Organism/$Strain/concatenated/concatenated.bam)
    GeneModelName="$Organism"_"$Strain"_braker
    rm -r /home/armita/prog/augustus-3.1/config/species/"$Organism"_"$Strain"_braker
    ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/braker1
    qsub $ProgDir/sub_braker_fungi.sh $Assembly $OutDir $AcceptedHits $GeneModelName
  done
```


## Supplimenting Braker gene models with CodingQuary genes

Additional genes were added to Braker gene predictions, using CodingQuary in
pathogen mode to predict additional regions.

Firstly, aligned RNAseq data was assembled into transcripts using Cufflinks.

Note - cufflinks doesn't always predict direction of a transcript and
therefore features can not be restricted by strand when they are intersected.

```bash
for Assembly in $(ls repeat_masked/*/*/filtered_contigs/*_contigs_softmasked_repeatmasker_TPSI_appended.fa | grep '650'); do
Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
OutDir=gene_pred/cufflinks/$Organism/$Strain/concatenated
mkdir -p $OutDir
AcceptedHits=$(ls alignment/star/$Organism/$Strain/concatenated/concatenated.bam)
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/RNAseq
qsub $ProgDir/sub_cufflinks.sh $AcceptedHits $OutDir
done
```

Secondly, genes were predicted using CodingQuary:

```bash
for Assembly in $(ls repeat_masked/*/*/filtered_contigs/*_contigs_softmasked_repeatmasker_TPSI_appended.fa | grep '650'); do
Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
OutDir=gene_pred/codingquary/$Organism/$Strain
CufflinksGTF=$(ls gene_pred/cufflinks/$Organism/$Strain/concatenated/transcripts.gtf)
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/codingquary
qsub $ProgDir/sub_CodingQuary.sh $Assembly $CufflinksGTF $OutDir
done
```

The next step had problems with the masked pacbio genome. Bioperl could not read in
the fasta sequences. This was overcome by wrapping the unmasked genome and using this
fasta file.

```bash
for Assembly in $(ls repeat_masked/*/*/filtered_contigs/*_contigs_unmasked.fa); do
NewName=$(echo $Assembly | sed 's/_unmasked.fa/_unmasked_wrapped.fa/g')
cat $Assembly | fold > $NewName
done
```

Then, additional transcripts were added to Braker gene models, when CodingQuary
genes were predicted in regions of the genome, not containing Braker gene
models:

```bash
for BrakerGff in $(ls gene_pred/braker/*/*_braker/*/augustus.gff3 | grep '650'); do
Strain=$(echo $BrakerGff| rev | cut -d '/' -f3 | rev | sed 's/_braker_new//g' | sed 's/_braker_pacbio//g' | sed 's/_braker//g')
Organism=$(echo $BrakerGff | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
Assembly=$(ls repeat_masked/$Organism/$Strain/filtered_contigs/*_contigs_softmasked_repeatmasker_TPSI_appended.fa)
# Assembly=$(ls repeat_masked/$Organism/$Strain/filtered_contigs/*_contigs_unmasked_wrapped.fa)
CodingQuaryGff=$(ls gene_pred/codingquary/$Organism/$Strain/out/PredictedPass.gff3)
PGNGff=$(ls gene_pred/codingquary/$Organism/$Strain/out/PGN_predictedPass.gff3)
AddDir=gene_pred/codingquary/$Organism/$Strain/additional
FinalDir=gene_pred/final/$Organism/$Strain/final
AddGenesList=$AddDir/additional_genes.txt
AddGenesGff=$AddDir/additional_genes.gff
FinalGff=$AddDir/combined_genes.gff
mkdir -p $AddDir
mkdir -p $FinalDir

for x in $CodingQuaryGff $PGNGff; do
  bedtools intersect -v -a $x -b $BrakerGff | grep 'gene'| cut -f2 -d'=' | cut -f1 -d';'
done > $AddGenesList

for y in $CodingQuaryGff $PGNGff; do
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation
  $ProgDir/gene_list_to_gff.pl $AddGenesList $y CodingQuarry_v2.0 ID CodingQuary
done > $AddGenesGff
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/codingquary
# -
# This section is edited
$ProgDir/add_CodingQuary_features.pl $AddGenesGff $Assembly > $AddDir/add_genes_CodingQuary_unspliced.gff3
$ProgDir/correct_CodingQuary_splicing.py --inp_gff $AddDir/add_genes_CodingQuary_unspliced.gff3 > $FinalDir/final_genes_CodingQuary.gff3
# -
$ProgDir/gff2fasta.pl $Assembly $FinalDir/final_genes_CodingQuary.gff3 $FinalDir/final_genes_CodingQuary
cp $BrakerGff $FinalDir/final_genes_Braker.gff3
$ProgDir/gff2fasta.pl $Assembly $FinalDir/final_genes_Braker.gff3 $FinalDir/final_genes_Braker
cat $FinalDir/final_genes_Braker.pep.fasta $FinalDir/final_genes_CodingQuary.pep.fasta | sed -r 's/\*/X/g' > $FinalDir/final_genes_combined.pep.fasta
cat $FinalDir/final_genes_Braker.cdna.fasta $FinalDir/final_genes_CodingQuary.cdna.fasta > $FinalDir/final_genes_combined.cdna.fasta
cat $FinalDir/final_genes_Braker.gene.fasta $FinalDir/final_genes_CodingQuary.gene.fasta > $FinalDir/final_genes_combined.gene.fasta
cat $FinalDir/final_genes_Braker.upstream3000.fasta $FinalDir/final_genes_CodingQuary.upstream3000.fasta > $FinalDir/final_genes_combined.upstream3000.fasta

GffBraker=$FinalDir/final_genes_CodingQuary.gff3
GffQuary=$FinalDir/final_genes_Braker.gff3
GffAppended=$FinalDir/final_genes_appended.gff3
cat $GffBraker $GffQuary > $GffAppended
done
```

The final number of genes per isolate was observed using:
```bash
  for DirPath in $(ls -d gene_pred/final/*/*/final); do
    Strain=$(echo $DirPath| rev | cut -d '/' -f2 | rev)
    Organism=$(echo $DirPath | rev | cut -d '/' -f3 | rev)
    echo "$Organism - $Strain"
    cat $DirPath/final_genes_Braker.pep.fasta | grep '>' | wc -l;
    cat $DirPath/final_genes_CodingQuary.pep.fasta | grep '>' | wc -l;
    cat $DirPath/final_genes_combined.pep.fasta | grep '>' | wc -l;
    echo "";
  done
```

The number of genes predicted by Braker, supplimented by CodingQuary and in the
final combined dataset was shown:

```
A.alternata_ssp_tenuissima - 1166
12761
904
13665

A.gaisen - 650
12023
851
12874
```


In preperation for submission to ncbi, gene models were renamed and duplicate gene features were identified and removed.
 * no duplicate genes were identified



```bash
for GffAppended in $(ls gene_pred/final/*/*/final/final_genes_appended.gff3); do
Strain=$(echo $GffAppended | rev | cut -d '/' -f3 | rev)
Organism=$(echo $GffAppended | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
FinalDir=gene_pred/final/$Organism/$Strain/final
GffFiltered=$FinalDir/filtered_duplicates.gff
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/codingquary
$ProgDir/remove_dup_features.py --inp_gff $GffAppended --out_gff $GffFiltered
GffRenamed=$FinalDir/final_genes_appended_renamed.gff3
LogFile=$FinalDir/final_genes_appended_renamed.log
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/codingquary
$ProgDir/gff_rename_genes.py --inp_gff $GffFiltered --conversion_log $LogFile > $GffRenamed
rm $GffFiltered
Assembly=$(ls repeat_masked/$Organism/$Strain/filtered_contigs/*_softmasked_repeatmasker_TPSI_appended.fa)
$ProgDir/gff2fasta.pl $Assembly $GffRenamed gene_pred/final/$Organism/$Strain/final/final_genes_appended_renamed
# The proteins fasta file contains * instead of Xs for stop codons, these should
# be changed
sed -i 's/\*/X/g' gene_pred/final/$Organism/$Strain/final/final_genes_appended_renamed.pep.fasta
done
```

```
A.alternata_ssp_tenuissima - 1166
Identifiied the following duplicated transcripts:
CUFF_8299_1_996.t2
CUFF_10096_1_77.t2
A.gaisen - 650
Identifiied the following duplicated transcripts:
CUFF_2980_1_77.t2
CUFF_1001_1_1582.t2
CUFF_2980_1_80.t2
CUFF_3334_1_60.t2
CUFF_720_1_86.t2
CUFF_7702_1_36.t2
CUFF_7380_1_535.t2
CUFF_9418_1_13.t2
```

```bash
for File in $(ls gene_pred/final/*/*/final/final_genes_appended_renamed.pep.fasta); do
  echo $File | rev | cut -f3 -d '/' |  rev
  cat $File | grep '>' | wc -l
done
```

```
1166
13663
650
12866
```

## Checking gene prediction accruacy using BUSCO

```bash
 for Assembly in $(ls gene_pred/final/*/*/final/final_genes_appended_renamed.gene.fasta); do
   Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
   Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
   echo "$Organism - $Strain"
   ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
   BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/ascomycota_odb9)
   OutDir=gene_pred/busco/$Organism/$Strain/transcript
   qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
 done
```

```bash
 for File in $(ls gene_pred/busco/*/*/transcript/*/short_summary_*.txt); do
 Strain=$(echo $File| rev | cut -d '/' -f4 | rev)
 Organism=$(echo $File | rev | cut -d '/' -f5 | rev)
 Complete=$(cat $File | grep "(C)" | cut -f2)
 Single=$(cat $File | grep "(S)" | cut -f2)
 Fragmented=$(cat $File | grep "(F)" | cut -f2)
 Missing=$(cat $File | grep "(M)" | cut -f2)
 Total=$(cat $File | grep "Total" | cut -f2)
 echo -e "$Organism\t$Strain\t$Complete\t$Single\t$Fragmented\t$Missing\t$Total"
 done
```


#Functional annotation

## A) Interproscan

Interproscan was used to give gene models functional annotations.
Annotation was run using the commands below:

Note: This is a long-running script. As such, these commands were run using
'screen' to allow jobs to be submitted and monitored in the background.
This allows the session to be disconnected and reconnected over time.

Screen ouput detailing the progress of submission of interporscan jobs
was redirected to a temporary output file named interproscan_submission.log .

```bash
	screen -a
	ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/interproscan
	for Genes in $(ls gene_pred/final/*/*/final/final_genes_appended_renamed.pep.fasta | grep -e '1166' -e '650'); do
	echo $Genes
	$ProgDir/sub_interproscan.sh $Genes
	done 2>&1 | tee -a interproscan_submisison.log
```

Following interproscan annotation split files were combined using the following
commands:

```bash
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/interproscan
for Proteins in $(ls gene_pred/final/*/*/final/final_genes_appended_renamed.pep.fasta | grep -e '1166' -e '650'); do
Strain=$(echo $Proteins | rev | cut -d '/' -f3 | rev)
Organism=$(echo $Proteins | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
echo $Strain
InterProRaw=gene_pred/interproscan/$Organism/$Strain/raw
$ProgDir/append_interpro.sh $Proteins $InterProRaw
done
```


## B) SwissProt


```bash
for Proteome in $(ls gene_pred/final/*/*/final/final_genes_appended_renamed.pep.fasta | grep -e '1166' -e '650'); do
Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
OutDir=gene_pred/swissprot/$Organism/$Strain
SwissDbDir=../../../../home/groups/harrisonlab/uniprot/swissprot
SwissDbName=uniprot_sprot
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/swissprot
qsub $ProgDir/sub_swissprot.sh $Proteome $OutDir $SwissDbDir $SwissDbName
done
```

## Small secreted proteins

Putative effectors identified within Augustus gene models using a number
of approaches:

 * A) From Braker gene models - Signal peptide & small cystein rich protein
 <!-- * B) From Augustus gene models - Hmm evidence of WY domains  
 * C) From Augustus gene models - Hmm evidence of RxLR effectors
 * D) From Augustus gene models - Hmm evidence of CRN effectors  
 * E) From ORF fragments - Signal peptide & RxLR motif  
 * F) From ORF fragments - Hmm evidence of WY domains  
 * G) From ORF fragments - Hmm evidence of RxLR effectors -->


### A) From Augustus gene models - Identifying secreted proteins

Required programs:
 * SigP
 * biopython
 * TMHMM


Proteins that were predicted to contain signal peptides were identified using
the following commands:

```bash
for Proteome in $(ls gene_pred/final/*/*/final/final_genes_appended_renamed.pep.fasta | grep -e '1166' -e '650'); do
SplitfileDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/signal_peptides
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/signal_peptides
Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
SplitDir=gene_pred/braker_split/$Organism/$Strain
mkdir -p $SplitDir
BaseName="$Organism""_$Strain"_braker_preds
$SplitfileDir/splitfile_500.py --inp_fasta $Proteome --out_dir $SplitDir --out_base $BaseName
for File in $(ls $SplitDir/*_braker_preds_*); do
Jobs=$(qstat | grep 'pred_sigP' | wc -l)
while [ $Jobs -gt '20' ]; do
sleep 10
printf "."
Jobs=$(qstat | grep 'pred_sigP' | wc -l)
done
printf "\n"
echo $File
qsub $ProgDir/pred_sigP.sh $File signalp-4.1
done
done
```

The batch files of predicted secreted proteins needed to be combined into a
single file for each strain. This was done with the following commands:
```bash
  for SplitDir in $(ls -d gene_pred/braker_split/*/*); do
    Strain=$(echo $SplitDir | cut -d '/' -f4)
    Organism=$(echo $SplitDir | cut -d '/' -f3)
    InStringAA=''
    InStringNeg=''
    InStringTab=''
    InStringTxt=''
    SigpDir=braker_signalp-4.1
    echo "$Organism - $Strain"
    for GRP in $(ls -l $SplitDir/*_braker_preds_*.fa | rev | cut -d '_' -f1 | rev | sort -n); do  
      InStringAA="$InStringAA gene_pred/$SigpDir/$Organism/$Strain/split/"$Organism"_"$Strain"_braker_preds_$GRP""_sp.aa";  
      InStringNeg="$InStringNeg gene_pred/$SigpDir/$Organism/$Strain/split/"$Organism"_"$Strain"_braker_preds_$GRP""_sp_neg.aa";  
      InStringTab="$InStringTab gene_pred/$SigpDir/$Organism/$Strain/split/"$Organism"_"$Strain"_braker_preds_$GRP""_sp.tab";
      InStringTxt="$InStringTxt gene_pred/$SigpDir/$Organism/$Strain/split/"$Organism"_"$Strain"_braker_preds_$GRP""_sp.txt";  
    done
    cat $InStringAA > gene_pred/$SigpDir/$Organism/$Strain/"$Strain"_aug_sp.aa
    cat $InStringNeg > gene_pred/$SigpDir/$Organism/$Strain/"$Strain"_aug_neg_sp.aa
    tail -n +2 -q $InStringTab > gene_pred/$SigpDir/$Organism/$Strain/"$Strain"_aug_sp.tab
    cat $InStringTxt > gene_pred/$SigpDir/$Organism/$Strain/"$Strain"_aug_sp.txt
  done
```

Some proteins that are incorporated into the cell membrane require secretion.
Therefore proteins with a transmembrane domain are not likely to represent
cytoplasmic or apoplastic effectors.

Proteins containing a transmembrane domain were identified:

```bash
  for Proteome in $(ls gene_pred/final/*/*/final/final_genes_appended_renamed.pep.fasta | grep -e '1166' -e '650'); do
    Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/transmembrane_helices
    qsub $ProgDir/submit_TMHMM.sh $Proteome
  done
```

Those proteins with transmembrane domains were removed from lists of Signal
peptide containing proteins

```bash
for File in $(ls gene_pred/trans_mem/*/*/*_TM_genes_neg.txt); do
Strain=$(echo $File | rev | cut -f2 -d '/' | rev)
Organism=$(echo $File | rev | cut -f3 -d '/' | rev)
# echo "$Organism - $Strain"
NonTmHeaders=$(echo "$File" | sed 's/neg.txt/neg_headers.txt/g')
cat $File | cut -f1 > $NonTmHeaders
SigP=$(ls gene_pred/braker_signalp-4.1/$Organism/$Strain/"$Strain"_aug_sp.aa)
OutDir=$(dirname $SigP)
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/extract_from_fasta.py --fasta $SigP --headers $NonTmHeaders > $OutDir/"$Strain"_final_sp_no_trans_mem.aa
# echo "Number of SigP proteins:"
TotalProts=$(cat $SigP | grep '>' | wc -l)
# echo "Number without transmembrane domains:"
SecProt=$(cat $OutDir/"$Strain"_final_sp_no_trans_mem.aa | grep '>' | wc -l)
# echo "Number of gene models:"
SecGene=$(cat $OutDir/"$Strain"_final_sp_no_trans_mem.aa | grep '>' | cut -f1 -d't' | sort | uniq |wc -l)
# A text file was also made containing headers of proteins testing +ve
PosFile=$(ls gene_pred/trans_mem/$Organism/$Strain/"$Strain"_TM_genes_pos.txt)
TmHeaders=$(echo $PosFile | sed 's/.txt/_headers.txt/g')
cat $PosFile | cut -f1 > $TmHeaders
printf "$Organism\t$Strain\t$TotalProts\t$SecProt\t$SecGene\n"
done
```

```
A.alternata_ssp_tenuissima	1166	1511	1253	1251
A.gaisen	650	1429	1186	1184
```

### C) From Augustus gene models - Effector identification using EffectorP

Required programs:
 * EffectorP.py

```bash
  for Proteome in $(ls gene_pred/final/*/*/final/final_genes_appended_renamed.pep.fasta | grep -e '1166' -e '650'); do
    Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
    BaseName="$Organism"_"$Strain"_EffectorP
    OutDir=analysis/effectorP/$Organism/$Strain
    ProgDir=~/git_repos/emr_repos/tools/seq_tools/feature_annotation/fungal_effectors
    qsub $ProgDir/pred_effectorP.sh $Proteome $BaseName $OutDir
  done
```

Those genes that were predicted as secreted and tested positive by effectorP
were identified:

Note - this doesnt exclude proteins with TM domains or GPI anchors

```bash
  for File in $(ls analysis/effectorP/*/*/*_EffectorP.txt); do
    Strain=$(echo $File | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $File | rev | cut -f3 -d '/' | rev)
    echo "$Organism - $Strain"
    Headers=$(echo "$File" | sed 's/_EffectorP.txt/_EffectorP_headers.txt/g')
    cat $File | grep 'Effector' | grep -v 'Effector probability:' | cut -f1 > $Headers
    printf "EffectorP headers:\t"
    cat $Headers | wc -l
    Secretome=$(ls gene_pred/braker_signalp-4.1/$Organism/$Strain/"$Strain"_final_sp_no_trans_mem.aa)
    OutFile=$(echo "$File" | sed 's/_EffectorP.txt/_EffectorP_secreted.aa/g')
    ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
    $ProgDir/extract_from_fasta.py --fasta $Secretome --headers $Headers > $OutFile
    OutFileHeaders=$(echo "$File" | sed 's/_EffectorP.txt/_EffectorP_secreted_headers.txt/g')
    cat $OutFile | grep '>' | tr -d '>' > $OutFileHeaders
    printf "Secreted effectorP headers:\t"
    cat $OutFileHeaders | wc -l
    Gff=$(ls gene_pred/final/$Organism/$Strain/*/final_genes_appended_renamed.gff3)
    EffectorP_Gff=$(echo "$File" | sed 's/_EffectorP.txt/_EffectorP_secreted.gff/g')
    ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
    $ProgDir/extract_gff_for_sigP_hits.pl $OutFileHeaders $Gff effectorP ID > $EffectorP_Gff
  done
```

```
A.alternata_ssp_tenuissima - 1166
EffectorP headers:	1970
Secreted effectorP headers:	248
A.gaisen - 650
EffectorP headers:	1921
Secreted effectorP headers:	247
```

## SSCP

Small secreted cysteine rich proteins were identified within secretomes. These
proteins may be identified by EffectorP, but this approach allows direct control
over what constitutes a SSCP.

```bash
for Secretome in $(ls gene_pred/braker_signalp-4.1/*/*/*_final_sp_no_trans_mem.aa); do
Strain=$(echo $Secretome| rev | cut -f2 -d '/' | rev)
Organism=$(echo $Secretome | rev | cut -f3 -d '/' | rev)
echo "$Organism - $Strain"
OutDir=analysis/sscp/$Organism/$Strain
mkdir -p $OutDir
ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/sscp
$ProgDir/sscp_filter.py --inp_fasta $Secretome --max_length 300 --threshold 3 --out_fasta $OutDir/"$Strain"_sscp_all_results.fa
cat $OutDir/"$Strain"_sscp_all_results.fa | grep 'Yes' > $OutDir/"$Strain"_sscp.fa
printf "number of SSC-rich genes:\t"
cat $OutDir/"$Strain"_sscp.fa | grep '>' | tr -d '>' | cut -f1 -d '.' | sort | uniq | wc -l
done
```

```
A.alternata_ssp_tenuissima - 1166
% cysteine content threshold set to:	3
maximum length set to:	300
No. short-cysteine rich proteins in input fasta:	202
number of SSC-rich genes:	202
A.gaisen - 650
% cysteine content threshold set to:	3
maximum length set to:	300
No. short-cysteine rich proteins in input fasta:	197
number of SSC-rich genes:	196
```

## CAZY proteins

Carbohydrte active enzymes were idnetified using CAZYfollowing recomendations
at http://csbl.bmb.uga.edu/dbCAN/download/readme.txt :

```bash
for Proteome in $(ls gene_pred/final/*/*/final/final_genes_appended_renamed.pep.fasta | grep -e '1166' -e '650'); do
Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
OutDir=gene_pred/CAZY/$Organism/$Strain
mkdir -p $OutDir
Prefix="$Strain"_CAZY
CazyHmm=../../../../home/groups/harrisonlab/dbCAN/dbCAN-fam-HMMs.txt
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/HMMER
qsub $ProgDir/sub_hmmscan.sh $CazyHmm $Proteome $Prefix $OutDir
done
```

The Hmm parser was used to filter hits by an E-value of E1x10-5 or E1x10-e3 if they had a hit over a length of X %.

Those proteins with a signal peptide were extracted from the list and gff files
representing these proteins made.

```bash
for File in $(ls gene_pred/CAZY/*/*/*CAZY.out.dm); do
Strain=$(echo $File | rev | cut -f2 -d '/' | rev)
Organism=$(echo $File | rev | cut -f3 -d '/' | rev)
OutDir=$(dirname $File)
# echo "$Organism - $Strain"
ProgDir=/home/groups/harrisonlab/dbCAN
$ProgDir/hmmscan-parser.sh $OutDir/"$Strain"_CAZY.out.dm > $OutDir/"$Strain"_CAZY.out.dm.ps
CazyHeaders=$(echo $File | sed 's/.out.dm/_headers.txt/g')
cat $OutDir/"$Strain"_CAZY.out.dm.ps | cut -f3 | sort | uniq > $CazyHeaders
# echo "number of CAZY proteins identified:"
TotalProts=$(cat $CazyHeaders | wc -l)
# Gff=$(ls gene_pred/codingquary/$Organism/$Strain/final/final_genes_appended_renamed.gff3)
Gff=$(ls gene_pred/final/$Organism/$Strain/final/final_genes_appended_renamed.gff3)
CazyGff=$OutDir/"$Strain"_CAZY.gff
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/extract_gff_for_sigP_hits.pl $CazyHeaders $Gff CAZyme ID > $CazyGff
# echo "number of CAZY genes identified:"
TotalGenes=$(cat $CazyGff | grep -w 'gene' | wc -l)

SecretedProts=$(ls gene_pred/braker_signalp-4.1/$Organism/$Strain/"$Strain"_final_sp_no_trans_mem.aa)
SecretedHeaders=$(echo $SecretedProts | sed 's/.aa/_headers.txt/g')
cat $SecretedProts | grep '>' | tr -d '>' > $SecretedHeaders
CazyGffSecreted=$OutDir/"$Strain"_CAZY_secreted.gff
$ProgDir/extract_gff_for_sigP_hits.pl $SecretedHeaders $CazyGff Secreted_CAZyme ID > $CazyGffSecreted
# echo "number of Secreted CAZY proteins identified:"
cat $CazyGffSecreted | grep -w 'mRNA' | cut -f9 | tr -d 'ID=' | cut -f1 -d ';' > $OutDir/"$Strain"_CAZY_secreted_headers.txt
SecProts=$(cat $OutDir/"$Strain"_CAZY_secreted_headers.txt | wc -l)
# echo "number of Secreted CAZY genes identified:"
SecGenes=$(cat $CazyGffSecreted | grep -w 'gene' | wc -l)
# cat $OutDir/"$Strain"_CAZY_secreted_headers.txt | cut -f1 -d '.' | sort | uniq | wc -l
printf "$Organism\t$Strain\t$TotalProts\t$TotalGenes\t$SecProts\t$SecGenes\n"
done
```

```
A.alternata_ssp_tenuissima	1166	777	777	389	389
A.gaisen	650	764	764	374	374
```

Note - the CAZY genes identified may need further filtering based on e value and
cuttoff length - see below:

Cols in yourfile.out.dm.ps:
1. Family HMM
2. HMM length
3. Query ID
4. Query length
5. E-value (how similar to the family HMM)
6. HMM start
7. HMM end
8. Query start
9. Query end
10. Coverage

* For fungi, use E-value < 1e-17 and coverage > 0.45

* The best threshold varies for different CAZyme classes (please see http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4132414/ for details). Basically to annotate GH proteins, one should use a very relax coverage cutoff or the sensitivity will be low (Supplementary Tables S4 and S9); (ii) to annotate CE families a very stringent E-value cutoff and coverage cutoff should be used; otherwise the precision will be very low due to a very high false positive rate (Supplementary Tables S5 and S10)


### Summary of CAZY families by organism


```bash
for CAZY in $(ls gene_pred/CAZY/*/*/*_CAZY.out.dm.ps); do
Strain=$(echo $CAZY | rev | cut -f2 -d '/' | rev)
Organism=$(echo $CAZY | rev | cut -f3 -d '/' | rev)
OutDir=$(dirname $CAZY)
echo "$Organism - $Strain"
Secreted=$(ls gene_pred/braker_signalp-4.1/$Organism/$Strain/*_final_sp_no_trans_mem_headers.txt)
Gff=$(ls gene_pred/final/$Organism/$Strain/final/final_genes_appended_renamed.gff3)
ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/CAZY
$ProgDir/summarise_CAZY.py --cazy $CAZY --inp_secreted $Secreted --inp_gff $Gff --summarise_family --trim_gene_id 2 --kubicek_2014
done
```

```
A.alternata_ssp_tenuissima - 1166
B-Galactosidases - 4
A-Galactosidases - 3
Polygalacturonase - 8
A-Arabinosidases - 14
Xylanases - 14
Polygalacturonate lyases - 21
B-Glucuronidases - 2
B-Glycosidases - 12
Cellulases - 30
Xyloglucanases - 1
other - 279


A.gaisen - 650
B-Galactosidases - 3
A-Galactosidases - 2
Polygalacturonase - 11
A-Arabinosidases - 14
Xylanases - 14
Polygalacturonate lyases - 20
B-Glucuronidases - 4
B-Glycosidases - 11
Cellulases - 29
Xyloglucanases - 1
other - 265
```

## D) Secondary metabolites (Antismash and SMURF)

Antismash was run to identify clusters of secondary metabolite genes within
the genome. Antismash was run using the weserver at:
http://antismash.secondarymetabolites.org


Results of web-annotation of gene clusters within the assembly were downloaded to
the following directories:

```bash
  for Assembly in $(ls repeat_masked/*/*/filtered_contigs/*_contigs_softmasked_repeatmasker_TPSI_appended.fa); do
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    OutDir=gene_pred/secondary_metabolites/antismash/$Organism/$Strain
    mkdir -p $OutDir
  done
```

```bash
  for Zip in $(ls gene_pred/secondary_metabolites/antismash/*/*/*.zip); do
    OutDir=$(dirname $Zip)
    unzip -d $OutDir $Zip
  done
```

```bash
for AntiSmash in $(ls gene_pred/secondary_metabolites/antismash/*/*/*/*.final.gbk); do
    Organism=$(echo $AntiSmash | rev | cut -f4 -d '/' | rev)
    Strain=$(echo $AntiSmash | rev | cut -f3 -d '/' | rev)
    echo "$Organism - $Strain"
    OutDir=gene_pred/secondary_metabolites/antismash/$Organism/$Strain
    Prefix=$OutDir/WT_antismash
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/secondary_metabolites
    $ProgDir/antismash2gff.py --inp_antismash $AntiSmash --out_prefix $Prefix

    # Identify secondary metabolites within predicted clusters
    printf "Number of secondary metabolite detected:\t"
    cat "$Prefix"_secmet_clusters.gff | wc -l
    GeneGff=gene_pred/final/$Organism/$Strain/final/final_genes_appended_renamed.gff3
    bedtools intersect -u -a $GeneGff -b "$Prefix"_secmet_clusters.gff > "$Prefix"_secmet_genes.gff
    cat "$Prefix"_secmet_genes.gff | grep -w 'mRNA' | cut -f9 | cut -f2 -d '=' | cut -f1 -d ';' > "$Prefix"_antismash_secmet_genes.txt
    bedtools intersect -wo -a $GeneGff -b "$Prefix"_secmet_clusters.gff | grep 'mRNA' | cut -f9,10,12,18 | sed "s/ID=//g" | perl -p -i -e "s/;Parent=g\w+//g" | perl -p -i -e "s/;Notes=.*//g" > "$Prefix"_secmet_genes.tsv
    printf "Number of predicted proteins in secondary metabolite clusters:\t"
    cat "$Prefix"_secmet_genes.tsv | wc -l
    printf "Number of predicted genes in secondary metabolite clusters:\t"
    cat "$Prefix"_secmet_genes.gff | grep -w 'gene' | wc -l

      # Identify cluster finder additional non-secondary metabolite clusters
      printf "Number of cluster finder non-SecMet clusters detected:\t"
      cat "$Prefix"_clusterfinder_clusters.gff | wc -l
      GeneGff=gene_pred/final/$Organism/$Strain/final/final_genes_appended_renamed.gff3
      bedtools intersect -u -a $GeneGff -b "$Prefix"_clusterfinder_clusters.gff > "$Prefix"_clusterfinder_genes.gff
      cat "$Prefix"_clusterfinder_genes.gff | grep -w 'mRNA' | cut -f9 | cut -f2 -d '=' | cut -f1 -d ';' > "$Prefix"_clusterfinder_genes.txt

      printf "Number of predicted proteins in cluster finder non-SecMet clusters:\t"
      cat "$Prefix"_clusterfinder_genes.txt | wc -l
      printf "Number of predicted genes in cluster finder non-SecMet clusters:\t"
      cat "$Prefix"_clusterfinder_genes.gff | grep -w 'gene' | wc -l
  done
```


These clusters represented the following genes. Note that these numbers just
show the number of intersected genes with gff clusters and are not confirmed by
function

```
A.alternata_ssp_tenuissima - 1166
Number of secondary metabolite detected:	34
Number of predicted proteins in secondary metabolite clusters:	946
Number of predicted genes in secondary metabolite clusters:	888
Number of cluster finder non-SecMet clusters detected:	102
Number of predicted proteins in cluster finder non-SecMet clusters:	2280
Number of predicted genes in cluster finder non-SecMet clusters:	2266
A.gaisen - 650
Number of secondary metabolite detected:	29
Number of predicted proteins in secondary metabolite clusters:	765
Number of predicted genes in secondary metabolite clusters:	763
Number of cluster finder non-SecMet clusters detected:	96
Number of predicted proteins in cluster finder non-SecMet clusters:	2145
Number of predicted genes in cluster finder non-SecMet clusters:	2136
```

SMURF was also run to identify secondary metabolite gene clusters.

Genes needed to be parsed into a specific tsv format prior to submission on the
SMURF webserver.

```bash
  for Gff in $(ls gene_pred/final/*/*/final/final_genes_appended_renamed.gff3); do
    Strain=$(echo $Gff | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Gff | rev | cut -f4 -d '/' | rev)
    OutDir=analysis/secondary_metabolites/smurf/$Organism/$Strain
    mkdir -p $OutDir
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/secondary_metabolites
    $ProgDir/gff2smurf.py --gff $Gff > $OutDir/"$Strain"_genes_smurf.tsv
  done
```

SMURF output was received by email and downloaded to the cluster in the output
directory above.

Output files were parsed into gff format:

```bash
  for OutDir in $(ls -d analysis/secondary_metabolites/smurf/*/*); do
    Strain=$(echo $OutDir | rev | cut -f1 -d '/' | rev)
    Organism=$(echo $OutDir | rev | cut -f2 -d '/' | rev)
    GeneGff=$(ls gene_pred/final/$Organism/$Strain/final/final_genes_appended_renamed.gff3)
    SmurfClusters=$(ls $OutDir/Secondary-Metabolite-Clusters.txt)
    SmurfBackbone=$(ls $OutDir/Backbone-genes.txt)
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/secondary_metabolites
    $ProgDir/smurf2gff.py --smurf_clusters $SmurfClusters --smurf_backbone $SmurfBackbone > $OutDir/Smurf_clusters.gff
   bedtools intersect -wo -a $GeneGff -b $OutDir/Smurf_clusters.gff | grep 'mRNA' | cut -f9,10,12,18 | sed "s/ID=//g" | perl -p -i -e "s/;Parent=g\w+//g" | perl -p -i -e "s/;Notes=.*//g" > $OutDir/"$Strain"_smurf_secmet_genes.tsv
  done
```

# Genes with transcription factor annotations:


A list of PFAM domains, superfamily annotations used as part of the DBD database
and a further set of interproscan annotations listed by Shelest et al 2017 were made
http://www.transcriptionfactor.org/index.cgi?Domain+domain:all
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5415576/

```bash
for Interpro in $(ls gene_pred/interproscan/*/*/*_interproscan.tsv); do
Organism=$(echo $Interpro | rev | cut -f3 -d '/' | rev)
Strain=$(echo $Interpro | rev | cut -f2 -d '/' | rev)
# echo "$Organism - $Strain"
OutDir=analysis/transcription_factors/$Organism/$Strain
mkdir -p $OutDir
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/transcription_factors
$ProgDir/interpro2TFs.py --InterPro $Interpro > $OutDir/"$Strain"_TF_domains.tsv
# echo "total number of transcription factors"
cat $OutDir/"$Strain"_TF_domains.tsv | cut -f1 | sort | uniq > $OutDir/"$Strain"_TF_gene_headers.txt
NumTF=$(cat $OutDir/"$Strain"_TF_gene_headers.txt | wc -l)
printf "$Organism\t$Strain\t$NumTF\n"
done
```

```
A.alternata_ssp_tenuissima	1166	690
A.gaisen	650	633
```


# Genomic analysis
The first analysis was based upon BLAST searches for genes known to be involved in toxin production


## Genes with homology to PHIbase
Predicted gene models were searched against the PHIbase database using tBLASTx.

```bash
qlogin -pe smp 12
cd /data/scratch/armita/alternaria
dbFasta=$(ls /home/groups/harrisonlab/phibase/v4.4/phi_accessions.fa)
dbType="prot"
for QueryFasta in $(ls gene_pred/final/*/*/final/final_genes_appended_renamed.cds.fasta); do
Organism=$(echo $QueryFasta | rev | cut -f4 -d '/' | rev)
Strain=$(echo $QueryFasta | rev | cut -f3 -d '/' | rev)
echo "$Organism - $Strain"
Prefix="${Strain}_phi_accessions"
Eval="1e-30"
OutDir=analysis/blast_homology/$Organism/$Strain
mkdir -p $OutDir
makeblastdb -in $dbFasta -input_type fasta -dbtype $dbType -title $Prefix.db -parse_seqids -out $OutDir/$Prefix.db
blastx -num_threads 12 -db $OutDir/$Prefix.db -query $QueryFasta -outfmt 6 -num_alignments 1 -out $OutDir/${Prefix}_hits.txt -evalue $Eval
cat $OutDir/${Prefix}_hits.txt | grep 'effector' | cut -f1,2 | sort | uniq > $OutDir/${Prefix}_hits_headers.txt
done
```

## Presence and genes with homology to Alternaria toxins

The first analysis was based upon BLAST searches for genes known to be involved in toxin production


```bash
qlogin -pe smp 4
cd /data/scratch/armita/alternaria
dbFasta=$(ls /home/groups/harrisonlab/project_files/alternaria/analysis/blast_homology/CDC_genes/A.alternata_CDC_genes.fa)
dbType="nucl"
for QueryFasta in $(ls gene_pred/final/*/*/final/final_genes_appended_renamed.cds.fasta); do
Organism=$(echo $QueryFasta | rev | cut -f4 -d '/' | rev)
Strain=$(echo $QueryFasta | rev | cut -f3 -d '/' | rev)
echo "$Organism - $Strain"
Prefix="${Strain}_CDC_genes"
Eval="1e-100"
OutDir=analysis/blast_homology/$Organism/$Strain
mkdir -p $OutDir
makeblastdb -in $dbFasta -input_type fasta -dbtype $dbType -title $Prefix.db -parse_seqids -out $OutDir/$Prefix.db
tblastx -num_threads 4 -db $OutDir/$Prefix.db -query $QueryFasta -outfmt 6 -num_alignments 1 -out $OutDir/${Prefix}_hits.txt -evalue $Eval
cat $OutDir/${Prefix}_hits.txt | cut -f1,2 | sort | uniq > $OutDir/${Prefix}_hits_headers.txt
# cat $OutDir/${Prefix}_hits.txt | grep 'effector' | cut -f1,2 | sort | uniq > $OutDir/${Prefix}_hits_headers.txt
done
```


Blast searches were also performed in comparison to the genome sequence

```bash
  for Subject in $(ls repeat_masked/*/*/filtered_contigs/*_contigs_softmasked_repeatmasker_TPSI_appended.fa); do
    ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
    Query=../../../../home/groups/harrisonlab/project_files/alternaria/analysis/blast_homology/CDC_genes/A.alternata_CDC_genes.fa
    qsub $ProgDir/blast_pipe.sh $Query dna $Subject
  done
```

```bash
  HitsList=""
  StrainList=""
  for BlastHits in $(ls analysis/blast_homology/*/*/*_A.alternata_CDC_genes.fa_homologs.csv); do
    sed -i "s/\t\t/\t/g" $BlastHits
    sed -i "s/Grp1\t//g" $BlastHits
    Strain=$(echo $BlastHits | rev | cut -f2 -d '/' | rev)
    HitsList="$HitsList $BlastHits"
    StrainList="$StrainList $Strain"
    ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
    $ProgDir/blast_parse.py --blast_csv $HitsList --headers $StrainList --identity 0.7 --evalue 1e-30
  done
```

```
Query ID's	1166	650
AB015351_-_AKT1_gene	0	2
AB015352_-_AKT2_gene	0	1
AB034586_-_ACTT1_gene	0	1
AB035491_-_AKTR_gene	0	2
AB035492_-_AKT3_gene	0	2
AB070711_-_AFT1-1_gene	0	2
AB070712_-_AFTR-1_gene	0	2
AB070713_-_AFT3-1_gene	0	2
AB119280_-_AFTS1_gene	1	0
AB176941_ACTT3_gene	0	2
AB176941_ACTTR_gene	0	2
AB179766_-_AFT3-2_gene	0	2
AB179766_-_AFT9-1_gene	0	2
AB179766_-_AFT10-1_gene	0	2
AB179766_-_AFT11-1_gene	0	2
AB179766_-_AFT12-1_gene	0	2
AB179766_-_AFTR-2_gene	0	2
AB432914_ACTT2_gene	0	1
AB444613_ACTT5_gene	0	1
AB444614_ACTT6_gene	0	2
AB465676_-_ALT1_gene	0	0
AB525198_-_AMT1_gene	2	0
AB525198_-_AMT2_gene	3	1
AB525198_-_AMT3_gene	2	0
AB525198_-_AMT4_gene	2	0
AB525198_-_AMT5_gene	2	0
AB525198_-_AMT6_gene	2	0
AB525198_-_AMT7_gene	2	0
AB525198_-_AMT8_gene	2	0
AB525198_-_AMT9_gene	2	0
AB525198_-_AMT10_gene	1	0
AB525198_-_AMT11_gene	0	0
AB525198_-_AMT12_gene	2	0
AB525198_-_AMT13_gene	1	0
AB525198_-_AMT14_gene	1	1
AB525198_-_AMT15_gene	0	0
AB525198_-_AMT16_gene	3	2
AB525198_-_AMTR1_gene	2	2
AB688098_-_ACRTS1_gene	0	0
AB725683_-_ACRTS2_gene	0	0
```

Once blast searches had completed, the BLAST hits were converted to GFF
annotations:

```bash
  for BlastHits in $(ls analysis/blast_homology/*/*/*_A.alternata_CDC_genes.fa_homologs.csv | grep '650'); do
    ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
    Strain=$(echo $BlastHits | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $BlastHits | rev | cut -f3 -d '/' | rev)
    HitsGff=analysis/blast_homology/$Organism/$Strain/"$Strain"_A.alternata_CDC_genes.fa_homologs.gff
    Column2=toxin_homolog
    NumHits=1
    $ProgDir/blast2gff.pl $Column2 $NumHits $BlastHits > $HitsGff
    cat $HitsGff | grep 'AKT' > analysis/blast_homology/$Organism/$Strain/"$Strain"_A.alternata_CDC_genes.fa_homologs_AKT.gff
    cat $HitsGff | grep 'AMT' > analysis/blast_homology/$Organism/$Strain/"$Strain"_A.alternata_CDC_genes.fa_homologs_AMT.gff
  done
```


Extracted gff files of the BLAST hit locations were intersected with gene models
to identify if genes were predicted for these homologs:

```bash
for HitsGff in $(ls analysis/blast_homology/*/*/*_A.alternata_CDC_genes.fa_homologs.gff); do
Strain=$(echo $HitsGff | rev | cut -f2 -d '/' | rev)
Organism=$(echo $HitsGff | rev | cut -f3 -d '/' | rev)
Proteins=$(ls gene_pred/final/$Organism/$Strain/*/final_genes_appended_renamed.gff3)
OutDir=analysis/blast_homology/$Organism/$Strain/"$Strain"_A.alternata_CDC_genes
IntersectBlast=$OutDir/"$Strain"_A.alternata_CDC_genes_Intersect.gff
# NoIntersectBlast=$OutDir/"$Strain"_A.alternata_CDC_genes_NoIntersect.gff
mkdir -p $OutDir
echo "$Organism - $Strain"
echo "The number of BLAST hits in the gff file were:"
cat $HitsGff | wc -l
echo "The number of blast hits intersected were:"
bedtools intersect -wao -a $HitsGff -b $Proteins > $IntersectBlast
cat $IntersectBlast | grep -w 'gene' | wc -l
echo "The number of blast hits not intersecting gene models were:"
# bedtools intersect -v -a $HitsGff -b $Proteins > $NoIntersectBlast
# rm $NoIntersectBlast
cat $IntersectBlast | grep -w -E '0$' | wc -l
done
```

```
A.alternata_ssp_tenuissima - 1166
The number of BLAST hits in the gff file were:
34
The number of blast hits intersected were:
32
The number of blast hits not intersecting gene models were:
8
A.gaisen - 650
The number of BLAST hits in the gff file were:
37
The number of blast hits intersected were:
33
The number of blast hits not intersecting gene models were:
10
```



# Build Annotation Tables


```bash
for GeneGff in $(ls gene_pred/final/*/*/final/final_genes_appended_renamed.gff3 | grep -e '1166' -e '650'); do
  Strain=$(echo $GeneGff | rev | cut -f3 -d '/' | rev)
  Organism=$(echo $GeneGff | rev | cut -f4 -d '/' | rev)
  Assembly=$(ls repeat_masked/$Organism/$Strain/filtered_contigs/*_contigs_unmasked.fa)
  Antismash=$(ls gene_pred/secondary_metabolites/antismash/$Organism/$Strain/*_antismash_secmet_genes.tsv)
  # Smurf=$(ls analysis/secondary_metabolites/smurf/$Organism/$Strain/*_smurf_secmet_genes.tsv)
  TFs=$(ls analysis/transcription_factors/$Organism/$Strain/"$Strain"_TF_domains.tsv)
  SigP=$(ls gene_pred/braker_signalp-4.1/$Organism/$Strain/"$Strain"_aug_sp.aa)
  TM_out=$(ls gene_pred/trans_mem/$Organism/$Strain/"$Strain"_TM_genes_pos.txt)
  # GPI_out=$(ls gene_pred/trans_mem/$Organism/$Strain/GPIsom/GPI_pos.fa)
  EffP_list=$(ls analysis/effectorP/$Organism/$Strain/"$Organism"_"$Strain"_EffectorP_headers.txt)
  CAZY_list=$(ls gene_pred/CAZY/$Organism/$Strain/"$Strain"_CAZY.out.dm.ps)
  PhiHits=$(ls analysis/blast_homology/$Organism/$Strain/"$Strain"_phi_accessions_hits_headers.txt)
  ToxinHits=$(ls analysis/blast_homology/$Organism/$Strain/"$Strain"_CDC_genes_hits_headers.txt)
  InterPro=$(ls gene_pred/interproscan/$Organism/$Strain/*_interproscan.tsv)
  SwissProt=$(ls gene_pred/swissprot/$Organism/$Strain/swissprot_vMar2018_tophit_parsed.tbl)
  # Orthology=$(ls analysis/orthology/orthomcl/At_Aa_Ag_all_isolates/formatted/Results_Apr10/Orthogroups.txt)
  Orthology=$(ls analysis/orthology/orthomcl/At_Aa_Ag_all_isolates/formatted/Results_May04/Orthogroups.txt)
  if [[ $Strain == '648' ]]; then
    OrthoStrainID='At_1'
    echo $OrthoStrainID
  elif [[ $Strain == '1082' ]]; then
    OrthoStrainID='At_2'
    echo $OrthoStrainID
  elif [[ $Strain == '1164' ]]; then
    OrthoStrainID='At_3'
    echo $OrthoStrainID
  elif [[ $Strain == '24350' ]]; then
    OrthoStrainID='At_4'
    echo $OrthoStrainID
  elif [[ $Strain == '635' ]]; then
    OrthoStrainID='At_5'
    echo $OrthoStrainID
  elif [[ $Strain == '743' ]]; then
    OrthoStrainID='At_6'
    echo $OrthoStrainID
  elif [[ $Strain == '1166' ]]; then
    OrthoStrainID='At_7'
    echo $OrthoStrainID
  elif [[ $Strain == '1177' ]]; then
    OrthoStrainID='At_8'
    echo $OrthoStrainID
  elif [[ $Strain == '675' ]]; then
    OrthoStrainID='Aa_1'
    echo $OrthoStrainID
  elif [[ $Strain == '97.0013' ]]; then
    OrthoStrainID='Aa_2'
    echo $OrthoStrainID
  elif [[ $Strain == '97.0016' ]]; then
    OrthoStrainID='Aa_3'
    echo $OrthoStrainID
  elif [[ $Strain == '650' ]]; then
    OrthoStrainID='Ag_1'
    echo $OrthoStrainID
  fi
  OrthoStrainAll='At_1 At_2 At_3 At_4 At_5 At_6 At_6 At_7 At_8 Aa_1 Aa_2 Aa_3 Ag_1'
  OutDir=gene_pred/annotation/$Organism/$Strain
  mkdir -p $OutDir
  ProgDir=/home/armita/git_repos/emr_repos/scripts/alternaria/pathogen/annotation
  # $ProgDir/build_annot_table_Alt2.py --genes_gff $GeneGff --SigP $SigP --TM_list $TM_out --EffP_list $EffP_list --CAZY_list $CAZY_list --TFs $TFs --PhiHits $PhiHits --ToxinHits $ToxinHits --Antismash $Antismash --Smurf $Smurf --InterPro $InterPro --Swissprot $SwissProt --orthogroups $Orthology --strain_id $OrthoStrainID --OrthoMCL_all $OrthoStrainAll > $OutDir/"$Strain"_annotation_ncbi.tsv
  # SMURF results were exluceded due to them being considered untrustworthy
  $ProgDir/build_annot_table_Alt2.py --genes_gff $GeneGff --SigP $SigP --TM_list $TM_out --EffP_list $EffP_list --CAZY_list $CAZY_list --TFs $TFs --PhiHits $PhiHits --ToxinHits $ToxinHits --Antismash $Antismash --InterPro $InterPro --Swissprot $SwissProt --orthogroups $Orthology --strain_id $OrthoStrainID --OrthoMCL_all $OrthoStrainAll > $OutDir/"$Strain"_annotation_ncbi.tsv
done
```


LS regions of the apple and pear pathotype were investigated:


```bash
# 1166
Isolate='1166'
AnnotTab=$(ls gene_pred/annotation/A.*/$Isolate/${Isolate}_annotation_ncbi.tsv)
GenesCDC=$(echo $AnnotTab | sed 's/_annotation_ncbi.tsv/_CDC_genes.tsv/g')
cat $AnnotTab | grep -e 'contig_14' -e 'contig_15' -e 'contig_18' -e 'contig_19' -e 'contig_20' -e 'contig_21' > $GenesCDC
TotalGenes=$(cat $GenesCDC | wc -l)
Secreted=$(cat $GenesCDC | cut -f1,8 | grep 'Yes' | wc -l)
EffectorP=$(cat $GenesCDC | cut -f1,8,9 | grep "Yes.Yes" | wc -l)
CAZyme=$(cat $GenesCDC | cut -f1,8,10 | grep 'CAZY' | grep 'Yes' | wc -l)
SecMet=$(cat $GenesCDC | grep 'AS_' | wc -l)
printf "$Isolate\t$TotalGenes\t$Secreted\t$EffectorP\t$CAZyme\t$SecMet\n"

# 650
Isolate='650'
AnnotTab=$(ls gene_pred/annotation/A.*/$Isolate/${Isolate}_annotation_ncbi.tsv)
GenesCDC=$(echo $AnnotTab | sed 's/_annotation_ncbi.tsv/_CDC_genes.tsv/g')
cat $AnnotTab | grep -e 'contig_12' -e 'contig_14' -e 'contig_19' -e 'contig_20' -e 'contig_21' -e 'contig_22' -e 'contig_23' > $GenesCDC
TotalGenes=$(cat $GenesCDC | wc -l)
Secreted=$(cat $GenesCDC | cut -f1,8 | grep 'Yes' | wc -l)
EffectorP=$(cat $GenesCDC | cut -f1,8,9 | grep "Yes.Yes" | wc -l)
CAZyme=$(cat $GenesCDC | cut -f1,8,10 | grep 'CAZY' | grep 'Yes' | wc -l)
SecMet=$(cat $GenesCDC | grep 'AS_' | wc -l)
printf "$Isolate\t$TotalGenes\t$Secreted\t$EffectorP\t$CAZyme\t$SecMet\n"
```

```
1166	624	32	12	6	153
650	504	38	12	7	140
```

# Genomic analysis


## Alignment vs minion genomes

Illumina sequence data from each of the isoaltes was aligned against the minion
assemblies.


### Read alignemnt

```bash
for Reference in $(ls repeat_masked/*/*/filtered_contigs/*_contigs_unmasked.fa); do
RefStrain=$(echo $Reference | rev | cut -f3 -d '/' | rev)
for StrainPath in $(ls -d ../../../../home/groups/harrisonlab/project_files/alternaria/qc_dna/paired/*/* | grep 'arborescens'); do
Jobs=$(qstat | grep 'sub_bowtie' | grep 'qw'| wc -l)
while [ $Jobs -gt 1 ]; do
sleep 1m
printf "."
Jobs=$(qstat | grep 'sub_bowtie' | grep 'qw'| wc -l)
done
printf "\n"
echo $StrainPath
Strain=$(echo $StrainPath | rev | cut -f1 -d '/' | rev)
Organism=$(echo $StrainPath | rev | cut -f2 -d '/' | rev)
echo "$Organism - $Strain"
F_Read=$(ls $StrainPath/F/*_trim.fq.gz | head -n1 | tail -n1)
R_Read=$(ls $StrainPath/R/*_trim.fq.gz | head -n1 | tail -n1)
echo $F_Read
echo $R_Read
OutDir=analysis/genome_alignment/bowtie/$Organism/$Strain/vs_${RefStrain}
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/genome_alignment
qsub $ProgDir/bowtie/sub_bowtie.sh $Reference $F_Read $R_Read $OutDir $Strain
done
done
```

### Read coverage

Identify read coverage over each bp

```bash
  for Bam in $(ls analysis/genome_alignment/bowtie/*/*/vs_*/*_aligned_sorted.bam); do
    Target=$(echo $Bam | rev | cut -f2 -d '/' | rev)
    Strain=$(echo $Bam | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Bam | rev | cut -f4 -d '/' | rev)
    echo "$Organism - $Strain - $Target"
    OutDir=$(dirname $Bam)
    samtools depth -aa $Bam > $OutDir/${Organism}_${Strain}_${Target}_depth.tsv
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/genome_alignment/coverage_analysis
    $ProgDir/cov_by_window.py --cov $OutDir/${Organism}_${Strain}_${Target}_depth.tsv > $OutDir/${Organism}_${Strain}_${Target}_depth_10kb.tsv
    sed -i "s/$/\t$Strain/g" $OutDir/${Organism}_${Strain}_${Target}_depth_10kb.tsv
  done
  for Target in "vs_650" "vs_1166"; do
    OutDir=analysis/genome_alignment/bowtie/grouped_${Target}
    mkdir -p $OutDir
    cat analysis/genome_alignment/bowtie/*/*/*/*_*_${Target}_depth_10kb.tsv > $OutDir/${Target}_grouped_depth.tsv
  done
```

### Plot read coverage


```R
library(readr)
setwd("~/Downloads/Aalt/coverage2")

appended_df <- read_delim("~/Downloads/Aalt/coverage2/vs_1166_grouped_depth.tsv", "\t", escape_double = FALSE, col_names = FALSE, col_types = cols(X4 = col_factor(levels = c("675", "97.0013", "97.0016", "650", "648", "24350", "1082", "1164", "635", "743", "1166", "1177"))), trim_ws = TRUE)


colnames(appended_df) <- c("contig","position", "depth", "strain")

appended_df$treatment <- paste(appended_df$strain , appended_df$contig)
tapply(appended_df$depth, appended_df$treatment, myFun)

df2 <- cbind(do.call(rbind, tapply(appended_df$depth, appended_df$treatment, myFun)))
write.csv(df2, '1166_contig_coverage.csv')

appended_df$depth <- ifelse(appended_df$depth > 100, 100, appended_df$depth)

# install.packages("ggplot2")
library(ggplot2)
require(scales)

for (i in 1:22){
contig = paste("contig", i, sep = "_")
p0 <- ggplot(data=appended_df[appended_df$contig == contig, ], aes(x=`position`, y=`depth`, group=1)) +
geom_line() +
labs(x = "Position (bp)", y = "Coverage") +
scale_y_continuous(breaks=seq(0,100,25), limits=c(0,100)) +
facet_wrap(~strain, nrow = 12, ncol = 1, strip.position = "left")
outfile = paste("1166_contig", i, "cov.jpg", sep = "_")
ggsave(outfile , plot = p0, device = 'jpg', path = NULL,
scale = 1, width = 500, height = 500, units = 'mm',
dpi = 150, limitsize = TRUE)
}

myFun <- function(x) {
  c(min = min(x), max = max(x),
    mean = mean(x), median = median(x),
    std = sd(x))
}


```


```R
library(readr)
setwd("~/Downloads/Aalt/coverage2")

df_650 <- read_delim("~/Downloads/Aalt/coverage2/vs_650_grouped_depth.tsv", "\t", escape_double = FALSE, col_names = FALSE, col_types = cols(X4 = col_factor(levels = c("675", "97.0013", "97.0016", "650", "648", "24350", "1082", "1164", "635", "743", "1166", "1177"))), trim_ws = TRUE)

colnames(df_650) <- c("contig","position", "depth", "strain")

df_650$treatment <- paste(df_650$strain , df_650$contig)
tapply(df_650$depth, df_650$treatment, myFun)

df2 <- cbind(do.call(rbind, tapply(df_650$depth, df_650$treatment, myFun)))
write.csv(df2, '650_contig_coverage.csv')

df_650$depth <- ifelse(df_650$depth > 100, 100, df_650$depth)

# install.packages("ggplot2")
library(ggplot2)
require(scales)

for (i in 1:23){
contig = paste("contig", i, sep = "_")
p0 <- ggplot(data=df_650[df_650$contig == contig, ], aes(x=`position`, y=`depth`, group=1)) +
geom_line() +
labs(x = "Position (bp)", y = "Coverage") +
scale_y_continuous(breaks=seq(0,100,25), limits=c(0,100)) +
facet_wrap(~strain, nrow = 12, ncol = 1, strip.position = "left")
outfile = paste("650_contig", i, "cov.jpg", sep = "_")
ggsave(outfile , plot = p0, device = 'jpg', path = NULL,
scale = 1, width = 500, height = 500, units = 'mm',
dpi = 150, limitsize = TRUE)
}
```

## 4.0 Characterising CDCs

### 4.1.b Repetative content

The % repeatmasking was identified for each contig:

```bash
for Assembly in $(ls repeat_masked/*/*/filtered_contigs/*_contigs_softmasked_repeatmasker_TPSI_appended.fa); do
  Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
  Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
  OutDir=analysis/CDC_contigs/$Organism/$Strain
  mkdir -p $OutDir
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/repeat_masking
  $ProgDir/count_Ns.py --inp_fasta $Assembly --out_txt $OutDir/N_content.txt
done
```

### 4.1.c Gene density and % content

```bash
for Assembly in $(ls repeat_masked/*/*/filtered_contigs/*_contigs_softmasked_repeatmasker_TPSI_appended.fa); do
  Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
  Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
  OutDir=analysis/CDC_contigs/$Organism/$Strain
  mkdir -p $OutDir
  Genes=$(ls gene_pred/final/$Organism/$Strain/final/final_genes_appended_renamed.gff3)
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation
  $ProgDir/feature_density.py --inp_fasta $Assembly --inp_gff $Genes --features gene
done
```
```
# 1166
contig, contig_lgth, feat_count, feat_density, feat_lgth, perc_feat
contig_1	3902980	1563	400.46	2412918	61.82
contig_2	3781932	1413	373.62	2203815	58.27
contig_3	3000390	1154	384.62	1733294	57.77
contig_4	2851745	1098	385.03	1685701	59.11
contig_5	2693844	1055	391.63	1663091	61.74
contig_6	2583941	986	381.59	1540588	59.62
contig_7	2502671	954	381.19	1498182	59.86
contig_8	2455819	960	390.91	1481084	60.31
contig_9	2451092	930	379.42	1359530	55.47
contig_10	2402550	872	362.95	1436796	59.8
contig_11	2194253	818	372.79	1268004	57.79
contig_12	1303685	504	386.6	753202	57.77
contig_13	767820	282	367.27	408936	53.26
contig_14	549494	184	334.85	237154	43.16
contig_15	435297	166	381.35	162962	37.44
contig_16	433285	165	380.81	236067	54.48
contig_17	391795	169	431.35	227914	58.17
contig_18	368761	109	295.58	154186	41.81
contig_19	225742	69	305.66	81190	35.97
contig_20	154051	35	227.2	57127	37.08
contig_21	143777	48	333.85	62379	43.39
contig_22	109256	42	384.42	63072	57.73
# 650
contig_1	6257968	2443	390.38	3763556	60.14
contig_2	4971742	1914	384.98	2952852	59.39
contig_3	2925786	1160	396.47	1682886	57.52
contig_4	2776593	1084	390.41	1684148	60.66
contig_5	2407617	910	377.97	1333497	55.39
contig_6	2321443	890	383.38	1354817	58.36
contig_7	2110033	790	374.4	1226331	58.12
contig_8	1975041	733	371.13	1184196	59.96
contig_9	1696734	660	388.98	937666	55.26
contig_10	1416541	541	381.92	845372	59.68
contig_11	1110254	434	390.9	650661	58.6
contig_12	629970	201	319.06	242443	38.48
contig_13	621992	220	353.7	329368	52.95
contig_14	547262	203	370.94	267573	48.89
contig_15	463030	168	362.83	253364	54.72
contig_16	353559	132	373.35	177450	50.19
contig_17	350965	129	367.56	183181	52.19
contig_18	288921	104	359.96	143279	49.59
contig_19	206182	63	305.56	84853	41.15
contig_20	89891	26	289.24	36102	40.16
contig_21	25203	3	119.03	3216	12.76
contig_23	19593	1	51.04	389	1.99
```

### 4.1.c
The number of putative effectors was identified in Core, lineage specific and species specific contigs:
