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
#   	Species=F.oxysporum_fsp_fragariae
#   	Strain=FeChina
#     mkdir -p raw_dna/paired/fusarium_ex_strawberry/FeChina/F
#     mkdir -p raw_dna/paired/fusarium_ex_strawberry/FeChina/R
#     RawDat=/home/groups/harrisonlab/raw_data/raw_seq/raw_reads/150925_M01678_0029_AC669
#     cp $RawDat/FeChina_S1_L001_R1_001.fastq.gz raw_dna/paired/fusarium_ex_strawberry/FeChina/F/.
#     cp $RawDat/FeChina_S1_L001_R1_001.fastq.gz raw_dna/paired/fusarium_ex_strawberry/FeChina/R/.
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
	for StrainPath in $(ls -d raw_dna/paired/*/*); do
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/rna_qc
		IlluminaAdapters=/home/armita/git_repos/emr_repos/tools/seq_tools/ncbi_adapters.fa
		ReadsF=$(ls $StrainPath/F/*.fastq*)
		ReadsR=$(ls $StrainPath/R/*.fastq*)
		echo $ReadsF
		echo $ReadsR
		qsub $ProgDir/rna_qc_fastq-mcf.sh $ReadsF $ReadsR $IlluminaAdapters DNA
	done
```
Data quality was visualised once again following trimming:

```bash
  for RawData in $(ls qc_dna/paired/*/*/*/*.fq.gz); do
    echo $RawData;
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc;
    qsub $ProgDir/run_fastqc.sh $RawData;
  done
```



Sequencing coverage was estimated:

```bash
for RawData in $(ls qc_dna/paired/*/*/*/*fq.gz | grep -e 'appended' -e 's_6' | grep -e 's_6'); do
echo $RawData;
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc;
qsub $ProgDir/run_fastqc.sh $RawData;
# GenomeSz=65
GenomeSz=33
OutDir=$(dirname $RawData)
qsub $ProgDir/sub_count_nuc.sh $GenomeSz $RawData $OutDir
done
```

Find predicted coverage for these isolates:

```bash
for StrainDir in $(ls -d qc_dna/paired/*/*); do
Strain=$(basename $StrainDir)
printf "$Strain\t"
for File in $(ls qc_dna/paired/*/"$Strain"/*/*.txt); do
echo $(basename $File);
cat $File | tail -n1 | rev | cut -f2 -d ' ' | rev;
done | grep -v '.txt' | awk '{ SUM += $1} END { print SUM }'
done
```

```
675	39.11
97.0013	40.1
97.0016	30.29
650	39.24
1082	24.6
1164	32.82
1166	ls: cannot access qc_dna/paired/*/1166/*/*.txt: No such file or directory

1177	76.19
24350	38.17
635	33.77
648	57.59
743	50.08
```

# Assembly
Assembly was performed using: Velvet / Abyss / Spades

## Spades Assembly



```bash
  for StrainPath in $(ls -d qc_dna/paired/*/*); do
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/spades
    Strain=$(echo $StrainPath | rev | cut -f1 -d '/' | rev)
    Organism=$(echo $StrainPath | rev | cut -f2 -d '/' | rev)
    F_Read=$(ls $StrainPath/F/*.fq.gz)
    R_Read=$(ls $StrainPath/R/*.fq.gz)
    OutDir=assembly/spades/$Organism/$Strain
    echo $F_Read
    echo $R_Read
    qsub $ProgDir/submit_SPAdes.sh $F_Read $R_Read $OutDir correct 10
  done
```

Quast

```bash
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
    for Assembly in $(ls assembly/spades/*/*/ncbi_edits/contigs_min_500bp.fasta); do
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
    OutDir=assembly/spades/$Organism/$Strain/ncbi_edits
    qsub $ProgDir/sub_quast.sh $Assembly $OutDir
  done
```

The results of quast were shown using the following commands:

```bash
  for Assembly in $(ls assembly/spades/*/*/ncbi_edits/report.txt); do
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev);
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev);
    echo;
    echo $Organism;
    echo $Strain;
    cat $Assembly;
  done > assembly/quast_results.txt
```

The output of this analysis is in the assembly/quast_results.txt file of this
git repository.


Contigs were renamed in accordance with ncbi recomendations.

```bash
  ProgDir=~/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
  touch tmp.csv
  for Assembly in $(ls assembly/spades/*/*/ncbi_edits/contigs_min_500bp.fasta); do
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
    OutDir=assembly/spades/$Organism/$Strain/ncbi_edits
    $ProgDir/remove_contaminants.py --inp $Assembly --out $OutDir/contigs_min_500bp_renamed.fasta --coord_file tmp.csv
  done
  rm tmp.csv
```


A Bioproject and Biosample was made with NCBI genbank for submission of genomes.
Following the creation of these submissions, the .fasta assembly was uploaded
through the submission portal. A note was provided requesting that the assembly
be run through the contamination screen to aid a more detailed resubmission in
future. The returned FCSreport.txt was downloaded from the NCBI webportal and
used to correct the assembly to NCBI standards.

NCBI reports (FCSreport.txt) were manually downloaded to the following loactions:

```bash
  for Assembly in $(ls assembly/spades/*/*/*/contigs_min_500bp_renamed.fasta | grep -v 'ncbi_edits'); do
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
    NCBI_report_dir=genome_submission/$Organism/$Strain/initial_submission
    mkdir -p $NCBI_report_dir
  done
```


These downloaded files were used to correct assemblies:

```bash
for Assembly in $(ls assembly/spades/*/*/ncbi_edits/contigs_min_500bp_renamed.fasta); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
NCBI_report=$(ls genome_submission/$Organism/$Strain/initial_submission/Contamination*.txt)
OutDir=assembly/spades/$Organism/$Strain/ncbi_edits
mkdir -p $OutDir
ProgDir=~/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
$ProgDir/remove_contaminants.py --keep_mitochondria --inp $Assembly --out $OutDir/contigs_min_500bp_renamed.fasta --coord_file $NCBI_report > $OutDir/log.txt
done
```

Quast

```bash
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
  for Assembly in $(ls assembly/spades/*/*/ncbi_edits/contigs_min_500bp_renamed.fasta); do
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
    OutDir=assembly/spades/$Organism/$Strain/ncbi_edits
    qsub $ProgDir/sub_quast.sh $Assembly $OutDir
  done
```

The results of quast were shown using the following commands:

```bash
  for Assembly in $(ls assembly/spades/*/*/ncbi_edits/report.txt); do
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev);
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev);
    echo;
    echo $Organism;
    echo $Strain;
    cat $Assembly;
  done > assembly/quast_results.txt
```


# Repeat masking
Repeat masking was performed and used the following programs: Repeatmasker Repeatmodeler

The renamed assembly was used to perform repeatmasking.

```bash
for Assembly in $(ls assembly/spades/*/*/ncbi_edits/contigs_min_500bp_renamed.fasta); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
echo "$Organism"
echo "$Strain"
OutDir=repeat_masked/$Organism/$Strain/ncbi_edits_repmask
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/repeat_masking
qsub $ProgDir/rep_modeling.sh $Assembly $OutDir
qsub $ProgDir/transposonPSI.sh $Assembly $OutDir
done
```

The TransposonPSI masked bases were used to mask additional bases from the
repeatmasker / repeatmodeller softmasked and hardmasked files.

```bash
for File in $(ls repeat_masked/*/*/ncbi_edits_repmask/*_contigs_softmasked.fa | grep -v -e '1166' -e '650'); do
OutDir=$(dirname $File)
TPSI=$(ls $OutDir/*_contigs_unmasked.fa.TPSI.allHits.chains.gff3)
OutFile=$(echo $File | sed 's/_contigs_softmasked.fa/_contigs_softmasked_repeatmasker_TPSI_appended.fa/g')
echo "$OutFile"
bedtools maskfasta -soft -fi $File -bed $TPSI -fo $OutFile
echo "Number of masked bases:"
cat $OutFile | grep -v '>' | tr -d '\n' | awk '{print $0, gsub("[a-z]", ".")}' | cut -f2 -d ' '
done
# The number of N's in hardmasked sequence are not counted as some may be present within the assembly and were therefore not repeatmasked.
for File in $(ls repeat_masked/*/*/ncbi_edits_repmask/*_contigs_softmasked.fa | grep -v -e '1166' -e '650'); do
OutDir=$(dirname $File)
TPSI=$(ls $OutDir/*_contigs_unmasked.fa.TPSI.allHits.chains.gff3)
OutFile=$(echo $File | sed 's/_contigs_hardmasked.fa/_contigs_hardmasked_repeatmasker_TPSI_appended.fa/g')
echo "$OutFile"
bedtools maskfasta -fi $File -bed $TPSI -fo $OutFile
done
```

```bash
for RepDir in $(ls -d repeat_masked/*/*/ncbi_edits_repmask); do
Strain=$(echo $RepDir | rev | cut -f2 -d '/' | rev)
Organism=$(echo $RepDir | rev | cut -f3 -d '/' | rev)  
RepMaskGff=$(ls $RepDir/*_contigs_hardmasked.gff)
TransPSIGff=$(ls $RepDir/*_contigs_unmasked.fa.TPSI.allHits.chains.gff3)
# printf "The number of bases masked by RepeatMasker:\t"
RepMaskerBp=$(sortBed -i $RepMaskGff | bedtools merge | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}')
# printf "The number of bases masked by TransposonPSI:\t"
TpsiBp=$(sortBed -i $TransPSIGff | bedtools merge | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}')
# printf "The total number of masked bases are:\t"
Total=$(cat $RepMaskGff $TransPSIGff | sortBed | bedtools merge | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}')
printf "$Organism\t$Strain\t$RepMaskerBp\t$TpsiBp\t$Total\n"
done
```

```
A.alternata_ssp._arborescens	675	764901	330735	853273
A.alternata_ssp._arborescens	97.0013	846080	346346	955052
A.alternata_ssp._arborescens	97.0016	653247	348895	783649
A.alternata_ssp._gaisen	650	347610	211034	480031
A.alternata_ssp._tenuissima	1082	485588	245543	627818
A.alternata_ssp._tenuissima	1164	845729	268881	983226
A.alternata_ssp._tenuissima	1166	563743	251126	684063
A.alternata_ssp._tenuissima	1177	699137	269647	865317
A.alternata_ssp._tenuissima	24350	406564	160521	466485
A.alternata_ssp._tenuissima	635	747976	334677	953149
A.alternata_ssp._tenuissima	648	538913	159935	657864
A.alternata_ssp._tenuissima	743	703937	329935	868641
```

The breakdown of repeats were shown using the following commands

```bash
printf "Organism\tStrain\tDDE_1\tGypsy\tHAT\tTY1_Copia\tMariner\tCacta\tLINE\tMuDR_A_B\tHelitronORF\tMariner_ant1\tISC1316\tCrypton\n"
for File in $(ls /home/groups/harrisonlab/project_files/alternaria/repeat_masked/A.*/*/ncbi_edits_repmask/*_contigs_unmasked.fa.TPSI.allHits.chains.chains.gff3); do
Strain=$(echo $File| rev | cut -d '/' -f3 | rev)
Organism=$(echo $File | rev | cut -d '/' -f4 | rev)
# echo "$Organism - $Strain"
DDE_1=$(cat $File | cut -f9 | cut -f2 -d';' | sort| uniq -c | grep 'DDE_1' | sed "s/^\s*//g" | cut -f1 -d ' ')
Gypsy=$(cat $File | cut -f9 | cut -f2 -d';' | sort| uniq -c | grep 'gypsy' | sed "s/^\s*//g" | cut -f1 -d ' ')
HAT=$(cat $File | cut -f9 | cut -f2 -d';' | sort| uniq -c | grep 'hAT' | sed "s/^\s*//g" | cut -f1 -d ' ')
TY1_Copia=$(cat $File | cut -f9 | cut -f2 -d';' | sort| uniq -c | grep 'TY1_Copia' | sed "s/^\s*//g" | cut -f1 -d ' ')
Mariner=$(cat $File | cut -f9 | cut -f2 -d';' | sort| uniq -c | grep -w 'mariner' | sed "s/^\s*//g" | cut -f1 -d ' ')
Cacta=$(cat $File | cut -f9 | cut -f2 -d';' | sort| uniq -c | grep 'cacta' | sed "s/^\s*//g" | cut -f1 -d ' ')
LINE=$(cat $File | cut -f9 | cut -f2 -d';' | sort| uniq -c | grep 'LINE' | sed "s/^\s*//g" | cut -f1 -d ' ')
MuDR_A_B=$(cat $File | cut -f9 | cut -f2 -d';' | sort| uniq -c | grep 'MuDR_A_B' | sed "s/^\s*//g" | cut -f1 -d ' ')
HelitronORF=$(cat $File | cut -f9 | cut -f2 -d';' | sort| uniq -c | grep 'helitronORF' | sed "s/^\s*//g" | cut -f1 -d ' ')
Mariner_ant1=$(cat $File | cut -f9 | cut -f2 -d';' | sort| uniq -c | grep 'mariner_ant1' | sed "s/^\s*//g" | cut -f1 -d ' ')
ISC1316=$(cat $File | cut -f9 | cut -f2 -d';' | sort| uniq -c | grep 'ISC1316' | sed "s/^\s*//g" | cut -f1 -d ' ')
Crypton=$(cat $File | cut -f9 | cut -f2 -d';' | sort| uniq -c | grep 'Crypton' | sed "s/^\s*//g" | cut -f1 -d ' ')
printf "$Organism\t$Strain\t$DDE_1\t$Gypsy\t$HAT\t$TY1_Copia\t$Mariner\t$Cacta\t$LINE\t$MuDR_A_B\t$HelitronORF\t$Mariner_ant1\t$ISC1316\t$Crypton\n"
done
```

This table was summaried in an excel spreadsheet and imported into R for further analysis.

# Gene Prediction
Gene prediction followed two steps:
Pre-gene prediction - Quality of genome assemblies were assessed using Cegma to see how many core eukaryotic genes can be identified.
Gene models were used to predict genes in the Alternaria genomes.

## Pre-gene prediction
Quality of genome assemblies was assessed by looking for the gene space in the assemblies.
<!--
```bash
  for Assembly in $(ls repeat_masked/*/*/ncbi_edits_repmask/*_contigs_unmasked.fa); do
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
    echo "$Organism"
    echo "$Strain"
  	ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/cegma
  	qsub $ProgDir/sub_cegma.sh $Assembly dna
  done
``` -->

Busco has replaced CEGMA and was run to check gene space in assemblies

```bash
  for Assembly in $(ls repeat_masked/*/*/ncbi_edits_repmask/*_contigs_unmasked.fa); do
    Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
    Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
    echo "$Organism - $Strain"
    ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
    BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/ascomycota_odb9)
    OutDir=gene_pred/busco/$Organism/$Strain/assembly
    qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
  done
```

```bash
  for File in $(ls gene_pred/busco/*/*/assembly/*/short_summary_*.txt); do
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

```
A.alternata_ssp._arborescens	675	1299	1296	6	10	1315
A.alternata_ssp._arborescens	97.0013	1299	1297	6	10	1315
A.alternata_ssp._arborescens	97.0016	1301	1298	3	11	1315
A.alternata_ssp._gaisen	650	1298	1296	4	13	1315
A.alternata_ssp._tenuissima	1082	1299	1296	4	12	1315
A.alternata_ssp._tenuissima	1164	1302	1299	3	10	1315
A.alternata_ssp_tenuissima	1166	1302	1301	3	10	1315
A.alternata_ssp._tenuissima	1166	1302	1299	5	8	1315
A.alternata_ssp._tenuissima	1177	1303	1299	3	9	1315
A.alternata_ssp._tenuissima	24350	1304	1302	2	9	1315
A.alternata_ssp._tenuissima	635	1303	1299	3	9	1315
A.alternata_ssp._tenuissima	648	1304	1302	3	8	1315
A.alternata_ssp._tenuissima	743	1303	1300	4	8	1315
A.gaisen	650	1298	1296	5	12	1315
Alternaria_alternata	ATCC11680	1300	1298	5	10	1315
Alternaria_alternata	ATCC66891	1299	1297	7	9	1315
Alternaria_alternata	BMP0270	1302	1300	4	9	1315
Alternaria_arborescens	BMP0308	1013	1011	203	99	1315
Alternaria_brassicicola	ATCC96836	1132	1125	111	72	1315
Alternaria_capsici	BMP0180	1290	1289	9	16	1315
Alternaria_carthami	BMP1963	1261	1259	25	29	1315
Alternaria_citriarbusti	BMP2343	1264	1262	25	26	1315
Alternaria_crassa	BMP0172	1286	1284	13	16	1315
Alternaria_dauci	BMP0167	1011	1007	133	171	1315
Alternaria_destruens	BMP0317	439	438	480	396	1315
Alternaria_fragariae	BMP3062	1268	1266	24	23	1315
Alternaria_gaisen	BMP2338	1052	1051	157	106	1315
Alternaria_limoniasperae	BMP2335	1261	1259	26	28	1315
Alternaria_longipes	BMP0313	1299	1297	6	10	1315
Alternaria_macrospora	BMP1949	1212	1211	57	46	1315
Alternaria_mali	BMP3063	1153	1152	87	75	1315
Alternaria_mali	BMP3064	1233	1231	41	41	1315
Alternaria_porri	BMP0178	859	858	244	212	1315
Alternaria_solani	BMP0185	1178	1177	94	43	1315
Alternaria_tagetica	BMP0179	1294	1291	9	12	1315
Alternaria_tangelonis	BMP2327	1221	1219	48	46	1315
Alternaria_tenuissima	BMP0304	1302	1300	4	9	1315
Alternaria_tomatophila	BMP2032	1217	1215	68	30	1315
Alternaria_turkisafria	BMP3436	1211	1209	50	54	1315
```

#Gene prediction

Gene prediction was performed for Alternaria genomes. Two gene prediction
approaches were used:

Gene prediction using Braker1
Prediction of all putative ORFs in the genome using the ORF finder (atg.pl)
approach.


## Gene prediction 1 - Braker1 gene model training and prediction

Gene prediction was performed using Braker1.

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


#### Aligning


```bash
for Assembly in $(ls repeat_masked/*/*/ncbi_edits_repmask/*_contigs_unmasked.fa); do
Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
for FileF in $(ls qc_rna/paired/*/*/F/*.fastq.gz); do
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
for OutDir in $(ls -d alignment/star/*/* | grep -v -e '1166' -e '650'); do
  Strain=$(echo $OutDir | rev | cut -d '/' -f1 | rev)
  Organism=$(echo $OutDir | rev | cut -d '/' -f2 | rev)
  echo "$Organism - $Strain"
  # For all alignments
  BamFiles=$(ls $OutDir/treatment/*/*.sortedByCoord.out.bam | tr -d '\n' | sed 's/.bam/.bam /g')
  mkdir -p $OutDir/concatenated
  samtools merge -@ 4 -f $OutDir/concatenated/concatenated.bam $BamFiles
done
```


#### Braker prediction

```bash
for Assembly in $(ls repeat_masked/*/*/ncbi_edits_repmask/*_contigs_unmasked.fa | grep -v -e '1166' -e '650'); do
Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
mkdir -p alignment/$Organism/$Strain/concatenated
OutDir=gene_pred/braker/$Organism/"$Strain"_braker
AcceptedHits=$(ls alignment/star/$Organism/$Strain/concatenated/concatenated.bam)
GeneModelName="$Organism"_"$Strain"_braker
rm -r /home/armita/prog/augustus-3.1/config/species/"$Organism"_"$Strain"_braker
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/braker1
qsub $ProgDir/sub_braker_fungi.sh $Assembly $OutDir $AcceptedHits $GeneModelName
done
```

** Number of genes predicted:  **
<!--
Prediction of V.inequalis gene models for tom
```bash
for Assembly in $(ls ../venturia/repeat_masked/*/*/filtered_contigs_repmask/*_contigs_unmasked.fa | grep '172_pacbio'); do
Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
mkdir -p alignment/$Organism/$Strain/concatenated
OutDir=gene_pred/braker/$Organism/"$Strain"_braker
# AcceptedHits=$(ls alignment/star/$Organism/$Strain/concatenated/concatenated.bam)
AcceptedHits=$(ls ../venturia/alignment/repeat_masked/v.inaequalis/172_pacbio/concatenated/concatenated.bam)
GeneModelName="$Organism"_"$Strain"_braker
rm -r /home/armita/prog/augustus-3.1/config/species/"$Organism"_"$Strain"_braker
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/braker1
qsub $ProgDir/sub_braker_fungi.sh $Assembly $OutDir $AcceptedHits $GeneModelName
done
``` -->


## Supplimenting Braker gene models with CodingQuary genes

Additional genes were added to Braker gene predictions, using CodingQuary in
pathogen mode to predict additional regions.

Firstly, aligned RNAseq data was assembled into transcripts using Cufflinks.

Note - cufflinks doesn't always predict direction of a transcript and
therefore features can not be restricted by strand when they are intersected.

```bash
for Assembly in $(ls repeat_masked/*/*/ncbi_edits_repmask/*_contigs_unmasked.fa | grep -v -e '1166' -e '650'); do
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
  for Assembly in $(ls repeat_masked/*/*/ncbi_edits_repmask/*_contigs_unmasked.fa | grep -v -e '1166' -e '650'); do
    Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
    Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
    echo "$Organism - $Strain"
    OutDir=gene_pred/codingquary/$Organism/$Strain
    CufflinksGTF=$(ls gene_pred/cufflinks/$Organism/$Strain/concatenated/transcripts.gtf)
    ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/codingquary
    qsub $ProgDir/sub_CodingQuary.sh $Assembly $CufflinksGTF $OutDir
  done
```


Then, additional transcripts were added to Braker gene models, when CodingQuary
genes were predicted in regions of the genome, not containing Braker gene
models:

Note - Ensure that the "TPSI_appended.fa" assembly file is correct.

```bash
for BrakerGff in $(ls gene_pred/braker/*/*_braker/*/augustus.gff3 | grep -v -e '1166' -e '650'); do
Strain=$(echo $BrakerGff| rev | cut -d '/' -f3 | rev | sed 's/_braker//g')
Organism=$(echo $BrakerGff | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
Assembly=$(ls repeat_masked/$Organism/$Strain/ncbi_edits_repmask/*_softmasked_repeatmasker_TPSI_appended.fa)
CodingQuaryGff=gene_pred/codingquary/$Organism/$Strain/out/PredictedPass.gff3
PGNGff=gene_pred/codingquary/$Organism/$Strain/out/PGN_predictedPass.gff3
AddDir=gene_pred/codingquary/$Organism/$Strain/additional
FinalDir=gene_pred/final/$Organism/$Strain/final
AddGenesList=$AddDir/additional_genes.txt
AddGenesGff=$AddDir/additional_genes.gff
FinalGff=$AddDir/combined_genes.gff
mkdir -p $AddDir
mkdir -p $FinalDir

bedtools intersect -v -a $CodingQuaryGff -b $BrakerGff | grep 'gene'| cut -f2 -d'=' | cut -f1 -d';' > $AddGenesList
bedtools intersect -v -a $PGNGff -b $BrakerGff | grep 'gene'| cut -f2 -d'=' | cut -f1 -d';' >> $AddGenesList
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation
$ProgDir/gene_list_to_gff.pl $AddGenesList $CodingQuaryGff CodingQuarry_v2.0 ID CodingQuary > $AddGenesGff
$ProgDir/gene_list_to_gff.pl $AddGenesList $PGNGff PGNCodingQuarry_v2.0 ID CodingQuary >> $AddGenesGff
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


GffBraker=$FinalDir/final_genes_Braker.gff3
GffQuary=$FinalDir/final_genes_CodingQuary.gff3
GffAppended=$FinalDir/final_genes_appended.gff3
cat $GffBraker $GffQuary > $GffAppended
done
```

In preperation for submission to ncbi, gene models were renamed and duplicate gene features were identified and removed.
 * no duplicate genes were identified

Codingquary was noted to predict a gene that went beyond the end of contig 47 in
isolate 1177.

As such this gene was removed manually:

```bash
GffAppended=$(ls gene_pred/final/*/*/final/final_genes_appended.gff3 | grep '1177')
cp $GffAppended tmp.gff
cat tmp.gff | grep -v 'CUFF_8208_1_74' > $GffAppended
```


```bash
  for GffAppended in $(ls gene_pred/final/*/*/final/final_genes_appended.gff3 | grep -v -e '1166' -e '650'); do
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
    Assembly=$(ls repeat_masked/$Organism/$Strain/ncbi_edits_repmask/*_softmasked_repeatmasker_TPSI_appended.fa)
    $ProgDir/gff2fasta.pl $Assembly $GffRenamed gene_pred/final/$Organism/$Strain/final/final_genes_appended_renamed
    # The proteins fasta file contains * instead of Xs for stop codons, these should
    # be changed
    sed -i 's/\*/X/g' gene_pred/final/$Organism/$Strain/final/final_genes_appended_renamed.pep.fasta
  done
```

No duplicate genes were found.
```bash
  for GffAppended in $(ls gene_pred/final/*/*/final/final_genes_appended.gff3); do
    Strain=$(echo $GffAppended | rev | cut -d '/' -f3 | rev)
    Organism=$(echo $GffAppended | rev | cut -d '/' -f4 | rev)
    Gene=$(cat $GffAppended | grep -w 'gene' | wc -l)
    Protein=$(cat $GffAppended | grep -w 'mRNA' | wc -l)
    Augustus=$(cat $GffAppended | grep -w 'gene' | grep 'AUGUSTUS' | wc -l)
    CodingQuary=$(cat $GffAppended | grep -w 'gene' | grep 'CodingQuarry_v2.0' | wc -l)
    printf "$Organism\t$Strain\t$Gene\t$Protein\t$Augustus\t$CodingQuary\n"
  done
```

```
A.alternata_ssp._arborescens	675	12896	12936	12031	865
A.alternata_ssp._arborescens	97.0013	12766	12813	11974	792
A.alternata_ssp._arborescens	97.0016	12820	12863	11988	832
A.alternata_ssp._gaisen	650	13176	13233	12409	767
A.alternata_ssp._tenuissima	1082	13028	13091	12345	683
A.alternata_ssp._tenuissima	1164	13114	13169	12330	784
A.alternata_ssp._tenuissima	1166	13524	13592	12900	624
A.alternata_ssp._tenuissima	1177	13580	13647	12707	873
A.alternata_ssp._tenuissima	24350	12806	12856	12118	688
A.alternata_ssp._tenuissima	635	13733	13812	13049	684
A.alternata_ssp._tenuissima	648	12757	12798	12036	721
A.alternata_ssp._tenuissima	743	13707	13776	12991	716
```

## Checking gene prediction accuracy using BUSCO



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

```
A.alternata_ssp._arborescens	675	1287	1285	20	8	1315
A.alternata_ssp._arborescens	97.0013	1291	1290	17	7	1315
A.alternata_ssp._arborescens	97.0016	1289	1287	18	8	1315
A.alternata_ssp._gaisen	650	1292	1291	18	5	1315
A.alternata_ssp._tenuissima	1082	1293	1291	17	5	1315
A.alternata_ssp._tenuissima	1164	1297	1295	12	6	1315
A.alternata_ssp._tenuissima	1166	1291	1289	18	6	1315
A.alternata_ssp._tenuissima	1177	1288	1285	19	8	1315
A.alternata_ssp._tenuissima	24350	1295	1294	12	8	1315
A.alternata_ssp._tenuissima	635	1293	1290	13	9	1315
A.alternata_ssp._tenuissima	648	1293	1292	13	9	1315
A.alternata_ssp._tenuissima	743	1293	1291	15	7	1315
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
	for Genes in $(ls gene_pred/final/*/*/final/final_genes_appended_renamed.pep.fasta | grep '1177'); do
	echo $Genes
	$ProgDir/sub_interproscan.sh $Genes
	done 2>&1 | tee -a interproscan_submisison.log
```

Following interproscan annotation split files were combined using the following
commands:

```bash
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/interproscan
for Proteins in $(ls gene_pred/final/*/*/final/final_genes_appended_renamed.pep.fasta | grep '1177'); do
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
for Proteome in $(ls gene_pred/final/*/*/final/final_genes_appended_renamed.pep.fasta | grep '1177'); do
Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
OutDir=gene_pred/swissprot/$Organism/$Strain
SwissDbDir=../../uniprot/swissprot
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
  for Proteome in $(ls gene_pred/final/*/*/final/final_genes_appended_renamed.pep.fasta); do
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
  for Proteome in $(ls gene_pred/final/*/*/final/final_genes_appended_renamed.pep.fasta); do
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
  A.alternata_ssp._arborescens	675	1475	1202	1199
  A.alternata_ssp._arborescens	97.0013	1466	1191	1189
  A.alternata_ssp._arborescens	97.0016	1448	1187	1185
  A.alternata_ssp._gaisen	650	1467	1216	1214
  A.alternata_ssp._tenuissima	1082	1510	1248	1246
  A.alternata_ssp._tenuissima	1164	1503	1238	1235
  A.alternata_ssp._tenuissima	1166	1510	1256	1251
  A.alternata_ssp._tenuissima	1177	1537	1272	1270
  A.alternata_ssp._tenuissima	24350	1485	1227	1225
  A.alternata_ssp._tenuissima	635	1539	1267	1261
  A.alternata_ssp._tenuissima	648	1483	1229	1228
  A.alternata_ssp._tenuissima	743	1525	1254	1247
```

### C) From Augustus gene models - Effector identification using EffectorP

Required programs:
 * EffectorP.py

```bash
  for Proteome in $(ls gene_pred/final/*/*/final/final_genes_appended_renamed.pep.fasta); do
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
A.alternata_ssp._arborescens - 675
EffectorP headers:	3844
Secreted effectorP headers:	236
A.alternata_ssp._arborescens - 97.0013
EffectorP headers:	3776
Secreted effectorP headers:	229
A.alternata_ssp._arborescens - 97.0016
EffectorP headers:	3812
Secreted effectorP headers:	229
A.alternata_ssp._gaisen - 650
EffectorP headers:	4012
Secreted effectorP headers:	247
A.alternata_ssp._tenuissima - 1082
EffectorP headers:	3814
Secreted effectorP headers:	252
A.alternata_ssp._tenuissima - 1164
EffectorP headers:	3824
Secreted effectorP headers:	252
A.alternata_ssp._tenuissima - 1166
EffectorP headers:	4032
Secreted effectorP headers:	257
A.alternata_ssp._tenuissima - 1177
EffectorP headers:	4170
Secreted effectorP headers:	268
A.alternata_ssp._tenuissima - 24350
EffectorP headers:	3644
Secreted effectorP headers:	241
A.alternata_ssp._tenuissima - 635
EffectorP headers:	4122
Secreted effectorP headers:	264
A.alternata_ssp._tenuissima - 648
EffectorP headers:	3598
Secreted effectorP headers:	233
A.alternata_ssp._tenuissima - 743
EffectorP headers:	4124
Secreted effectorP headers:	263
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
# printf "Number of effectors predicted by EffectorP:\t"
# EffectorP=$(ls analysis/effectorP/$Organism/$Strain/*_EffectorP_secreted_headers.txt)
# cat $EffectorP | wc -l
# printf "Number of SSCPs predicted by both effectorP and this approach: \t"
# cat $OutDir/"$Strain"_sscp.fa | grep '>' | tr -d '>' > $OutDir/"$Strain"_sscp_headers.txt
# cat $OutDir/"$Strain"_sscp_headers.txt $EffectorP | cut -f1 | sort | uniq -d | wc -l
# echo ""
done
```

```
A.alternata_ssp._arborescens - 675
% cysteine content threshold set to:	3
maximum length set to:	300
No. short-cysteine rich proteins in input fasta:	196
number of SSC-rich genes:	196
A.alternata_ssp._arborescens - 97.0013
% cysteine content threshold set to:	3
maximum length set to:	300
No. short-cysteine rich proteins in input fasta:	187
number of SSC-rich genes:	187
A.alternata_ssp._arborescens - 97.0016
% cysteine content threshold set to:	3
maximum length set to:	300
No. short-cysteine rich proteins in input fasta:	190
number of SSC-rich genes:	190
A.alternata_ssp._gaisen - 650
% cysteine content threshold set to:	3
maximum length set to:	300
No. short-cysteine rich proteins in input fasta:	201
number of SSC-rich genes:	201
A.alternata_ssp._tenuissima - 1082
% cysteine content threshold set to:	3
maximum length set to:	300
No. short-cysteine rich proteins in input fasta:	199
number of SSC-rich genes:	199
A.alternata_ssp._tenuissima - 1164
% cysteine content threshold set to:	3
maximum length set to:	300
No. short-cysteine rich proteins in input fasta:	202
number of SSC-rich genes:	202
A.alternata_ssp._tenuissima - 1166
% cysteine content threshold set to:	3
maximum length set to:	300
No. short-cysteine rich proteins in input fasta:	204
number of SSC-rich genes:	204
A.alternata_ssp._tenuissima - 1177
% cysteine content threshold set to:	3
maximum length set to:	300
No. short-cysteine rich proteins in input fasta:	217
number of SSC-rich genes:	217
A.alternata_ssp._tenuissima - 24350
% cysteine content threshold set to:	3
maximum length set to:	300
No. short-cysteine rich proteins in input fasta:	190
number of SSC-rich genes:	190
A.alternata_ssp._tenuissima - 635
% cysteine content threshold set to:	3
maximum length set to:	300
No. short-cysteine rich proteins in input fasta:	212
number of SSC-rich genes:	212
A.alternata_ssp._tenuissima - 648
% cysteine content threshold set to:	3
maximum length set to:	300
No. short-cysteine rich proteins in input fasta:	184
number of SSC-rich genes:	184
A.alternata_ssp._tenuissima - 743
% cysteine content threshold set to:	3
maximum length set to:	300
No. short-cysteine rich proteins in input fasta:	211
number of SSC-rich genes:	211
```

## CAZY proteins

Carbohydrte active enzymes were idnetified using CAZYfollowing recomendations
at http://csbl.bmb.uga.edu/dbCAN/download/readme.txt :

```bash
for Proteome in $(ls gene_pred/final/*/*/final/final_genes_appended_renamed.pep.fasta); do
Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
OutDir=gene_pred/CAZY/$Organism/$Strain
mkdir -p $OutDir
Prefix="$Strain"_CAZY
CazyHmm=../../dbCAN/dbCAN-fam-HMMs.txt
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
A.alternata_ssp._arborescens	675	766	766	375	375
A.alternata_ssp._arborescens	97.0013	753	753	372	372
A.alternata_ssp._arborescens	97.0016	767	767	382	382
A.alternata_ssp._gaisen	650	780	780	383	383
A.alternata_ssp._tenuissima	1082	787	787	401	401
A.alternata_ssp._tenuissima	1164	774	774	385	385
A.alternata_ssp._tenuissima	1166	783	783	384	384
A.alternata_ssp._tenuissima	1177	791	791	397	397
A.alternata_ssp._tenuissima	24350	770	770	391	391
A.alternata_ssp._tenuissima	635	791	791	397	397
A.alternata_ssp._tenuissima	648	773	773	392	392
A.alternata_ssp._tenuissima	743	788	788	390	390
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
Gff=gene_pred/final/$Organism/$Strain/final/final_genes_appended_renamed.gff3
ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/CAZY
$ProgDir/summarise_CAZY.py --cazy $CAZY --inp_secreted $Secreted --inp_gff $Gff --summarise_family --trim_gene_id 2 --kubicek_2014
done
```

```
A.alternata_ssp._arborescens - 675
B-Galactosidases - 4
B-Glucuronidases - 3
Polygalacturonase - 8
A-Arabinosidases - 12
Xylanases - 13
Polygalacturonate lyases - 22
A-Galactosidases - 2
B-Glycosidases - 12
Cellulases - 29
other - 268
Xyloglucanases - 1

A.alternata_ssp._arborescens - 97.0013
B-Galactosidases - 4
B-Glucuronidases - 3
Polygalacturonase - 10
A-Arabinosidases - 12
Xylanases - 12
Polygalacturonate lyases - 20
A-Galactosidases - 2
B-Glycosidases - 10
Cellulases - 31
other - 266
Xyloglucanases - 1

A.alternata_ssp._arborescens - 97.0016
B-Galactosidases - 4
A-Galactosidases - 3
Polygalacturonase - 8
A-Arabinosidases - 13
Xylanases - 14
Polygalacturonate lyases - 22
B-Glucuronidases - 3
B-Glycosidases - 13
Cellulases - 29
other - 271
Xyloglucanases - 1

A.alternata_ssp._gaisen - 650
B-Galactosidases - 4
B-Glucuronidases - 4
Polygalacturonase - 12
A-Arabinosidases - 14
Xylanases - 14
Polygalacturonate lyases - 21
A-Galactosidases - 3
B-Glycosidases - 11
Cellulases - 29
other - 269
Xyloglucanases - 1

A.alternata_ssp._tenuissima - 1082
B-Galactosidases - 3
B-Glucuronidases - 3
Polygalacturonase - 8
A-Arabinosidases - 15
Xylanases - 14
Polygalacturonate lyases - 21
A-Galactosidases - 3
B-Glycosidases - 13
Cellulases - 31
other - 288
Xyloglucanases - 1

A.alternata_ssp._tenuissima - 1164
B-Galactosidases - 3
B-Glucuronidases - 3
Polygalacturonase - 10
A-Arabinosidases - 15
Xylanases - 14
Polygalacturonate lyases - 21
A-Galactosidases - 2
B-Glycosidases - 13
Cellulases - 31
other - 272

A.alternata_ssp._tenuissima - 1166
B-Galactosidases - 4
A-Galactosidases - 3
Polygalacturonase - 8
A-Arabinosidases - 15
Xylanases - 14
Polygalacturonate lyases - 20
B-Glucuronidases - 3
B-Glycosidases - 12
Cellulases - 31
other - 271
Xyloglucanases - 1

A.alternata_ssp._tenuissima - 1177
B-Galactosidases - 4
A-Galactosidases - 2
Polygalacturonase - 9
A-Arabinosidases - 14
Xylanases - 14
Polygalacturonate lyases - 22
B-Glucuronidases - 4
B-Glycosidases - 12
Cellulases - 30
other - 284
Xyloglucanases - 1

A.alternata_ssp._tenuissima - 24350
B-Galactosidases - 3
A-Galactosidases - 2
Polygalacturonase - 9
A-Arabinosidases - 15
Xylanases - 14
Polygalacturonate lyases - 21
B-Glucuronidases - 4
B-Glycosidases - 12
B-Mannase - 1
Cellulases - 31
other - 276
Xyloglucanases - 1

A.alternata_ssp._tenuissima - 635
B-Galactosidases - 4
A-Galactosidases - 2
Polygalacturonase - 8
A-Arabinosidases - 14
Xylanases - 14
Polygalacturonate lyases - 22
B-Glucuronidases - 3
B-Glycosidases - 13
Cellulases - 29
other - 285

A.alternata_ssp._tenuissima - 648
B-Galactosidases - 4
B-Glucuronidases - 3
Polygalacturonase - 8
A-Arabinosidases - 15
Xylanases - 13
Polygalacturonate lyases - 22
A-Galactosidases - 3
B-Glycosidases - 12
Cellulases - 31
other - 279
Xyloglucanases - 1

A.alternata_ssp._tenuissima - 743
B-Galactosidases - 4
A-Galactosidases - 2
Polygalacturonase - 8
A-Arabinosidases - 14
Xylanases - 14
Polygalacturonate lyases - 22
B-Glucuronidases - 3
B-Glycosidases - 13
Cellulases - 29
other - 278
```


## D) Secondary metabolites (Antismash and SMURF)

Antismash was run to identify clusters of secondary metabolite genes within
the genome. Antismash was run using the webserver at:
http://antismash.secondarymetabolites.org


Results of web-annotation of gene clusters within the assembly were downloaded to
the following directories:

```bash
  for Assembly in $(ls repeat_masked/*/*/ncbi_edits_repmask*/*_contigs_softmasked_repeatmasker_TPSI_appended.fa); do
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    OutDir=analysis/secondary_metabolites/antismash/$Organism/$Strain
    mkdir -p $OutDir
  done
```

```bash
  for Zip in $(ls analysis/secondary_metabolites/antismash/*/*/*.zip); do
    OutDir=$(dirname $Zip)
    unzip -d $OutDir $Zip
  done
```

```bash
  for AntiSmash in $(ls analysis/secondary_metabolites/antismash/*/*/*/*.final.gbk); do
    Organism=$(echo $AntiSmash | rev | cut -f4 -d '/' | rev)
    Strain=$(echo $AntiSmash | rev | cut -f3 -d '/' | rev)
    echo "$Organism - $Strain"
    OutDir=analysis/secondary_metabolites/antismash/$Organism/$Strain
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
```

SMURF was also run to identify secondary metabolite gene clusters.

Genes needed to be parsed into a specific tsv format prior to submission on the
SMURF webserver.

```bash
for Gff in $(ls gene_pred/final/*/*/final/final_genes_appended_renamed.gff3); do
	Organism=$(echo $Gff | rev | cut -f4 -d '/' | rev)
	Strain=$(echo $Gff | rev | cut -f3 -d '/' | rev)
	echo "$Organism - $Strain"
  OutDir=analysis/secondary_metabolites/smurf/$Organism/$Strain
  mkdir -p $OutDir
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/secondary_metabolites
  $ProgDir/gff2smurf.py --gff $Gff > $OutDir/"$Strain"_genes_smurf.tsv
  cp $OutDir/"$Strain"_genes_smurf.tsv download/.
done
```

SMURF results were observed and found to be untrustworthy as the highly fragmented
assembly meant that the order of gene models was not interpreted correctly.

<!--
SMURF output was received by email and downloaded to the cluster in the output
directory above.

Output files were parsed into gff format:

```bash
  OutDir=analysis/secondary_metabolites/smurf/F.venenatum/WT
  Prefix="WT"
  GeneGff=gene_pred/final/F.venenatum/WT/final/final_genes_appended_renamed.gff3
  SmurfClusters=$OutDir/Secondary-Metabolite-Clusters.txt
  SmurfBackbone=$OutDir/Backbone-genes.txt
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/secondary_metabolites
  $ProgDir/smurf2gff.py --smurf_clusters $SmurfClusters --smurf_backbone $SmurfBackbone > $OutDir/Smurf_clusters.gff
  bedtools intersect -wo -a $GeneGff -b $OutDir/Smurf_clusters.gff | grep 'mRNA' | cut -f9,10,12,18 | sed "s/ID=//g" | perl -p -i -e "s/;Parent=g\w+//g" | perl -p -i -e "s/;Notes=.*//g" > $OutDir/"$Prefix"_smurf_secmet_genes.tsv
```

Total number of secondary metabolite clusters:

```bash
for Assembly in $(ls repeat_masked/*/*/illumina_assembly_ncbi/*_contigs_softmasked_repeatmasker_TPSI_appended.fa | grep -w 'WT'); do
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
OutDir=analysis/secondary_metabolites/antismash/$Organism/$Strain
mkdir -p $OutDir
GeneGff=$(ls gene_pred/final/$Organism/$Strain/final/final_genes_appended_renamed.gff3)
AntismashClusters=$(ls analysis/secondary_metabolites/antismash/$Organism/$Strain/*_secmet_clusters.gff)
SmurfClusters=$(ls analysis/secondary_metabolites/smurf/$Organism/$Strain/Smurf_clusters.gff)
echo "Total number of Antismash clusters"
cat $AntismashClusters | wc -l
echo "Total number of SMURF clusters"
cat $SmurfClusters | wc -l
echo "number of Antismash clusters intersecting Smurf clusters"
bedtools intersect -a $AntismashClusters -b $SmurfClusters | wc -l
echo "number of Antismash clusters not intersecting Smurf clusters"
bedtools intersect -v -a $AntismashClusters -b $SmurfClusters | wc -l
echo "number of smurf clusters intersecting antismash clusters"
bedtools intersect -a $SmurfClusters -b $AntismashClusters | wc -l
echo "number of smurf clusters not intersecting antismash clusters"
bedtools intersect -v -a $SmurfClusters -b $AntismashClusters | wc -l
done
``` -->



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
A.alternata_ssp._arborescens	675	605
A.alternata_ssp._arborescens	97.0013	601
A.alternata_ssp._arborescens	97.0016	593
A.alternata_ssp._gaisen	650	640
A.alternata_ssp._tenuissima	1082	640
A.alternata_ssp._tenuissima	1164	627
A.alternata_ssp._tenuissima	1166	660
A.alternata_ssp._tenuissima	1177	647
A.alternata_ssp._tenuissima	24350	619
A.alternata_ssp._tenuissima	635	677
A.alternata_ssp._tenuissima	648	606
A.alternata_ssp._tenuissima	743	676
```




# Genomic analysis
The first analysis was based upon BLAST searches for genes known to be involved in toxin production


## Genes with homology to PHIbase
Predicted gene models were searched against the PHIbase database using tBLASTx.

<!-- ```bash
  for Subject in $(ls repeat_masked/*/*/ncbi_edits_repmask/*_contigs_unmasked.fa); do
    ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
    Query=../../phibase/v4.2/PHI_accessions.fa
    qsub $ProgDir/blast_pipe.sh $Query protein $Subject
  done
``` -->


```bash
qlogin -l h=blacklace11.blacklace
cd /home/groups/harrisonlab/project_files/alternaria
dbFasta=$(ls ../../phibase/v4.4/phi_accessions.fa)
dbType="prot"
for QueryFasta in $(ls gene_pred/final/*/*/final/final_genes_appended_renamed.cds.fasta); do
Organism=$(echo $QueryFasta | rev | cut -f4 -d '/' | rev)
Strain=$(echo $QueryFasta | rev | cut -f3 -d '/' | rev)
echo "$Organism - $Strain"
Prefix="${Strain}_phi_accessions"
Eval="1e-30"
# WorkDir=$TMPDIR
OutDir=analysis/blast_homology/$Organism/$Strain
mkdir -p $OutDir
makeblastdb -in $dbFasta -input_type fasta -dbtype $dbType -title $Prefix.db -parse_seqids -out $OutDir/$Prefix.db
blastx -num_threads 20 -db $OutDir/$Prefix.db -query $QueryFasta -outfmt 6 -num_alignments 1 -out $OutDir/${Prefix}_hits.txt -evalue $Eval
cat $OutDir/${Prefix}_hits.txt | grep 'effector' | cut -f1,2 | sort | uniq > $OutDir/${Prefix}_hits_headers.txt
done
```



## Presence and genes with homology to Alternaria toxins

The first analysis was based upon BLAST searches for genes known to be involved in toxin production



```bash
  for Subject in $(ls repeat_masked/*/*/ncbi_edits_repmask/*_contigs_unmasked.fa); do
    ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
    Query=analysis/blast_homology/CDC_genes/A.alternata_CDC_genes.fa
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
  done
  ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
  $ProgDir/blast_parse.py --blast_csv $HitsList --headers $StrainList --identity 0.7 --evalue 1e-30
```


Once blast searches had completed, the BLAST hits were converted to GFF
annotations:

```bash
  for BlastHits in $(ls analysis/blast_homology/*/*/*_A.alternata_CDC_genes.fa_homologs.csv); do
    ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
    Strain=$(echo $BlastHits | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $BlastHits | rev | cut -f3 -d '/' | rev)
    HitsGff=analysis/blast_homology/$Organism/$Strain/"$Strain"_A.alternata_CDC_genes.fa_homologs.gff
    Column2=toxin_homolog
    NumHits=2
    $ProgDir/blast2gff.pl $Column2 $NumHits $BlastHits > $HitsGff
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
NoIntersectBlasthost
=$OutDir/"$Strain"_A.alternata_CDC_genes_NoIntersect.gff
mkdir -p $OutDir
echo "$Organism - $Strain"
echo "The number of BLAST hits in the gff file were:"
cat $HitsGff | wc -l
echo "The number of blast hits intersected were:"
bedtools intersect -wao -a $HitsGff -b $Proteins > $IntersectBlast
cat $IntersectBlast | grep -w 'gene' | wc -l
echo "The number of blast hits not intersecting gene models were:"
bedtools intersect -v -a $HitsGff -b $Proteins > $NoIntersectBlast
# rm $NoIntersectBlast
cat $IntersectBlast | grep -w -E '0$' | wc -l
done
```

As a proof of concept the gene intersecting the AMT2 and AMT4 homolog was
searched for in the results of the orthology analysis. The AMT2 homolog
was present in a large orthology group, containing genes from the core genome -
present in all A. tenuissima isolates. This explains why this locus may not be
approapriate as a marker for molecular identification of apple pathotypes. The
AMT2 group did not appear to be expanded in pathotypic isolates. The AMT4 locus
was identified in an ortholog group specific to apple pathotype isolates:

```bash
# AMT2
cat analysis/orthology/orthomcl/At_Aa_Ag_all_isolates/At_Aa_Ag_all_isolates_orthogroups.txt | grep 'At_6|g4087' | less
cat analysis/orthology/orthomcl/At_Aa_Ag_all_isolates/At_Aa_Ag_all_isolates_orthogroups.txt | grep 'At_6|g4087' | sed 's/\s/\n/g' | grep -v 'orthogroup' |  cut -f1 -d '|' | sort | uniq -c | less
# AMT4
cat analysis/orthology/orthomcl/At_Aa_Ag_all_isolates/At_Aa_Ag_all_isolates_orthogroups.txt | grep 'At_6|g4038' | less
```

<!--
```bash
qlogin -l h=blacklace01.blacklace -pe smp 12
cd /home/groups/harrisonlab/project_files/alternaria
dbFasta=$(ls analysis/blast_homology/CDC_genes/A.alternata_CDC_genes.fa)
dbType="nucl"
for QueryFasta in $(ls gene_pred/final/*/*/final/final_genes_appended_renamed.cds.fasta | grep -v '650'); do
Organism=$(echo $QueryFasta | rev | cut -f4 -d '/' | rev)
Strain=$(echo $QueryFasta | rev | cut -f3 -d '/' | rev)
echo "$Organism - $Strain"
Prefix="${Strain}_CDC_genes"
Eval="1e-30"
# WorkDir=$TMPDIR
OutDir=analysis/blast_homology/$Organism/$Strain
mkdir -p $OutDir
makeblastdb -in $dbFasta -input_type fasta -dbtype $dbType -title $Prefix.db -parse_seqids -out $OutDir/$Prefix.db
tblastx -num_threads 12 -db $OutDir/$Prefix.db -query $QueryFasta -outfmt 6 -num_alignments 10 -out $OutDir/${Prefix}_hits.txt -evalue $Eval
cat $OutDir/${Prefix}_hits.txt | cut -f1,2 | sort | uniq > $OutDir/${Prefix}_hits_headers.txt
done
``` -->




# Build Annotation Tables


```bash
for GeneGff in $(ls gene_pred/final/*/*/final/final_genes_appended_renamed.gff3); do
Strain=$(echo $GeneGff | rev | cut -f3 -d '/' | rev)
Organism=$(echo $GeneGff | rev | cut -f4 -d '/' | rev)
Assembly=$(ls repeat_masked/$Organism/$Strain/ncbi_edits_repmask/*_contigs_unmasked.fa)
# GeneConversions=$(ls genome_submission/$Organism/$Strain/gag/edited/genome_gene_conversions.tsv)
Antismash=$(ls analysis/secondary_metabolites/antismash/$Organism/$Strain/*/geneclusters.txt)
# Smurf=$(ls analysis/secondary_metabolites/smurf/F.venenatum/WT/WT_smurf_secmet_genes.tsv)
TFs=$(ls analysis/transcription_factors/$Organism/$Strain/"$Strain"_TF_domains.tsv)
SigP=$(ls gene_pred/braker_signalp-4.1/$Organism/$Strain/"$Strain"_aug_sp.aa)
TM_out=$(ls gene_pred/trans_mem/$Organism/$Strain/"$Strain"_TM_genes_pos.txt)
# GPI_out=$(ls gene_pred/trans_mem/$Organism/$Strain/GPIsom/GPI_pos.fa)
# MIMP_list=$(ls analysis/mimps/$Organism/$Strain/"$Strain"_genes_in_2kb_mimp.txt)
EffP_list=$(ls analysis/effectorP/$Organism/$Strain/"$Organism"_"$Strain"_EffectorP_headers.txt)
CAZY_list=$(ls gene_pred/CAZY/$Organism/$Strain/"$Strain"_CAZY.out.dm.ps)
PhiHits=$(ls analysis/blast_homology/$Organism/$Strain/"$Strain"_phi_accessions_hits_headers.txt)
ToxinHits=$(ls analysis/blast_homology/$Organism/$Strain/"$Strain"_CDC_genes_hits_headers.txt)
InterPro=$(ls gene_pred/interproscan/$Organism/$Strain/*_interproscan.tsv)
SwissProt=$(ls gene_pred/swissprot/$Organism/$Strain/swissprot_vJul2016_tophit_parsed.tbl)
# Orthology=$(ls analysis/orthology/orthomcl/At_Aa_Ag_all_isolates/At_Aa_Ag_all_isolates_orthogroups.txt)
Orthology=$(ls /data/scratch/armita/alternaria/analysis/orthology/orthomcl/At_Aa_Ag_all_isolates/formatted/Results_May31/Orthogroups.txt)
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
$ProgDir/build_annot_table_Alt.py \
  --genome $Assembly \
  --genes_gff $GeneGff \
  --Antismash $Antismash \
  --TFs $TFs \
  --SigP $SigP \
  --TM_list $TM_out \
  --EffP_list $EffP_list \
  --CAZY_list $CAZY_list \
  --PhiHits $PhiHits \
  --ToxinHits $ToxinHits \
  --InterPro $InterPro \
  --Swissprot $SwissProt \
  --orthogroups $Orthology \
  --strain_id $OrthoStrainID \
  --OrthoMCL_all $OrthoStrainAll \
  > $OutDir/"$Strain"_annotation_ncbi.tsv
# $ProgDir/FoN_build_gene_annot_table.py \
#   --genome $Assembly \
#   --genes_gff $GeneGff \
#   --genes_renamed $GeneConversions \
#   --Antismash $Antismash \
#   --TFs $TFs \
#   --SigP $SigP \
#   --TM_list $TM_out \
#   --GPI_list $GPI_out \
#   --MIMP_list $MIMP_list \
#   --EffP_list $EffP_list \
#   --CAZY_list $CAZY_list \
#   --InterPro $InterPro \
#   --Swissprot $SwissProt \
#   --orthogroups $Orthology \
#   --strain_id $OrthoStrainID \
#   --OrthoMCL_all $OrthoStrainAll \
#   > $OutDir/"$Strain"_annotation_ncbi.tsv
done
```


### NLPs

```bash
for AnnotTab in $(ls gene_pred/annotation/A.*/*/*_annotation_ncbi.tsv); do
echo "$AnnotTab"
cat $AnnotTab | grep 'Necrosis' | wc -l
done
```

Two NLP proteins were identified in each genome

### CAZY proteins:

```bash
for AnnotTab in $(ls gene_pred/annotation/A.*/*/*_annotation_ncbi.tsv); do
echo "$AnnotTab"
Strain=$(echo $AnnotTab| rev | cut -f2 -d '/' | rev)
Organism=$(echo $AnnotTab | rev | cut -f3 -d '/' | rev)
# OutDir=$(dirname $AnnotTab)/subset
OutDir=/data/scratch/armita/alternaria/$(dirname $AnnotTab)/subset
mkdir -p $OutDir
cat $AnnotTab | cut -f1,7,8,9,10,19,20,21 | grep 'CAZY' | wc -l
cat $AnnotTab | cut -f1,7,8,9,10,19,20,21| grep 'CAZY' > $OutDir/10300_gene_table_CAZY.tsv
cat $AnnotTab | cut -f1,7,8,9,10,19,20,21 | grep 'CAZY' | grep 'SigP' | wc -l
cat $AnnotTab | cut -f1,7,8,9,10,19,20,21 | grep 'CAZY' | grep 'SigP' > $OutDir/10300_gene_table_CAZY_secreted.tsv

cat $OutDir/10300_gene_table_CAZY_secreted.tsv | grep 'GBGX' | wc -l
cat $OutDir/10300_gene_table_CAZY_secreted.tsv | grep 'GBGX' | cut -f5 | sort | uniq -c | wc -l


cat $AnnotTab | cut -f1,7,8,9,10,19,20,21 | grep 'CAZY' | grep 'SigP' | cut -f5 | sort | uniq -c | sort -nr > $OutDir/10300_gene_table_CAZY_hmm_models.txt
# cat $OutDir/10300_gene_table_CAZY_hmm_models.txt | sed 's/.hmm//g' | sed 's/CAZY://g' | grep 'GH'
# cat $OutDir/10300_gene_table_CAZY_hmm_models.txt | sed 's/.hmm//g' | sed 's/CAZY://g' | grep -v 'GH' | grep 'CBM'
# cat $OutDir/10300_gene_table_CAZY_hmm_models.txt | sed 's/.hmm//g' | sed 's/CAZY://g' | grep 'AA'
# cat $OutDir/10300_gene_table_CAZY_hmm_models.txt | sed 's/.hmm//g' | sed 's/CAZY://g' | grep 'CE'
# cat $OutDir/10300_gene_table_CAZY_hmm_models.txt | sed 's/.hmm//g' | sed 's/CAZY://g' | grep -v 'GH' | grep 'GT'
# cat $OutDir/10300_gene_table_CAZY_hmm_models.txt | sed 's/.hmm//g' | sed 's/CAZY://g' | grep 'PL'
done
```



## Mating type genes:

```bash
qlogin -pe smp 12
cd /home/groups/harrisonlab/project_files/alternaria
QueryFasta=$(ls /data/scratch/armita/alternaria/analysis/blast/MAT_genes/MAT_genes.fasta)
dbType="nucl"
for dbFasta in $(ls repeat_masked/*/*/ncbi_edits_repmask//*_contigs_softmasked_repeatmasker_TPSI_appended.fa); do
Organism=$(echo $dbFasta | rev | cut -f4 -d '/' | rev)
Strain=$(echo $dbFasta | rev | cut -f3 -d '/' | rev)
echo "$Organism - $Strain"
Prefix="${Strain}_MAT"
Eval="1e-30"
OutDir=analysis/blast_homology/$Organism/$Strain
mkdir -p $OutDir
makeblastdb -in $dbFasta -input_type fasta -dbtype $dbType -title $Prefix.db -parse_seqids -out $OutDir/$Prefix.db
blastn -num_threads 12 -db $OutDir/$Prefix.db -query $QueryFasta -outfmt 6 -num_alignments 1 -out $OutDir/${Prefix}_hits.txt -evalue $Eval
done
```

```bash
for File in $(ls analysis/blast_homology/*/*/*_MAT_hits.txt); do
  Organism=$(echo $File | rev | cut -f3 -d '/' | rev)
  Strain=$(echo $File | rev | cut -f2 -d '/' | rev)
  for Hit in $(cat $File | cut -f1 | cut -f1 -d '_'); do
    printf "$Organism\t$Strain\t$Hit\n"
  done
done
```

```
A.alternata_ssp._arborescens	675	MAT1-2-1
A.alternata_ssp._arborescens	97.0013	MAT1-1-1
A.alternata_ssp._arborescens	97.0016	MAT1-2-1
A.alternata_ssp._gaisen	650	MAT1-1-1
A.alternata_ssp._tenuissima	1082	MAT1-2-1
A.alternata_ssp._tenuissima	1164	MAT1-2-1
A.alternata_ssp._tenuissima	1166	MAT1-2-1
A.alternata_ssp._tenuissima	1177	MAT1-1-1
A.alternata_ssp._tenuissima	24350	MAT1-1-1
A.alternata_ssp._tenuissima	635	MAT1-2-1
A.alternata_ssp._tenuissima	648	MAT1-1-1
A.alternata_ssp._tenuissima	743	MAT1-2-1
Alternaria_alternata	ATCC11680	MAT1-1-1
Alternaria_alternata	ATCC66891	MAT1-1-1
Alternaria_alternata	BMP0270	MAT1-1-1
Alternaria_arborescens	BMP0308	MAT1-2-1
Alternaria_brassicicola	ATCC96836	MAT1-2-1
Alternaria_capsici	BMP0180	MAT1-2-1
Alternaria_carthami	BMP1963	MAT1-1-1
Alternaria_citriarbusti	BMP2343	MAT1-1-1
Alternaria_crassa	BMP0172	MAT1-1-1
Alternaria_dauci	BMP0167	MAT1-2-1
Alternaria_destruens	BMP0317	MAT1-1-1
Alternaria_destruens	BMP0317	MAT1-1-1
Alternaria_fragariae	BMP3062	MAT1-1-1
Alternaria_gaisen	BMP2338	MAT1-1-1
Alternaria_limoniasperae	BMP2335	MAT1-2-1
Alternaria_longipes	BMP0313	MAT1-1-1
Alternaria_macrospora	BMP1949	MAT1-1-1
Alternaria_mali	BMP3063	MAT1-1-1
Alternaria_mali	BMP3064	MAT1-1-1
Alternaria_porri	BMP0178	MAT1-1-1
Alternaria_solani	BMP0185	MAT1-1-1
Alternaria_tagetica	BMP0179	MAT1-1-1
Alternaria_tangelonis	BMP2327	MAT1-2-1
Alternaria_tenuissima	BMP0304	MAT1-2-1
Alternaria_tomatophila	BMP2032	MAT1-2-1
Alternaria_turkisafria	BMP3436	MAT1-1-1
```

# Toxin genes

```bash
qlogin -pe smp 4
cd /home/groups/harrisonlab/project_files/alternaria
QueryFasta=$(ls analysis/blast_homology/CDC_genes/A.alternata_CDC_genes.fa)
dbType="nucl"
for dbFasta in $(ls repeat_masked/*/*/ncbi_edits_repmask/*_contigs_softmasked_repeatmasker_TPSI_appended.fa | grep -v -e '650' -e '1166'); do
Strain=$(echo $dbFasta | rev | cut -f3 -d '/' | rev)
Organism=$(echo $dbFasta | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
Prefix="${Strain}_toxgenes"
Eval="1e-30"
OutDir=analysis/blast_homology/$Organism/$Strain
mkdir -p $OutDir
makeblastdb -in $dbFasta -input_type fasta -dbtype $dbType -title $Prefix.db -parse_seqids -out $OutDir/$Prefix.db
blastn -num_threads 4 -db $OutDir/$Prefix.db -query $QueryFasta -outfmt 6 -num_alignments 1 -out $OutDir/${Prefix}_hits.txt -evalue $Eval
done
```

```bash
for File in $(ls analysis/blast_homology/*/*/*_toxgenes_hits.txt); do
  Organism=$(echo $File | rev | cut -f3 -d '/' | rev)
  Strain=$(echo $File | rev | cut -f2 -d '/' | rev)
  for Hit in $(cat $File | cut -f1 | cut -f1 -d ' '); do
    printf "$Organism\t$Strain\t$Hit\n"
  done
done > analysis/blast_homology/CDC_genes/ref_genome_summarised_hits.tsv
```
