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


kmer counting was performed using kmc.
This allowed estimation of sequencing depth and total genome size:

```bash
  for StrainPath in $(ls -d qc_dna/paired/*/*); do
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc
    TrimF=$(ls $StrainPath/F/*.fq.gz)
    TrimR=$(ls $StrainPath/R/*.fq.gz)
    echo $TrimF
    echo $TrimR
    qsub $ProgDir/kmc_kmer_counting.sh $Trim_F $Trim_R
  done
```

** Estimated Genome Size is: 48490129

** Esimated Coverage is: 35

# Assembly
Assembly was performed using: Velvet / Abyss / Spades

<!-- ## Velvet Assembly
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
``` -->

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
    for Assembly in $(ls assembly/spades/*/*/filtered_contigs/contigs_min_500bp.fasta); do
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
    OutDir=assembly/spades/$Organism/$Strain/filtered_contigs
    qsub $ProgDir/sub_quast.sh $Assembly $OutDir
  done
```

The results of quast were shown using the following commands:

```bash
  for Assembly in $(ls assembly/spades/*/*/filtered_contigs/report.txt); do
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
  for Assembly in $(ls assembly/spades/*/*/filtered_contigs/contigs_min_500bp.fasta); do
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
    OutDir=assembly/spades/$Organism/$Strain/filtered_contigs
    $ProgDir/remove_contaminants.py --inp $Assembly --out $OutDir/contigs_min_500bp_renamed.fasta --coord_file tmp.csv
  done
  rm tmp.csv
```



# Repeat masking
Repeat masking was performed and used the following programs: Repeatmasker Repeatmodeler

The renamed assembly was used to perform repeatmasking.

```bash
  for Assembly in $(ls assembly/spades/*/1177/filtered_contigs/contigs_min_500bp_renamed.fasta); do
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
    echo "$Organism"
    echo "$Strain"
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/repeat_masking
    qsub $ProgDir/rep_modeling.sh $Assembly
    qsub $ProgDir/transposonPSI.sh $Assembly
  done
 ```


 The number of bases masked by transposonPSI and Repeatmasker were summarised
 using the following commands:

 ```bash
  for RepDir in $(ls -d repeat_masked/*/*/filtered_contigs_repmask); do
    Strain=$(echo $RepDir | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $RepDir | rev | cut -f3 -d '/' | rev)  
    RepMaskGff=$(ls $RepDir/*_contigs_hardmasked.gff)
    TransPSIGff=$(ls $RepDir/*_contigs_unmasked.fa.TPSI.allHits.chains.gff3)
    printf "$Organism\t$Strain\n"
    printf "The number of bases masked by RepeatMasker:\t"
    sortBed -i $RepMaskGff | bedtools merge | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'
    printf "The number of bases masked by TransposonPSI:\t"
    sortBed -i $TransPSIGff | bedtools merge | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'
    printf "The total number of masked bases are:\t"
    cat $RepMaskGff $TransPSIGff | sortBed | bedtools merge | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'
    echo
  done
 ```

```bash
 A.alternata_ssp._arborescens	675
 The number of bases masked by RepeatMasker:	679814
 The number of bases masked by TransposonPSI:	325226
 The total number of masked bases are:	778492

 A.alternata_ssp._arborescens	97.0013
 The number of bases masked by RepeatMasker:	841134
 The number of bases masked by TransposonPSI:	355274
 The total number of masked bases are:	972725

 A.alternata_ssp._arborescens	97.0016
 The number of bases masked by RepeatMasker:	555661
 The number of bases masked by TransposonPSI:	337525
 The total number of masked bases are:	712560

 A.alternata_ssp._gaisen	650
 The number of bases masked by RepeatMasker:	329570
 The number of bases masked by TransposonPSI:	213597
 The total number of masked bases are:	467341

 A.alternata_ssp._tenuissima	1082
 The number of bases masked by RepeatMasker:	510714
 The number of bases masked by TransposonPSI:	244175
 The total number of masked bases are:	621605

 A.alternata_ssp._tenuissima	1164
 The number of bases masked by RepeatMasker:	744712
 The number of bases masked by TransposonPSI:	270404
 The total number of masked bases are:	906187

 A.alternata_ssp._tenuissima	1166
 The number of bases masked by RepeatMasker:	709427
 The number of bases masked by TransposonPSI:	262625
 The total number of masked bases are:	855593

 A.alternata_ssp._tenuissima	1177
 The number of bases masked by RepeatMasker:	855107
 The number of bases masked by TransposonPSI:	274429
 The total number of masked bases are:	1006851

 A.alternata_ssp._tenuissima	24350
 The number of bases masked by RepeatMasker:	318302
 The number of bases masked by TransposonPSI:	159189
 The total number of masked bases are:	424669

 A.alternata_ssp._tenuissima	635
 The number of bases masked by RepeatMasker:	849942
 The number of bases masked by TransposonPSI:	335963
 The total number of masked bases are:	1025475

 A.alternata_ssp._tenuissima	648
 The number of bases masked by RepeatMasker:	547898
 The number of bases masked by TransposonPSI:	160220
 The total number of masked bases are:	652885

 A.alternata_ssp._tenuissima	743
 The number of bases masked by RepeatMasker:	627923
 The number of bases masked by TransposonPSI:	333470
 The total number of masked bases are:	829749
```


# Gene Prediction
Gene prediction followed two steps:
Pre-gene prediction - Quality of genome assemblies were assessed using Cegma to see how many core eukaryotic genes can be identified.
Gene models were used to predict genes in the Alternaria genomes.

## Pre-gene prediction
Quality of genome assemblies was assessed by looking for the gene space in the assemblies.

```bash
  for Assembly in $(ls repeat_masked/*/*/filtered_contigs_repmask/*_contigs_unmasked.fa); do
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
    echo "$Organism"
    echo "$Strain"
  	ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/cegma
  	qsub $ProgDir/sub_cegma.sh $Assembly dna
  done
```

results were summarised:

```bash
  for File in $(ls gene_pred/cegma/A.alternata_ssp._*/*/*_dna_cegma.completeness_report); do
    basename $File
    cat $File | grep -w -e 'Complete' -e 'Partial' | grep -v 'Category' | sed -E 's/\s+/\t/g'| cut -f2,3,4
  done
```

```
  675_dna_cegma.completeness_report
  Complete	239	96.37
  Partial	241	97.18
  97.0013_dna_cegma.completeness_report
  Complete	238	95.97
  Partial	240	96.77
  97.0016_dna_cegma.completeness_report
  Complete	237	95.56
  Partial	240	96.77
  650_dna_cegma.completeness_report
  Complete	241	97.18
  Partial	244	98.39
  1082_dna_cegma.completeness_report
  Complete	241	97.18
  Partial	243	97.98
  1164_dna_cegma.completeness_report
  Complete	241	97.18
  Partial	243	97.98
  1166_dna_cegma.completeness_report
  Complete	239	96.37
  Partial	241	97.18
  1177_dna_cegma.completeness_report
  Complete	241	97.18
  Partial	243	97.98
  24350_dna_cegma.completeness_report
  Complete	241	97.18
  Partial	243	97.98
  635_dna_cegma.completeness_report
  Complete	240	96.77
  Partial	243	97.98
  648_dna_cegma.completeness_report
  Complete	240	96.77
  Partial	242	97.58
  743_dna_cegma.completeness_report
  Complete	241	97.18
  Partial	243	97.98
```

## Gene prediction 1 - Braker1 gene model training and prediction

Gene prediction was performed using Braker1.
 <!-- * The commands used to do this can be found in /gene_prediction/braker/braker_gene_pred.md -->


 Gene prediction was performed for A. alternata isolates.
 RNAseq reads were used as Hints for the location of CDS.
 A concatenated dataset of both ssp. tenuissima and ssp. gaisen RNAseq reads were used as hints for all strains.
 Genes were predicted for ssp. tenuissima using the gene model trained to ssp. tenunissima.
 Genes were predicted for ssp. gaisen using the gene model trained to ssp. gaisen.
 Genes were predicted for ssp. arborescens using the gene model trained to ssp. tenuisima.

<!--
Numbers of predicted genes were summarised, amino acid sequence extracted and
gff files extracted:

```bash
  for File in $(ls gene_pred/braker/A.alternata_*/*/*_braker/augustus.gff); do
    echo $File;
    tail -n1 $File;
    getAnnoFasta.pl $File
    OutDir=$(dirname $File)
    echo "##gff-version 3" > $OutDir/augustus_extracted.gff
    cat $File | grep -v '#' >> $OutDir/augustus_extracted.gff
  done
```

```
  gene_pred/braker/A.alternata_ssp._arborescens/675/A.alternata_ssp._arborescens_675_braker/augustus.gff
  # end gene g12046
  gene_pred/braker/A.alternata_ssp._arborescens/97.0013/A.alternata_ssp._arborescens_97.0013_braker/augustus.gff
  # end gene g11974
  gene_pred/braker/A.alternata_ssp._arborescens/97.0016/A.alternata_ssp._arborescens_97.0016_braker/augustus.gff
  # end gene g11951
  gene_pred/braker/A.alternata_ssp._gaisen/650/A.alternata_ssp._gaisen_650_braker/augustus.gff
  # end gene g12281
  gene_pred/braker/A.alternata_ssp._tenuissima/1082/A.alternata_ssp._tenuissima_1082_braker/augustus.gff
  # end gene g12346
  gene_pred/braker/A.alternata_ssp._tenuissima/1164/A.alternata_ssp._tenuissima_1164_braker/augustus.gff
  # end gene g12529
  gene_pred/braker/A.alternata_ssp._tenuissima/1166/A.alternata_ssp._tenuissima_1166_braker/augustus.gff
  # end gene g12628
  gene_pred/braker/A.alternata_ssp._tenuissima/1177/A.alternata_ssp._tenuissima_1177_braker/augustus.gff
  # end gene g12902
  gene_pred/braker/A.alternata_ssp._tenuissima/24350/A.alternata_ssp._tenuissima_24350_braker/augustus.gff
  # end gene g12074
  gene_pred/braker/A.alternata_ssp._tenuissima/635/A.alternata_ssp._tenuissima_635_braker/augustus.gff
  # end gene g12759
  gene_pred/braker/A.alternata_ssp._tenuissima/648/A.alternata_ssp._tenuissima_648_braker/augustus.gff
  # end gene g12173
  gene_pred/braker/A.alternata_ssp._tenuissima/743/A.alternata_ssp._tenuissima_743_braker/augustus.gff
  # end gene g12867
```


## Gene prediction 2 - atg.pl prediction of ORFs

Open reading frame predictions were made using the atg.pl script as part of the
path_pipe.sh pipeline. This pipeline also identifies open reading frames containing
Signal peptide sequences and RxLRs. This pipeline was run with the following commands:

```bash
  ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
  for Assembly in $(ls repeat_masked/*/*/filtered_contigs_repmask/*_contigs_unmasked.fa); do
    qsub $ProgDir/run_ORF_finder.sh $Assembly
  done
```

The Gff files from the the ORF finder are not in true Gff3 format. These were
corrected using the following commands:

```bash
  for ORF_Gff in $(ls gene_pred/ORF_finder/*/*/*_ORF.gff); do
    ProgDir=~/git_repos/emr_repos/tools/seq_tools/feature_annotation
    basename $ORF_Gff
    ORF_Gff_mod=$(echo $ORF_Gff | sed 's/ORF.gff/ORF_corrected.gff3/g')
    $ProgDir/gff_corrector.pl $ORF_Gff > $ORF_Gff_mod
  done
```
 -->


#### Aligning

 Insert sizes of the RNA seq library were unknown until a draft alignment could
 be made. To do this tophat and cufflinks were run, aligning the reads against a
 single genome. The fragment length and stdev were printed to stdout while
 cufflinks was running.

 ```bash
   for Assembly in $(ls repeat_masked/*/1166/*/*_contigs_softmasked.fa | head -n1); do
     Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
     Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
     echo "$Organism - $Strain"
     for RNADir in $(ls -d  qc_rna/paired/A.alternata_ssp._*/*); do
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

 Alignments were concatenated prior to running cufflinks:
 Cufflinks was run to produce the fragment length and stdev statistics:

 ```bash
  for Assembly in $(ls repeat_masked/*/1166/*/*_contigs_softmasked.fa | head -n1); do
    Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
    Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
    AcceptedHits=alignment/$Organism/$Strain/concatenated/concatenated.bam
    OutDir=gene_pred/cufflinks/$Organism/$Strain/concatenated_prelim
    echo "$Organism - $Strain"
    mkdir -p alignment/$Organism/$Strain/concatenated
    mkdir -p $OutDir
    samtools merge -f $AcceptedHits \
    alignment/$Organism/$Strain/1166/accepted_hits.bam \
    alignment/$Organism/$Strain/650/accepted_hits.bam
    cufflinks -o $OutDir/cufflinks -p 8 --max-intron-length 4000 $AcceptedHits 2>&1 | tee $OutDir/cufflinks/cufflinks.log
  done
 ```


 Output from stdout included:
 ```
 A.alternata_ssp._tenuissima - 1166
 You are using Cufflinks v2.2.1, which is the most recent release.
 [12:28:53] Inspecting reads and determining fragment length distribution.
 > Processed 21515 loci.                        [*************************] 100%
 > Map Properties:
 >	Normalized Map Mass: 10812067.92
 >	Raw Map Mass: 10812067.92
 >	Fragment Length Distribution: Empirical (learned)
 >	              Estimated Mean: 224.66
 >	           Estimated Std Dev: 61.11
 [12:30:13] Assembling transcripts and estimating abundances.
 > Processed 21530 loci.                        [*************************] 100%
 ```

 The Estimated Mean: 224.66 allowed calculation of of the mean insert gap to be
 -169bp 225-(197*2) where 197 was the mean read length. This was provided to tophat
 on a second run (as the -r option) along with the fragment length stdev to
 increase the accuracy of mapping.


 Then Rnaseq data was aligned to each genome assembly:

 ```bash
  for Assembly in $(ls repeat_masked/*/*/filtered_contigs_repmask/*_contigs_unmasked.fa); do
    Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
    Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
    echo "$Organism - $Strain"
    for RNADir in $(ls -d qc_rna/paired/*/*); do
      Timepoint=$(echo $RNADir | rev | cut -f1 -d '/' | rev)
      echo "$Timepoint"
      FileF=$(ls $RNADir/F/*_trim.fq.gz)
      FileR=$(ls $RNADir/R/*_trim.fq.gz)
      OutDir=alignment/$Organism/$Strain/$Timepoint
      InsertGap='-169'
      InsertStdDev='61'
      Jobs=$(qstat | grep 'tophat' | grep 'qw' | wc -l)
      while [ $Jobs -gt 1 ]; do
        sleep 10
        printf "."
        Jobs=$(qstat | grep 'tophat' | grep 'qw' | wc -l)
      done
      printf "\n"
      ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/RNAseq
      qsub $ProgDir/tophat_alignment.sh $Assembly $FileF $FileR $OutDir $InsertGap $InsertStdDev
    done
  done
 ```

 #### Braker prediction

 ```bash
  for Assembly in $(ls repeat_masked/*/*/filtered_contigs_repmask/*_contigs_unmasked.fa); do
    Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
    Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
    echo "$Organism - $Strain"
    mkdir -p alignment/$Organism/$Strain/concatenated
    samtools merge -f alignment/$Organism/$Strain/concatenated/concatenated.bam \
    alignment/$Organism/$Strain/1166/accepted_hits.bam \
    alignment/$Organism/$Strain/650/accepted_hits.bam
    OutDir=gene_pred/braker/$Organism/"$Strain"_braker
    AcceptedHits=alignment/$Organism/$Strain/concatenated/concatenated.bam
    GeneModelName="$Organism"_"$Strain"_braker
    rm -r /home/armita/prog/augustus-3.1/config/species/"$Organism"_"$Strain"_braker
    ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/braker1
    qsub $ProgDir/sub_braker_fungi.sh $Assembly $OutDir $AcceptedHits $GeneModelName
  done
 ```

 Fasta and gff files were extracted from Braker1 output.

 ```bash
 	for File in $(ls gene_pred/braker/F.*/*_braker/*/augustus.gff); do
 		getAnnoFasta.pl $File
 		OutDir=$(dirname $File)
 		echo "##gff-version 3" > $OutDir/augustus_extracted.gff
 		cat $File | grep -v '#' >> $OutDir/augustus_extracted.gff
 	done
 ```

<!--
 Cufflinks was run to compare the predicted genes to assembled transcripts:

 ```bash
 	for Assembly in $(ls repeat_masked/*/*/filtered_contigs_repmask/*_contigs_unmasked.fa); do
 		Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
 		Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
 		AcceptedHits=alignment/$Organism/$Strain/concatenated/concatenated.bam
 		OutDir=gene_pred/cufflinks/$Organism/$Strain/concatenated
 		echo "$Organism - $Strain"
 		mkdir -p $OutDir
 		samtools merge -f $AcceptedHits \
     alignment/$Organism/$Strain/1166/accepted_hits.bam \
 		alignment/$Organism/$Strain/650/accepted_hits.bam
 		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/RNAseq
 		qsub $ProgDir/sub_cufflinks.sh $AcceptedHits $OutDir/cuflfinks
 	# cufflinks -o $OutDir/cufflinks -p 8 --max-intron-length 4000 $AcceptedHits 2>&1 | tee $OutDir/cufflinks/cufflinks.log
 	done
 ```
 -->


 ## Supplimenting Braker gene models with CodingQuary genes

 Additional genes were added to Braker gene predictions, using CodingQuary in
 pathogen mode to predict additional regions.

 Fistly, aligned RNAseq data was assembled into transcripts using Cufflinks.

 Note - cufflinks doesn't always predict direction of a transcript and
 therefore features can not be restricted by strand when they are intersected.

 ```bash
  for Assembly in $(ls repeat_masked/*/*/filtered_contigs_repmask/*_contigs_unmasked.fa); do
    Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
    Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
    echo "$Organism - $Strain"
    OutDir=gene_pred/cufflinks/$Organism/$Strain/concatenated
    mkdir -p $OutDir
    AcceptedHits=alignment/$Organism/$Strain/concatenated/concatenated.bam
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/RNAseq
    qsub $ProgDir/sub_cufflinks.sh $AcceptedHits $OutDir
  done
 ```

 Secondly, genes were predicted using CodingQuary:

 ```bash
  for Assembly in $(ls repeat_masked/*/*/filtered_contigs_repmask/*_contigs_softmasked.fa); do
    Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
    Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
    echo "$Organism - $Strain"
    OutDir=gene_pred/codingquary/$Organism/$Strain
    CufflinksGTF=gene_pred/cufflinks/$Organism/$Strain/concatenated/transcripts.gtf
    ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/codingquary
    qsub $ProgDir/sub_CodingQuary.sh $Assembly $CufflinksGTF $OutDir
  done
 ```

 Then, additional transcripts were added to Braker gene models, when CodingQuary
 genes were predicted in regions of the genome, not containing Braker gene
 models:

 ```bash
  for BrakerGff in $(ls gene_pred/braker/*/*_braker/*/augustus.gff3); do
    Strain=$(echo $BrakerGff| rev | cut -d '/' -f3 | rev | sed 's/_braker//g')
    Organism=$(echo $BrakerGff | rev | cut -d '/' -f4 | rev)
    echo "$Organism - $Strain"
    Assembly=$(ls repeat_masked/$Organism/$Strain/filtered_contigs_repmask/"$Strain"_contigs_softmasked.fa)
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

    $ProgDir/add_CodingQuary_features.pl $AddGenesGff $Assembly > $FinalDir/final_genes_CodingQuary.gff3
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
    echo $DirPath;
    cat $DirPath/final_genes_Braker.pep.fasta | grep '>' | wc -l;
    cat $DirPath/final_genes_CodingQuary.pep.fasta | grep '>' | wc -l;
    cat $DirPath/final_genes_combined.pep.fasta | grep '>' | wc -l;
    echo "";
  done
 ```

```
  gene_pred/final/A.alternata_ssp._arborescens/675/final
  12196
  511
  12707

  gene_pred/final/A.alternata_ssp._arborescens/97.0013/final
  12310
  457
  12767

  gene_pred/final/A.alternata_ssp._arborescens/97.0016/final
  12178
  502
  12680

  gene_pred/final/A.alternata_ssp._gaisen/650/final
  12485
  667
  13152

  gene_pred/final/A.alternata_ssp._tenuissima/1082/final
  12514
  606
  13120

  gene_pred/final/A.alternata_ssp._tenuissima/1164/final
  12631
  643
  13274

  gene_pred/final/A.alternata_ssp._tenuissima/1166/final
  13048
  686
  13734

  gene_pred/final/A.alternata_ssp._tenuissima/1177/final
  13170
  712
  13882

  gene_pred/final/A.alternata_ssp._tenuissima/24350/final
  12289
  622
  12911

  gene_pred/final/A.alternata_ssp._tenuissima/635/final
  13163
  757
  13920

  gene_pred/final/A.alternata_ssp._tenuissima/648/final
  12530
  504
  13034

  gene_pred/final/A.alternata_ssp._tenuissima/743/final
  13330
  686
  14016
```

#Functional annotation

Interproscan was used to give gene models functional annotations.

```bash
  for Proteome in $(ls gene_pred/final/A.alternata_ssp._*/*/*/final_genes_combined.pep.fasta); do
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/interproscan
    $ProgDir/sub_interproscan.sh $Proteome
  done
```

Following interproscan annotation split files were combined using the following commands:

```bash
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/interproscan
  for StrainPath in $(ls -d gene_pred/interproscan/*/*); do
  Strain=$(basename $StrainPath)
  Organism=$(echo $StrainPath | rev | cut -d "/" -f2 | rev)
  echo $Strain
  PredGenes=$(ls gene_pred/final/"$Organism"/"$Strain"/*/final_genes_combined.pep.fasta)
  InterProRaw=gene_pred/interproscan/"$Organism"/"$Strain"/raw
  $ProgDir/append_interpro.sh $PredGenes $InterProRaw
  done
```

## Small secreted proteins

Putative RxLR genes were identified within Augustus gene models using a number
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
  for Proteome in $(ls gene_pred/final/A.alternata_ssp._*/*/*/final_genes_combined.pep.fasta); do
    SplitfileDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/signal_peptides
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/signal_peptides
    Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
    SplitDir=gene_pred/braker_split/$Organism/$Strain
    mkdir -p $SplitDir
    BaseName="$Organism""_$Strain"_braker_preds
    $SplitfileDir/splitfile_500.py --inp_fasta $Proteome --out_dir $SplitDir --out_base $BaseName
    for File in $(ls $SplitDir/*_braker_preds_*); do
    Jobs=$(qstat | grep 'pred_sigP' | grep 'qw' | wc -l)
    while [ $Jobs -gt '1' ]; do
    sleep 10
    printf "."
    Jobs=$(qstat | grep 'pred_sigP' | grep 'qw' | wc -l)
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
  for Proteome in $(ls gene_pred/final/*/*/*/final_genes_combined.pep.fasta); do
    Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/transmembrane_helices
    qsub $ProgDir/submit_TMHMM.sh $Proteome
  done
```


B) SwissProt

```bash
	for Proteome in $(ls gene_pred/final/A.alternata_ssp._*/*/*/final_genes_combined.pep.fasta); do
		Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
		Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
		OutDir=gene_pred/swissprot/$Organism/$Strain
		SwissDbDir=../../uniprot/swissprot
		SwissDbName=uniprot_sprot
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/swissprot
		qsub $ProgDir/sub_swissprot.sh $Proteome $OutDir $SwissDbDir $SwissDbName
	done
```

```bash
	for SwissTable in $(ls gene_pred/swissprot/*/*/swissprot_v2015_10_hits.tbl); do
		Strain=$(echo $SwissTable | rev | cut -f2 -d '/' | rev)
		Organism=$(echo $SwissTable | rev | cut -f3 -d '/' | rev)
		echo "$Organism - $Strain"
		OutTable=gene_pred/swissprot/$Organism/$Strain/swissprot_v2015_tophit_parsed.tbl
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/swissprot
		$ProgDir/swissprot_parser.py --blast_tbl $SwissTable --blast_db_fasta ../../uniprot/swissprot/uniprot_sprot.fasta > $OutTable
	done
```


### C) From Augustus gene models - Effector identification using EffectorP

Required programs:
 * EffectorP.py

```bash
  for Proteome in $(ls gene_pred/final/A.alternata_ssp._*/*/*/final_genes_combined.pep.fasta); do
    Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
    BaseName="$Organism"_"$Strain"_EffectorP
    OutDir=analysis/effectorP/$Organism/$Strain
    ProgDir=~/git_repos/emr_repos/tools/seq_tools/feature_annotation/fungal_effectors
    qsub $ProgDir/pred_effectorP.sh $Proteome $BaseName $OutDir
  done
```


# Genomic analysis
The first analysis was based upon BLAST searches for genes known to be involved in toxin production


## Genes with homology to PHIbase
Predicted gene models were searched against the PHIbase database using tBLASTx.

```bash
  for Subject in $(ls repeat_masked/A.alternata_ssp._*/*/filtered_contigs_repmask/*_contigs_unmasked.fa); do
    ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
    Query=../../phibase/v3.8/PHI_accessions.fa
    qsub $ProgDir/blast_pipe.sh $Query protein $Subject
  done
```


## Presence and genes with homology to Alternaria toxins

The first analysis was based upon BLAST searches for genes known to be involved in toxin production



```bash
  for Subject in $(ls repeat_masked/A.alternata_ssp._*/*/filtered_contigs_repmask/*_contigs_unmasked.fa); do
    ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
    Query=analysis/blast_homology/CDC_genes/A.alternata_CDC_genes.fa
    qsub $ProgDir/blast_pipe.sh $Query dna $Subject
  done
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
    NumHits=1
    $ProgDir/blast2gff.pl $Column2 $NumHits $BlastHits > $HitsGff
  done
```


Extracted gff files of the BLAST hit locations were intersected with gene models
to identify if genes were predicted for these homologs:

```bash
for HitsGff in $(ls analysis/blast_homology/*/*/*_A.alternata_CDC_genes.fa_homologs.gff); do
Strain=$(echo $HitsGff | rev | cut -f2 -d '/' | rev)
Organism=$(echo $HitsGff | rev | cut -f3 -d '/' | rev)
Proteins=$(ls gene_pred/braker/$Organism/$Strain/*/augustus.gff)
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
