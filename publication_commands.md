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
Gene models were used to predict genes in the Neonectria genome. This used results from CEGMA as hints for gene models.

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
 * The commands used to do this can be found in /gene_prediction/braker/braker_gene_pred.md



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


#Functional annotation

Interproscan was used to give gene models functional annotations.

```bash
  for Proteome in $(ls gene_pred/braker/A.alternata_ssp._*/*/*/augustus.aa); do
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/interproscan
    $ProgDir/sub_interproscan.sh $Proteome
  done
```

```bash
  screen -a
  for Proteome in $(ls gene_pred/braker/A.alternata_ssp._*/*/*/augustus.aa); do
    InFile=$(basename $Proteome)
    echo $InFile
    Organism=$(echo $Proteome  | rev | cut -d "/" -f3 | rev)
    Strain=$(echo $Proteome  | rev | cut -d "/" -f2 | rev)
    Split_Dir=gene_pred/interproscan/$Organism/$Strain
    mkdir -p $Split_Dir
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/signal_peptides
    $ProgDir/splitfile_500.py --inp_fasta $Proteome  --out_dir $Split_Dir --out_base "$InFile"_split
    for File in $(ls $Split_Dir/*_split*); do
      Jobs=$(qstat | grep 'run_interp' | grep 'qw' | wc -l)
      while [ $Jobs -gt 1 ]; do
        sleep 10
        printf "."
        Jobs=$(qstat | grep 'run_interp' | grep 'qw' | wc -l)
      done
      printf "\n"
      echo $File
      ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/interproscan
      qsub $ProgDir/run_interproscan.sh $File
    done
  done
```

Following interproscan annotation split files were combined using the following commands:

```bash
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/interproscan
  for StrainPath in $(ls -d gene_pred/interproscan/*/*); do
  Strain=$(basename $StrainPath)
  Organism=$(echo $StrainPath | rev | cut -d "/" -f2 | rev)
  echo $Strain
  PredGenes=gene_pred/augustus/"$Organism"/"$Strain"/"$Strain"_augustus_preds.aa
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


### A) From Augustus gene models - Signal peptide & RxLR motif

Required programs:
 * SigP
 * biopython


Proteins that were predicted to contain signal peptides were identified using
the following commands:

```bash
  for Proteome in $(ls gene_pred/braker/A.alternata_ssp._*/*/*/augustus.aa); do
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


B) SwissProt

```bash
  screen -a
  qlogin
  ProjDir=/home/groups/harrisonlab/project_files/alternaria
  cd $ProjDir
  for Proteome in $(ls gene_pred/braker/A.alternata_ssp._*/*/*/augustus.aa); do
    Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
    OutDir=$ProjDir/gene_pred/augustus/$Species/$Strain/swissplot
    mkdir -p $OutDir
    blastp \
    -db /home/groups/harrisonlab/uniprot/swissprot/uniprot_sprot \
    -query $ProjDir/$Proteome \
    -out $OutDir/swissprot_v2015_10_hits.tbl \
    -evalue 1e-100 \
    -outfmt 6 \
    -num_threads 16 \
    -num_alignments 10
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
