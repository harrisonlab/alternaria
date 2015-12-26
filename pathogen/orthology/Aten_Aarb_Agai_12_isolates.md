# Orthology analysis between 12 Alternaria spp. isolates
 Including 1 A.gaisen (pear pathotype), 3 A. arborescens isolates,
 8 A. tenuissima isolates (4 apple pathotype).


```bash
  ProjDir=/home/groups/harrisonlab/project_files/alternaria
  cd $ProjDir
  IsolateAbrv=At_Aa_Ag_all_isolates
  WorkDir=analysis/orthology/orthomcl/$IsolateAbrv
  mkdir -p $WorkDir
  mkdir -p $WorkDir/formatted
  mkdir -p $WorkDir/goodProteins
  mkdir -p $WorkDir/badProteins  
```

## 1.a) Format fasta files - A. tenuissima


### for A. tenuissima isolate 648
```bash
  Taxon_code=At_1
  Fasta_file=$(ls gene_pred/braker/A.alternata_ssp._tenuissima/648/*/*.aa)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```

### for A.alt_ssp.ten 1082
```bash
  Taxon_code=At_2
  Fasta_file=$(ls gene_pred/braker/A.alternata_ssp._tenuissima/1082/*/*.aa)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```

### for A.alt_ssp.ten 1164
```bash
  Taxon_code=At_3
  Fasta_file=$(ls gene_pred/braker/A.alternata_ssp._tenuissima/1164/*/*.aa)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```

### for A.alt_ssp.ten 24350
```bash
  Taxon_code=At_4
  Fasta_file=$(ls gene_pred/braker/A.alternata_ssp._tenuissima/24350/*/*.aa)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```

## 1.b) Format fasta files - A. tenuissima apple pathotype


### for A.alt_ssp.ten_apple_path 648
```bash
  Taxon_code=At_5
  Fasta_file=$(ls gene_pred/braker/A.alternata_ssp._tenuissima/635/*/*.aa)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```

### for A.alt_ssp.ten_apple_path 743
```bash
  Taxon_code=At_6
  Fasta_file=$(ls gene_pred/braker/A.alternata_ssp._tenuissima/743/*/*.aa)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```

### for A.alt_ssp.ten_apple_path 1166
```bash
  Taxon_code=At_7
  Fasta_file=$(ls gene_pred/braker/A.alternata_ssp._tenuissima/1166/*/*.aa)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```

### for A.alt_ssp.ten_apple_path 1177
```bash
  Taxon_code=At_8
  Fasta_file=$(ls gene_pred/braker/A.alternata_ssp._tenuissima/1177/*/*.aa)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```

## 1.c) Format fasta files - A. arborescens

### for A.alt_ssp.arb 675
```bash
  Taxon_code=Aa_1
  Fasta_file=$(ls gene_pred/braker/A.alternata_ssp._arborescens/675/*/*.aa)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```

### for A.alt_ssp.arb 97.0013
```bash
  Taxon_code=Aa_2
  Fasta_file=$(ls gene_pred/braker/A.alternata_ssp._arborescens/97.0013/*/*.aa)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```

### for A.alt_ssp.arb 97.0016
```bash
  Taxon_code=Aa_3
  Fasta_file=$(ls gene_pred/braker/A.alternata_ssp._arborescens/97.0016/*/*.aa)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```

## 1.d) Format fasta files - A. gaisen

### for A.alt_ssp.gai 650
```bash
  Taxon_code=Ag_1
  Fasta_file=$(ls gene_pred/braker/A.alternata_ssp._gaisen/650/*/*.aa)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```



## 2) Filter proteins into good and poor sets.

```bash
  Input_dir=$WorkDir/formatted
  Min_length=10
  Max_percent_stops=20
  Good_proteins_file=$WorkDir/goodProteins/goodProteins.fasta
  Poor_proteins_file=$WorkDir/badProteins/poorProteins.fasta
  orthomclFilterFasta $Input_dir $Min_length $Max_percent_stops $Good_proteins_file $Poor_proteins_file
```

## 3.1) Perform an all-vs-all blast of the proteins

```bash
BlastDB=$WorkDir/blastall/$IsolateAbrv.db

makeblastdb -in $Good_proteins_file -dbtype prot -out $BlastDB
BlastOut=$WorkDir/all-vs-all_results.tsv
mkdir -p $WorkDir/splitfiles

SplitDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/signal_peptides
$SplitDir/splitfile_500.py --inp_fasta $Good_proteins_file --out_dir $WorkDir/splitfiles --out_base goodProteins

ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/orthology  
for File in $(find $WorkDir/splitfiles); do
Jobs=$(qstat | grep 'blast_500' | grep 'qw' | wc -l)
while [ $Jobs -gt 1 ]; do
sleep 10
printf "."
Jobs=$(qstat | grep 'blast_500' | grep 'qw' | wc -l)
done
printf "\n"
echo $File
BlastOut=$(echo $File | sed 's/.fa/.tab/g')
qsub $ProgDir/blast_500.sh $BlastDB $File $BlastOut
done
```

<- progress
<!-- 
## 3.2) Merge the all-vs-all blast results  
```bash  
  MergeHits="$IsolateAbrv"_blast.tab
  printf "" > $MergeHits
  for Num in $(ls $WorkDir/splitfiles/*.tab | rev | cut -f1 -d '_' | rev | sort -n); do
    File=$(ls $WorkDir/splitfiles/*_$Num)
    cat $File
  done > $MergeHits
```

## 4) Perform ortholog identification

```bash
  ProgDir=~/git_repos/emr_repos/tools/pathogen/orthology/orthoMCL
  MergeHits="$IsolateAbrv"_blast.tab
  GoodProts=$WorkDir/goodProteins/goodProteins.fasta
  qsub $ProgDir/qsub_orthomcl.sh $MergeHits $GoodProts
```

## 5) Plot venn diagrams:

```bash
  ProgDir=~/git_repos/emr_repos/tools/pathogen/orthology/venn_diagrams
  $ProgDir/ven_diag_5_way.R --inp $WorkDir/"$IsolateAbrv"_orthogroups.tab --out $WorkDir/"$IsolateAbrv"_orthogroups.pdf
```

Output was a pdf file of the venn diagram.

The following additional information was also provided. The format of the
following lines is as follows:

Isolate name (total number of orthogroups)
number of unique singleton genes
number of unique groups of inparalogs

```
  [1] "Pcac (8814)"
  [1] 567
  [1] 118
  [1] "Pcap (7646)"
  [1] 333
  [1] 52
  [1] "Pinf (8335)"
  [1] 562
  [1] 100
  [1] "Ppar (8987)"
  [1] 695
  [1] 79
  [1] "Psoj (9156)"
  [1] 883
  [1] 388
```

# 6) Downstream analysis

Particular orthogroups were analysed for expansion in isolates.

This section details the commands used and the results observed.


### 6.1 ) P. cactotum unique gene families

The genes unique to P.cactorum were identified within the orthology analysis.

First variables were set:
```bash
  WorkDir=analysis/orthology/orthomcl/Pcac_Pinf_Ppar_Pcap_Psoj
  PcacUniqDir=$WorkDir/Pcac_unique
  Orthogroups=$WorkDir/Pcac_Pinf_Ppar_Pcap_Psoj_orthogroups.txt
  GoodProts=$WorkDir/goodProteins/goodProteins.fasta
  Braker_genes=gene_pred/braker/P.cactorum/10300/P.cactorum/augustus.aa
  Uniq_Pcac_groups=$PcacUniqDir/Pcac_uniq_orthogroups.txt
  mkdir -p $PcacUniqDir
```

Orthologroups only containing P.cactorum 10300 genes were extracted:

```bash
  cat $Orthogroups | grep -v 'Pinf' | grep -v 'Ppar' | grep -v 'Pcap' | grep -v 'Psoj' > $Uniq_Pcac_groups
  echo "The number of orthogroups unique to P. cactorum are:"
  cat $Uniq_Pcac_groups | wc -l
  echo "The following number genes are contained in these orthogorups:"
  cat $Uniq_Pcac_groups | grep -o 'Pcac|' | wc -l  
```

```
  The number of orthogroups unique to P. cactorum are:
  118
  The following number genes are contained in these orthogorups:
  328
```

### 6.2.a) P. cactorum unique RxLR families

P. cactorum strain 10300 RxLR genes were parsed to the same format as the gene
names used in the analysis:

```bash
  RxLR_Names_10300=analysis/RxLR_effectors/combined_evidence/P.cactorum/10300/10300_Aug_RxLR_EER_motif_hmm_headers.txt
  WorkDir=analysis/orthology/orthomcl/Pcac_Pinf_Ppar_Pcap_Psoj
  RxLR_Dir=$WorkDir/Pcac_RxLR
  Orthogroups=$WorkDir/Pcac_Pinf_Ppar_Pcap_Psoj_orthogroups.txt
  RxLR_ID_10300=$RxLR_Dir/10300_aug_RxLR_EER_IDs.txt
  mkdir -p $RxLR_Dir
  cat $RxLR_Names_10300 | sed 's/g/Pcac|g/g' > $RxLR_ID_10300
```

Ortholog groups containing RxLR proteins were identified using the following
commands:
```bash
  echo "The number of RxLRs searched for is:"
  cat $RxLR_ID_10300 | wc -l
  echo "Of these, the following number were found in orthogroups:"
  RxLR_Orthogroup_hits_10300=$RxLR_Dir/Pcac_RxLR_Orthogroups_hits.txt
  cat $Orthogroups | grep -o -w -f $RxLR_ID_10300 > $RxLR_Orthogroup_hits_10300
  cat $RxLR_Orthogroup_hits_10300 | wc -l
  echo "These were distributed through the following number of Orthogroups:"
  RxLR_Orthogroup_10300=$RxLR_Dir/Pcac_RxLR_Orthogroups.txt
  cat $Orthogroups | grep -w -f $RxLR_ID_10300 > $RxLR_Orthogroup_10300
  cat $RxLR_Orthogroup_10300 | wc -l
  echo "The following RxLRs were found in Pcac unique orthogroups:"
  RxLR_Pcac_uniq_groups=$RxLR_Dir/Pcac_uniq_RxLR_Orthogroups_hits.txt
  cat $RxLR_Orthogroup_10300 | grep -v 'Pinf' | grep -v 'Ppar' | grep -v 'Pcap' | grep -v 'Psoj' > $RxLR_Pcac_uniq_groups
  cat $RxLR_Pcac_uniq_groups | wc -l
  echo "The following RxLRs were found in Group1 unique orthogroups:"
  RxLR_Group1_uniq_groups=$RxLR_Dir/Group1_RxLR_Orthogroups_hits.txt
  cat $RxLR_Orthogroup_10300 | grep -v 'Pcap' | grep -v 'Psoj' > $RxLR_Group1_uniq_groups
  cat $RxLR_Group1_uniq_groups | wc -l
``` -->
