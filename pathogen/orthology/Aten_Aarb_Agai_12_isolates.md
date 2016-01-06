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


### for A.alt_ssp.ten_apple_path 635
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


<!--
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
-->


# 6) Downstream analysis

Particular orthogroups were analysed for expansion in isolates.

This section details the commands used and the results observed.


### 6.1 ) Clade unique gene families

The genes unique to A. tenuissima apple pathotypes were identified within the orthology analysis.

First variables were set:
```bash
  WorkDir=analysis/orthology/orthomcl/At_Aa_Ag_all_isolates
  Orthogroups=$WorkDir/At_Aa_Ag_all_isolates_orthogroups.txt
  GoodProts=$WorkDir/goodProteins/goodProteins.fasta
  Braker_genes=
```

#### 6.1.a ) Orthologroups only containing A. tenuissima genes were extracted:

```bash
for num in 1; do
  # AtenUniq
  AtenUniqDir=$WorkDir/A.tenuissima_unique
  Uniq_Aten_groups=$AtenPathUniqDir/A.tenuissima_unique.txt
  mkdir -p $AtenUniqDir
  cat $Orthogroups | grep -e 'At_1' | grep -e 'At_2' | grep -e 'At_3' | grep -e 'At_4' | grep -e 'At_5' | grep -e 'At_6' | grep -e 'At_7' | grep -e 'At_8' | grep -v -e 'Aa' -e 'Ag' > $Uniq_Aten_groups
  echo "The number of orthogroups unique to A.tenuissima apple pathotype are:"
  cat $Uniq_Aten_groups | wc -l
  echo "The following number genes from isolate 648 are contained in these orthogorups:"
  cat $Uniq_Aten_groups | grep -o -e 'At_1' | wc -l
  echo "The following number genes from isolate 1082 are contained in these orthogorups:"
  cat $Uniq_Aten_groups | grep -o -e 'At_2' | wc -l
  echo "The following number genes from isolate 1164 are contained in these orthogorups:"
  cat $Uniq_Aten_groups | grep -o -e 'At_3' | wc -l
  echo "The following number genes from isolate 24350 are contained in these orthogorups:"
  cat $Uniq_Aten_groups | grep -o -e 'At_4' | wc -l
  echo "The following number genes from isolate 648 are contained in these orthogorups:"
  cat $Uniq_Aten_groups | grep -o -e 'At_5' | wc -l
  echo "The following number genes from isolate 743 are contained in these orthogorups:"
  cat $Uniq_Aten_groups | grep -o -e 'At_6' | wc -l
  echo "The following number genes from isolate 1166 are contained in these orthogorups:"
  cat $Uniq_Aten_groups | grep -o -e 'At_7' | wc -l
  echo "The following number genes from isolate 1177 are contained in these orthogorups:"
  cat $Uniq_Aten_groups | grep -o -e 'At_8' | wc -l
done
```

```
The number of orthogroups unique to A.tenuissima apple pathotype are:
21
The following number genes from isolate 648 are contained in these orthogorups:
22
The following number genes from isolate 1082 are contained in these orthogorups:
25
The following number genes from isolate 1164 are contained in these orthogorups:
24
The following number genes from isolate 24350 are contained in these orthogorups:
23
The following number genes from isolate 648 are contained in these orthogorups:
23
The following number genes from isolate 743 are contained in these orthogorups:
25
The following number genes from isolate 1166 are contained in these orthogorups:
22
The following number genes from isolate 1177 are contained in these orthogorups:
22
```

#### 6.1.b ) Orthologroups only containing A. arborescens genes were extracted:

```bash
for num in 1; do
  # AarbUniq
  AarbUniqDir=$WorkDir/A.arborescens_unique
  Uniq_Aarb_groups=$AarbUniqDir/A.arborescens_unique.txt
  mkdir -p $AarbUniqDir
  cat $Orthogroups | grep -e 'Aa_1' | grep -e 'Aa_2' | grep -e 'Aa_3' | grep -v -e 'At' -e 'Ag' > $Uniq_Aarb_groups
  echo "The number of orthogroups unique to A.arborescens are:"
  cat $Uniq_Aarb_groups | wc -l
  echo "The following number genes from isolate 675 are contained in these orthogorups:"
  cat $Uniq_Aarb_groups | grep -o -e 'Aa_1' | wc -l
  echo "The following number genes from isolate 97.0013 are contained in these orthogorups:"
  cat $Uniq_Aarb_groups | grep -o -e 'Aa_2' | wc -l
  echo "The following number genes from isolate 97.0016 are contained in these orthogorups:"
  cat $Uniq_Aarb_groups | grep -o -e 'Aa_3' | wc -l
done
```

```
The number of orthogroups unique to A.arborescens are:
41
The following number genes from isolate 675 are contained in these orthogorups:
45
The following number genes from isolate 97.0013 are contained in these orthogorups:
45
The following number genes from isolate 97.0016 are contained in these orthogorups:
44
```

#### 6.1.c ) Orthologroups only containing A. gaisen genes were extracted:

```bash
for num in 1; do
  # AgaiPathUniq
  AgaiUniqDir=$WorkDir/A.gaisen_unique
  Uniq_Agai_groups=$AgaiUniqDir/A.gaisen_unique.txt
  mkdir -p $AgaiUniqDir
  cat $Orthogroups | grep -e 'Ag_1' | grep -v  -e 'At' -e 'Aa' > $Uniq_Agai_groups
  echo "The number of orthogroups unique to A.gaisen pear pathotype pathotype are:"
  cat $Uniq_Agai_groups | wc -l
  echo "The following number genes from isolate 650 are contained in these orthogorups:"
  cat $Uniq_Agai_groups | grep -o -e 'Ag_1' | wc -l
done
```

```
The number of orthogroups unique to A.gaisen pear pathotype pathotype are:
7
The following number genes from isolate 650 are contained in these orthogorups:
15
```

#### 6.1.d ) Orthologroups only containing A. tenuissima non pathotype genes were extracted:

```bash
for num in 1; do
  # AtenNonPathUniq
  AtenNonPathUniqDir=$WorkDir/A.tenuissima_non_pathotype_unique
  Uniq_AtenNonPath_groups=$AtenPathUniqDir/A.tenuissima_non_pathotype_unique.txt
  mkdir -p $AtenNonPathUniqDir
  cat $Orthogroups | grep -e 'At_1' | grep -e 'At_2' | grep -e 'At_3' | grep -e 'At_4' | grep -v -e 'At_5' -e 'At_6' -e 'At_7' -e 'At_8' -e 'Aa' -e 'Ag' > $Uniq_AtenNonPath_groups
  echo "The number of orthogroups unique to A.tenuissima apple pathotype are:"
  cat $Uniq_AtenNonPath_groups | wc -l
  echo "The following number genes from isolate 648 are contained in these orthogorups:"
  cat $Uniq_AtenNonPath_groups | grep -o -e 'At_1' | wc -l
  echo "The following number genes from isolate 1082 are contained in these orthogorups:"
  cat $Uniq_AtenNonPath_groups | grep -o -e 'At_2' | wc -l
  echo "The following number genes from isolate 1164 are contained in these orthogorups:"
  cat $Uniq_AtenNonPath_groups | grep -o -e 'At_3' | wc -l
  echo "The following number genes from isolate 24350 are contained in these orthogorups:"
  cat $Uniq_AtenNonPath_groups | grep -o -e 'At_4' | wc -l
done
```

```
The number of orthogroups unique to A.tenuissima apple pathotype are:
0
The following number genes from isolate 648 are contained in these orthogorups:
0
The following number genes from isolate 1082 are contained in these orthogorups:
0
The following number genes from isolate 1164 are contained in these orthogorups:
0
The following number genes from isolate 24350 are contained in these orthogorups:
0
```

#### 6.1.d ) Orthologroups only containing A. tenuissima apple pathotype genes were extracted:

```bash
for num in 1; do
  # AtenPathUniq
  AtenPathUniqDir=$WorkDir/A.tenuissima_apple_pathotype_unique
  Uniq_AtenPath_groups=$AtenPathUniqDir/A.tenuissima_apple_pathotype_unique.txt
  mkdir -p $AtenPathUniqDir
  cat $Orthogroups | grep -e 'At_5' | grep -e 'At_6' | grep -e 'At_7' | grep -e 'At_8' | grep -v -e 'At_1' -e 'At_2' -e 'At_3' -e 'At_4' -e 'Aa' -e 'Ag' > $Uniq_AtenPath_groups
  echo "The number of orthogroups unique to A.tenuissima apple pathotype are:"
  cat $Uniq_AtenPath_groups | wc -l
  echo "The following number genes from isolate 635 are contained in these orthogorups:"
  cat $Uniq_AtenPath_groups | grep -o -e 'At_5' | wc -l
  echo "The following number genes from isolate 743 are contained in these orthogorups:"
  cat $Uniq_AtenPath_groups | grep -o -e 'At_6' | wc -l
  echo "The following number genes from isolate 1166 are contained in these orthogorups:"
  cat $Uniq_AtenPath_groups | grep -o -e 'At_7' | wc -l
  echo "The following number genes from isolate 1177 are contained in these orthogorups:"
  cat $Uniq_AtenPath_groups | grep -o -e 'At_8' | wc -l
done
```

```
The number of orthogroups unique to A.tenuissima apple pathotype are:
9
The following number genes from isolate 635 are contained in these orthogorups:
10
The following number genes from isolate 743 are contained in these orthogorups:
10
The following number genes from isolate 1166 are contained in these orthogorups:
10
The following number genes from isolate 1177 are contained in these orthogorups:
10
```


#### 6.1.e ) Orthologroups only containing genes from apple and pear pathotypes were extracted:

```bash
for num in 1; do
  # PathUniq
  PathUniqDir=$WorkDir/Pathotype_unique
  Uniq_Path_groups=$PathUniqDir/Path_unique.txt
  mkdir -p $PathUniqDir
  cat $Orthogroups | grep -e 'Ag' | grep -e 'At_5' | grep -e 'At_6' | grep -e 'At_7' | grep -e 'At_8' | grep -v -e 'At_1' -e 'At_2' -e 'At_3' -e 'At_4' -e 'Aa' > $Uniq_Path_groups
  echo "The number of orthogroups unique to apple and pear pathotypes are:"
  cat $Uniq_Path_groups | wc -l
  echo "The following number genes from isolate 650 are contained in these orthogorups:"
  cat $Uniq_Path_groups | grep -o -e 'Ag_1' | wc -l
  echo "The following number genes from isolate 648 are contained in these orthogorups:"
  cat $Uniq_Path_groups | grep -o -e 'At_5' | wc -l
  echo "The following number genes from isolate 743 are contained in these orthogorups:"
  cat $Uniq_Path_groups | grep -o -e 'At_6' | wc -l
  echo "The following number genes from isolate 1166 are contained in these orthogorups:"
  cat $Uniq_Path_groups | grep -o -e 'At_7' | wc -l
  echo "The following number genes from isolate 1177 are contained in these orthogorups:"
  cat $Uniq_Path_groups | grep -o -e 'At_8' | wc -l
done
```

```
The number of orthogroups unique to apple and pear pathotypes are:
18
The following number genes from isolate 650 are contained in these orthogorups:
20
The following number genes from isolate 648 are contained in these orthogorups:
25
The following number genes from isolate 743 are contained in these orthogorups:
23
The following number genes from isolate 1166 are contained in these orthogorups:
18
The following number genes from isolate 1177 are contained in these orthogorups:
25
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
```
