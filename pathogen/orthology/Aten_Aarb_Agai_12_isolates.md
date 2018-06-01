# Orthology analysis between 12 Alternaria spp. isolates
 Including 1 A.gaisen (pear pathotype), 3 A. arborescens isolates,
 8 A. tenuissima isolates (4 apple pathotype).


```bash
  ProjDir=/data/scratch/armita/alternaria
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
  Fasta_file=$(ls /home/groups/harrisonlab/project_files/alternaria/gene_pred/final/A.alternata_ssp._tenuissima/648/*/final_genes_appended_renamed.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```

### for A.alt_ssp.ten 1082
```bash
  Taxon_code=At_2
  Fasta_file=$(ls /home/groups/harrisonlab/project_files/alternaria/gene_pred/final/A.alternata_ssp._tenuissima/1082/*/final_genes_appended_renamed.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```

### for A.alt_ssp.ten 1164
```bash
  Taxon_code=At_3
  Fasta_file=$(ls /home/groups/harrisonlab/project_files/alternaria/gene_pred/final/A.alternata_ssp._tenuissima/1164/*/final_genes_appended_renamed.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```

### for A.alt_ssp.ten 24350
```bash
  Taxon_code=At_4
  Fasta_file=$(ls /home/groups/harrisonlab/project_files/alternaria/gene_pred/final/A.alternata_ssp._tenuissima/24350/*/final_genes_appended_renamed.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```

## 1.b) Format fasta files - A. tenuissima apple pathotype


### for A.alt_ssp.ten_apple_path 635
```bash
  Taxon_code=At_5
  Fasta_file=$(ls /home/groups/harrisonlab/project_files/alternaria/gene_pred/final/A.alternata_ssp._tenuissima/635/*/final_genes_appended_renamed.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```

### for A.alt_ssp.ten_apple_path 743
```bash
  Taxon_code=At_6
  Fasta_file=$(ls /home/groups/harrisonlab/project_files/alternaria/gene_pred/final/A.alternata_ssp._tenuissima/743/*/final_genes_appended_renamed.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```

### for A.alt_ssp.ten_apple_path 1166
```bash
  Taxon_code=At_7
  Fasta_file=$(ls gene_pred/final/A.alternata_ssp_tenuissima/1166/final/final_genes_appended_renamed.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```

### for A.alt_ssp.ten_apple_path 1177
```bash
  Taxon_code=At_8
  Fasta_file=$(ls /home/groups/harrisonlab/project_files/alternaria/gene_pred/final/A.alternata_ssp._tenuissima/1177/*/final_genes_appended_renamed.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```

## 1.c) Format fasta files - A. arborescens

### for A.alt_ssp.arb 675
```bash
  Taxon_code=Aa_1
  Fasta_file=$(ls /home/groups/harrisonlab/project_files/alternaria/gene_pred/final/A.alternata_ssp._arborescens/675/*/final_genes_appended_renamed.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```

### for A.alt_ssp.arb 97.0013
```bash
  Taxon_code=Aa_2
  Fasta_file=$(ls /home/groups/harrisonlab/project_files/alternaria/gene_pred/final/A.alternata_ssp._arborescens/97.0013/*/final_genes_appended_renamed.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```

### for A.alt_ssp.arb 97.0016
```bash
  Taxon_code=Aa_3
  Fasta_file=$(ls /home/groups/harrisonlab/project_files/alternaria/gene_pred/final/A.alternata_ssp._arborescens/97.0016/*/final_genes_appended_renamed.pep.fasta)
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```

## 1.d) Format fasta files - A. gaisen

### for A.alt_ssp.gai 650
```bash
  Taxon_code=Ag_1
  Fasta_file=$(ls gene_pred/final/A.gaisen/650/final/final_genes_appended_renamed.pep.fasta)
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
  qsub $ProgDir/qsub_orthomcl.sh $MergeHits $GoodProts 5
```


Also try using orthofinder

```bash
qlogin -pe smp 16

#16 threads used
ProjDir=/data/scratch/armita/alternaria
cd $ProjDir
IsolateAbrv=At_Aa_Ag_all_isolates
WorkDir=analysis/orthology/orthomcl/$IsolateAbrv
orthofinder -f $WorkDir/formatted -t 3 -a 3
```

orthofinder results:

```
OrthoFinder assigned 157038 genes (99.2% of total) to 14187 orthogroups. Fifty percent of all genes were in orthogroups
with 12 or more genes (G50 was 12) and were contained in the largest 6160 orthogroups (O50 was 6160). There were 10669
orthogroups with all species present and 10016 of these consisted entirely of single-copy genes.
```

output files are in:
```bash
ls $WorkDir/formatted/Results_Apr10
```

## 5) Plot venn diagrams:

Orthofinder output:

```bash
  GoodProts=$WorkDir/goodProteins/goodProteins.fasta
  ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/orthology/orthoMCL
  $ProgDir/orthoMCLgroups2tab.py $GoodProts $WorkDir/formatted/Results_Apr10/Orthogroups.txt > $WorkDir/formatted/Results_Apr10/"$IsolateAbrv"_orthogroups.tab
  ProgDir=~/git_repos/emr_repos/scripts/phytophthora/pathogen/orthology
  $ProgDir/Pcac_Pinf_Ppar_Pcap_Psoj_venn.r --inp $WorkDir/formatted/Results_Apr10/"$IsolateAbrv"_orthogroups.tab --out $WorkDir/formatted/Results_Apr10/"$IsolateAbrv"_orthogroups.pdf
```

Orthomcl output:
```bash
  ProgDir=~/git_repos/emr_repos/scripts/phytophthora/pathogen/orthology
  $ProgDir/Pcac_Pinf_Ppar_Pcap_Psoj_venn.r --inp $WorkDir/"$IsolateAbrv"_orthogroups.tab --out $WorkDir/"$IsolateAbrv"_orthogroups.pdf
```

Output was a pdf file of the venn diagram.

The following additional information was also provided. The format of the
following lines is as follows:

Isolate name
total number of orthogroups
number of unique singleton genes
number of unique groups of inparalogs

```
```

# 6) Downstream analysis

Particular orthogroups were analysed for expansion in isolates.

This section details the commands used and the results observed.


## Extracting fasta files for orthogroups:

```bash
  ProgDir=~/git_repos/emr_repos/tools/pathogen/orthology/orthoMCL
  OrthogroupsTxt=$(ls $WorkDir/formatted/Results_Apr10/Orthogroups.txt)
  GoodProts=$WorkDir/goodProteins/goodProteins.fasta
  OutDir=$WorkDir/orthogroups_fasta
  mkdir -p $OutDir
  $ProgDir/orthoMCLgroups2fasta.py --orthogroups $OrthogroupsTxt --fasta $GoodProts --out_dir $OutDir
```

# 6) Downstream analysis

Particular orthogroups were analysed for expansion in isolates.

This section details the commands used and the results observed.


### 6.1 ) Clade unique gene families

The genes unique to A. tenuissima apple pathotypes were identified within the orthology analysis.

First variables were set:
```bash
  WorkDir=analysis/orthology/orthomcl/At_Aa_Ag_all_isolates
  Orthogroups=$(ls $WorkDir/formatted/Results_May04/Orthogroups.txt)
  GoodProts=$WorkDir/goodProteins/goodProteins.fasta
  # Braker_genes=
```

#### 6.1.a ) Orthologroups only containing A. tenuissima genes were extracted:

```bash
for num in 1; do
  # AtenUniq
AtenUniqDir=$WorkDir/A.tenuissima_unique
Uniq_Aten_groups=$AtenUniqDir/A.tenuissima_unique.txt
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
51
The following number genes from isolate 648 are contained in these orthogorups:
51
The following number genes from isolate 1082 are contained in these orthogorups:
51
The following number genes from isolate 1164 are contained in these orthogorups:
51
The following number genes from isolate 24350 are contained in these orthogorups:
52
The following number genes from isolate 648 are contained in these orthogorups:
51
The following number genes from isolate 743 are contained in these orthogorups:
51
The following number genes from isolate 1166 are contained in these orthogorups:
51
The following number genes from isolate 1177 are contained in these orthogorups:
51
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
150
The following number genes from isolate 675 are contained in these orthogorups:
155
The following number genes from isolate 97.0013 are contained in these orthogorups:
154
The following number genes from isolate 97.0016 are contained in these orthogorups:
156
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
327
The following number genes from isolate 650 are contained in these orthogorups:
329
```

#### 6.1.d ) Orthologroups only containing A. tenuissima non pathotype genes were extracted:

```bash
for num in 1; do
  # AtenNonPathUniq
  AtenNonPathUniqDir=$WorkDir/A.tenuissima_non_pathotype_unique
  Uniq_AtenNonPath_groups=$AtenNonPathUniqDir/A.tenuissima_non_pathotype_unique.txt
  mkdir -p $AtenNonPathUniqDir
  cat $Orthogroups | grep -e 'At_1' | grep -e 'At_2' | grep -e 'At_3' | grep -e 'At_4' | grep -v -e 'At_5' -e 'At_6' -e 'At_7' -e 'At_8' -e 'Aa' -e 'Ag' > $Uniq_AtenNonPath_groups
  echo "The number of orthogroups unique to A.tenuissima non-apple pathotype are:"
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
The number of orthogroups unique to A.tenuissima non-apple pathotype are:
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
49
The following number genes from isolate 635 are contained in these orthogorups:
56
The following number genes from isolate 743 are contained in these orthogorups:
57
The following number genes from isolate 1166 are contained in these orthogorups:
72
The following number genes from isolate 1177 are contained in these orthogorups:
55
```

```bash
cat $Uniq_AtenPath_groups | grep -o -P "At_7\|\S+" | cut -f2 -d '|' > tmp.txt
Isolate='1166'
AnnotTab=$(ls gene_pred/annotation/A.*/$Isolate/${Isolate}_annotation_ncbi.tsv)
OutFile=$(echo $AnnotTab | sed 's/_annotation_ncbi.tsv/_apple_pathotype_unique_orthogroups.tsv/g')
for Gene in $(cat tmp.txt); do
  cat $AnnotTab | grep -w "^$Gene"
done > $OutFile
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
48
The following number genes from isolate 650 are contained in these orthogorups:
62
The following number genes from isolate 648 are contained in these orthogorups:
58
The following number genes from isolate 743 are contained in these orthogorups:
58
The following number genes from isolate 1166 are contained in these orthogorups:
54
The following number genes from isolate 1177 are contained in these orthogorups:
60
```
