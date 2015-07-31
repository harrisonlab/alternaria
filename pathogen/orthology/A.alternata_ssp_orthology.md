# For a comparison between A. alternata, ssp. gaisen, spp. teniuissima, spp. tenuissima apple pathotype, ssp. arborescens


```bash
  ProjDir=/home/groups/harrisonlab/project_files/alternaria
  cd $ProjDir
  IsolateAbrv=Aalt_ssp_gai_ten_path_arb
  WorkDir=analysis/orthology/orthomcl/$IsolateAbrv
  mkdir -p $WorkDir
  mkdir -p $WorkDir/formatted
  mkdir -p $WorkDir/goodProteins
  mkdir -p $WorkDir/badProteins  
```

## Format fasta files

### for A.alt_ssp.gai 650
```bash
  Taxon_code=Agai
  Isolate1=gene_pred/augustus/A.alternata_ssp._gaisen/650/650_augustus_preds.aa
  Fasta_file=$WorkDir/"$Taxon_code"_concat.aa
  printf "" > $Fasta_file
  cat $Isolate1 | sed 's/>/>650_/g' >> $Fasta_file
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```

### for A.alt_ssp.ten 648, 1082, 1164, 24350
```bash
  Taxon_code=Aten
  Isolate1=gene_pred/augustus/A.alternata_ssp._tenuissima/648/648_augustus_preds.aa
  Isolate2=gene_pred/augustus/A.alternata_ssp._tenuissima/1082/1082_augustus_preds.aa
  Isolate3=gene_pred/augustus/A.alternata_ssp._tenuissima/1164/1164_augustus_preds.aa
  Isolate4=gene_pred/augustus/A.alternata_ssp._tenuissima/24350/24350_augustus_preds.aa
  Fasta_file=$WorkDir/"$Taxon_code"_concat.aa
  printf "" > $Fasta_file
  cat $Isolate1 | sed 's/>/>648_/g' >> $Fasta_file
  cat $Isolate2 | sed 's/>/>1082_/g' >> $Fasta_file
  cat $Isolate3 | sed 's/>/>1164_/g' >> $Fasta_file
  cat $Isolate4 | sed 's/>/>24350_/g' >> $Fasta_file
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```

### for A.alt_ssp.ten_apple_path 635, 743, 1166, 1177
```bash
  Taxon_code=AtAP
  Isolate1=gene_pred/augustus/A.alternata_ssp._tenuissima/635/635_augustus_preds.aa
  Isolate2=gene_pred/augustus/A.alternata_ssp._tenuissima/743/743_augustus_preds.aa
  Isolate3=gene_pred/augustus/A.alternata_ssp._tenuissima/1166/1166_augustus_preds.aa
  Isolate4=gene_pred/augustus/A.alternata_ssp._tenuissima/1177/1177_augustus_preds.aa
  Fasta_file=$WorkDir/"$Taxon_code"_concat.aa
  printf "" > $Fasta_file
  cat $Isolate1 | sed 's/>/>635_/g' >> $Fasta_file
  cat $Isolate2 | sed 's/>/>743_/g' >> $Fasta_file
  cat $Isolate3 | sed 's/>/>1166_/g' >> $Fasta_file
  cat $Isolate4 | sed 's/>/>1177_/g' >> $Fasta_file
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```

### for A.alt_ssp.arb 675, 97.0013, 97.0016
```bash
  Taxon_code=Aarb
  Isolate1=gene_pred/augustus/A.alternata_ssp._arborescens/675/675_augustus_preds.aa
  Isolate2=gene_pred/augustus/A.alternata_ssp._arborescens/97.0013/97.0013_augustus_preds.aa
  Isolate3=gene_pred/augustus/A.alternata_ssp._arborescens/97.0016/97.0016_augustus_preds.aa
  Fasta_file=$WorkDir/"$Taxon_code"_concat.aa
  printf "" > $Fasta_file
  cat $Isolate1 | sed 's/>/>675_/g' >> $Fasta_file
  cat $Isolate2 | sed 's/>/>97.0013_/g' >> $Fasta_file
  cat $Isolate3 | sed 's/>/>97.0016_/g' >> $Fasta_file
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```

## Filter proteins into good and poor sets.

```bash
  Input_dir=$WorkDir/formatted
  Min_length=10
  Max_percent_stops=20
  Good_proteins_file=$WorkDir/goodProteins/goodProteins.fasta
  Poor_proteins_file=$WorkDir/badProteins/poorProteins.fasta
  orthomclFilterFasta $Input_dir $Min_length $Max_percent_stops $Good_proteins_file $Poor_proteins_file
```

## Perform an all-vs-all blast of the proteins

```bash
  BlastDB=$WorkDir/blastall/$IsolateAbrv.db

  makeblastdb -in $Good_proteins_file -dbtype prot -out $BlastDB
  BlastOut=$WorkDir/all-vs-all_results.tsv
  mkdir -p $WorkDir/splitfiles

  SplitDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/signal_peptides
  $SplitDir/splitfile_500.py --inp_fasta $Good_proteins_file --out_dir $WorkDir/splitfiles --out_base goodProteins

  ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/orthology  
  for File in $(find $WorkDir/splitfiles); do
    Jobs=$(qstat | grep 'blast_500' | wc -l)
    while [ $Jobs -gt 32 ]; do
      sleep 10
      printf "."
      Jobs=$(qstat | grep 'blast_500' | wc -l)
    done
    printf "\n"
    echo $File
    BlastOut=$(echo $File | sed 's/.fa/.tab/g')
    qsub $ProgDir/blast_500.sh $BlastDB $File $BlastOut
  done
```

## Merge the all-vs-all blast results  
```bash  
  MergeHits="$IsolateAbrv"_blast.tab
  printf "" > $MergeHits
  for Num in $(ls $WorkDir/splitfiles/*.tab | rev | cut -f1 -d '_' | rev | sort -n); do
    File=$(ls $WorkDir/splitfiles/*_$Num)
    cat $File
  done > $MergeHits
```

## Perform ortholog identification

```bash
  ProgDir=~/git_repos/emr_repos/tools/pathogen/orthology/orthoMCL
  MergeHits="$IsolateAbrv"_blast.tab
  GoodProts=$WorkDir/goodProteins/goodProteins.fasta
  qsub $ProgDir/qsub_orthomcl.sh $MergeHits $GoodProts
```

## Plot venn diagrams:

```bash
  ProgDir=~/git_repos/emr_repos/tools/pathogen/orthology/venn_diagrams
  $ProgDir/venn_diag_4_way.r --inp $WorkDir/"$IsolateAbrv"_orthogroups.tab --out $WorkDir/"$IsolateAbrv"_orthogroups.pdf
```

Output was a pdf file of the venn diagram.

The following additional information was also provided. The format of the
following lines is as follows:

Isolate name (total number of orthogroups)
number of unique singleton genes
number of unique groups of inparalogs
