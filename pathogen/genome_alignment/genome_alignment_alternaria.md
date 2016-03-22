# Genome alignment

Genome alignment was performed using the pipeline recomended by Eva Stukenbrock
at the Max Planck Institute for Evolutionary Biology, Kiel.

## 1) Generation of a phylogenetic tree for isolates

A phylogenetic tree was provided to TBA

```bash
TreeString='"((648 1082 1164 24350 (635 743 1166 1177)) (675 97.0013 97.0016) (650))"'
```

## 2) Generate multiple alignment between genomes

The TBA package http://www.bx.psu.edu/miller_lab/ was used to generate multiple
alignments.

### 2.1) generating a series of pair-wise alignments

generating a series of pair-wise alignments to “seed” the multiple alignment process

```bash
  ProjDir=/home/groups/harrisonlab/project_files/alternaria
  WorkDir=analysis/genome_alignment/tmp
  mkdir -p $WorkDir
  cd $WorkDir

  for Assembly in $(ls $ProjDir/repeat_masked/*/*/filtered_contigs_repmask/*_contigs_softmasked.fa); do
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
    Strain=$(echo $Assembly| rev | cut -f3 -d '/' | rev)
    cat $Assembly | sed "s/>/>$Strain:/g" > $Strain.fa
    get_standard_headers $Strain.fa > "$Strain"_headers.txt
    ProgDir=~/git_repos/emr_repos/tools/pathogen/lineage_specific_regions
    $ProgDir/parse_tba_headers.py --inp_fasta $Strain.fa --new_headers "$Strain"_headers.txt > $Strain
  done

  all_bz - \
    "((648 1082 1164 24350 (635 743 1166 1177)) (675 97.0013 97.0016) (650))" >& all_bz.log

  cat all_bz.log | grep 'blastzWrapper' > commands_part1.log
  while read Commands; do
    qsub -S /bin/bash -b y -pe smp 1 -l virtual_free=1G -N blastzWrapper -cwd "$Commands"
  done < commands_part1.log

  cat all_bz.log | grep 'single_cov2' > commands_part2.log
  while read Commands; do
    qsub -S /bin/bash -b y -pe smp 1 -l virtual_free=1G -N blastzWrapper -cwd "$Commands"
  done < commands_part2.log

```

### 2.2) generating the multiple alignment

```bash
tba \
  "((648 1082 1164 24350 (635 743 1166 1177)) (675 97.0013 97.0016) (650))" \
  *.*.maf tba.maf >& tba.log
```

### 2.3) “projecting” the alignment onto a reference sequence

## 3) Process alignment file

The Alignment file was processed using Maffilter.
http://biopp.univ-montp2.fr/forge/maffilter
This also allowed calling of SNPs, extraction of CDS and preparation for
phylogenetic analysis.
Importantly, it allows identification of contigs that do not align to other
strains, facilitating identification of CDCs.

```bash
  # maffilter input.file=mydata.maf.gz input.file.compression=gzip input.format=Maf output.log=mydata.maffilter.log
  maffilter input.file=tba.maf input.file.compression=none input.format=Maf output.log=mydata.maffilter.log
```
