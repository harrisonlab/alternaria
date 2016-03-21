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
qlogin
WorkDir=/tmp/genome_alignment
mkdir -p $WorkDir
cd $WorkDir
ProjDir=/home/groups/harrisonlab/project_files/alternaria
Assembly_648=$(ls $ProjDir/repeat_masked/*/648/filtered_contigs_repmask/*_contigs_softmasked.fa)
Assembly_1082=$(ls $ProjDir/repeat_masked/*/1082/filtered_contigs_repmask/*_contigs_softmasked.fa)
Assembly_1164=$(ls $ProjDir/repeat_masked/*/1164/filtered_contigs_repmask/*_contigs_softmasked.fa)
Assembly_24350=$(ls $ProjDir/repeat_masked/*/24350/filtered_contigs_repmask/*_contigs_softmasked.fa)
Assembly_635=$(ls $ProjDir/repeat_masked/*/635/filtered_contigs_repmask/*_contigs_softmasked.fa)
Assembly_743=$(ls $ProjDir/repeat_masked/*/743/filtered_contigs_repmask/*_contigs_softmasked.fa)
Assembly_1166=$(ls $ProjDir/repeat_masked/*/1166/filtered_contigs_repmask/*_contigs_softmasked.fa)
Assembly_1177=$(ls $ProjDir/repeat_masked/*/1177/filtered_contigs_repmask/*_contigs_softmasked.fa)
Assembly_675=$(ls $ProjDir/repeat_masked/*/675/filtered_contigs_repmask/*_contigs_softmasked.fa)
Assembly_970013=$(ls $ProjDir/repeat_masked/*/97.0013/filtered_contigs_repmask/*_contigs_softmasked.fa)
Assembly_970016=$(ls $ProjDir/repeat_masked/*/97.0016/filtered_contigs_repmask/*_contigs_softmasked.fa)
Assembly_650=$(ls $ProjDir/repeat_masked/*/650/filtered_contigs_repmask/*_contigs_softmasked.fa)
cat $Assembly_648 | sed 's/>/>648:/g' | sed '/>/s/$/:1:+:1/' > 648
cat $Assembly_1082 | sed 's/>/>1082:/g' | sed '/>/s/$/:1:+:1/' > 1082
cat $Assembly_1164 | sed 's/>/>1164:/g' | sed '/>/s/$/:1:+:1/' > 1164
cat $Assembly_24350 | sed 's/>/>24350:/g' | sed '/>/s/$/:1:+:1/' > 24350
cat $Assembly_635 | sed 's/>/>635:/g' | sed '/>/s/$/:1:+:1/' > 635
cat $Assembly_743 | sed 's/>/>743:/g' | sed '/>/s/$/:1:+:1/' > 743
cat $Assembly_1166 | sed 's/>/>1166:/g' | sed '/>/s/$/:1:+:1/' > 1166
cat $Assembly_1177 | sed 's/>/>1177:/g' | sed '/>/s/$/:1:+:1/' > 1177
cat $Assembly_675 | sed 's/>/>675:/g' | sed '/>/s/$/:1:+:1/' > 675
cat $Assembly_970013 | sed 's/>/>97.0013:/g' | sed '/>/s/$/:1:+:1/' > 97.0013
cat $Assembly_970016 | sed 's/>/>97.0016:/g' | sed '/>/s/$/:1:+:1/' > 97.0016
all_bz - \
  "((648 1082 1164 24350 (635 743 1166 1177)) (675 97.0013 97.0016) (650))" >& all_bz.log

  #  all_bz - \
  #    "((648 1082 1164 24350 (635 743 1166 1177)) (675 97.0013 97.0016) (650))" \
  #     | sed 's/maf_sort/maf_sort.py/g' >& all_bz.log
# +(off) verbose
# -(off) output command only.
# b(2) 0: run post-process only 1: run blastzWrapper only, transform to maf 2: run both
# A(1) 0: toast 1: single_cov2 2: toast, following by chain and single cov on reference
# F(null) null: single coverage is done for both species; reference: single coverage is done for reference only, effective in single_cov2
# T(null): annotation file path and name, used for running toast and chaining procedure
# h(300) minimum chaining size, effective in toast
# q(600) minimum cluster size, effective in toast
# D(1) 0: run all_bz for roast 1: run all_bz for TBA.
# c(500): parameter transfered to blastz_clean, alignments closer than c are subjected to be cleaned.
# f(2) x% is used for determine in-paralogs, effective in toast.

```

### 2.2) generating the multiple alignment

### 2.3) “projecting” the alignment onto a reference sequence

## 3) Process alignment file

The Alignment file was processed using Maffilter.
http://biopp.univ-montp2.fr/forge/maffilter
This also allowed calling of SNPs, extraction of CDS and preparation for
phylogenetic analysis.
Importantly, it allows identification of contigs that do not align to other
strains, facilitating identification of CDCs.

```bash

```
