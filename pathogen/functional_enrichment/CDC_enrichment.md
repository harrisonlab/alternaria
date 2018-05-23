# Enrichment of annotation terms of Fus2 genes

## GO enrichment using topGO

### Chromosomal comparisons

GO enrichment of terms in PS contigs:

```bash
for Assembly in $(ls repeat_masked/*/*/filtered_contigs/*_contigs_softmasked_repeatmasker_TPSI_appended.fa); do
  Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
  Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
  OutDir=analysis/enrichment/$Organism/$Strain/CDC_vs_core
  mkdir -p $OutDir
  # Extract GO terms
  InterProTSV=gene_pred/interproscan/$Organism/$Strain/${Strain}_interproscan.tsv
  ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium/analysis/gene_enrichment
  $ProgDir/GO_prep_table.py --interpro $InterProTSV > $OutDir/${Strain}_gene_GO_annots.tsv

AnnotTable=$(ls gene_pred/annotation/$Organism/$Strain/${Strain}_annotation_ncbi.tsv)
AllGenes=$OutDir/${Strain}_all_genes.txt
cat $AnnotTable | tail -n+2 | cut -f1 > $AllGenes
Set1Genes=$OutDir/${Strain}_Set1_genes.txt
Set2Genes=$OutDir/${Strain}_Set2_genes.txt
AllGenes=$OutDir/${Strain}_all_genes.txt
if [ $Strain == 1166 ]; then
  cat $AnnotTable | tail -n+2 | grep -e 'contig_14' -e 'contig_15' -e 'contig_18' -e 'contig_19' -e 'contig_20' -e 'contig_21' | cut -f1 | sed -e 's/$/\t0.001/g'> $Set1Genes
  cat $AnnotTable | tail -n+2 | grep -v -e 'contig_14' -e 'contig_15' -e 'contig_18' -e 'contig_19' -e 'contig_20' -e 'contig_21' | cut -f1 | sed -e 's/$/\t1.00/g' > $Set2Genes
elif [ $Strain == 650 ]; then
  cat $AnnotTable | tail -n+2 | grep -e 'contig_12' -e 'contig_14' -e 'contig_19' -e 'contig_20' -e 'contig_21' -e 'contig_22' -e 'contig_23' | cut -f1 | sed -e 's/$/\t0.001/g'> $Set1Genes
  cat $AnnotTable | tail -n+2 | grep -v -e 'contig_12' -e 'contig_14' -e 'contig_19' -e 'contig_20' -e 'contig_21' -e 'contig_22' -e 'contig_23' | cut -f1 | sed -e 's/$/\t1.00/g' > $Set2Genes
fi
cat $Set1Genes $Set2Genes > $AllGenes

$ProgDir/GO_enrichment.r --all_genes $AllGenes --GO_annotations $OutDir/${Strain}_gene_GO_annots.tsv --out_dir $OutDir > $OutDir/output.txt
done
```

GO enrichment of terms in other LS contigs:

```bash
OutDir=analysis/enrichment/F.oxysporum_fsp_cepae/Fus2_canu_new/other_LS_vs_core
mkdir -p $OutDir
InterProTSV=gene_pred/interproscan/F.oxysporum_fsp_cepae/Fus2_canu_new/Fus2_canu_new_interproscan.tsv
ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium/analysis/gene_enrichment
$ProgDir/GO_prep_table.py --interpro $InterProTSV > $OutDir/Fus2_gene_GO_annots.tsv

AnnotTable=gene_pred/annotations/F.oxysporum_fsp_cepae/Fus2_canu_new/Fus2_canu_new_gene_annotations.tab
AllGenes=$OutDir/Fus2_all_genes.txt
cat $AnnotTable | tail -n+2  | cut -f1 > $AllGenes
Set1Genes=$OutDir/Fus2_Set1_genes.txt
Set2Genes=$OutDir/Fus2_Set2_genes.txt
AllGenes=$OutDir/Fus2_all_genes.txt
cat $AnnotTable | tail -n+2 | grep -e 'contig_14' -e 'contig_20' -e 'contig_22' | cut -f1 | sed -e 's/$/\t0.001/g'> $Set1Genes
cat $AnnotTable | tail -n+2 | grep -v -e 'contig_14' -e 'contig_20' -e 'contig_22' | cut -f1 | sed -e 's/$/\t1.00/g' > $Set2Genes
cat $Set1Genes $Set2Genes > $AllGenes

$ProgDir/GO_enrichment.r --all_genes $AllGenes --GO_annotations $OutDir/Fus2_gene_GO_annots.tsv --out_dir $OutDir > $OutDir/output.txt
```

## Interproscan term enrichment using Fishers exact test in R

```bash

for Assembly in $(ls repeat_masked/*/*/filtered_contigs/*_contigs_softmasked_repeatmasker_TPSI_appended.fa); do
  Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
  Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
  OutDir=analysis/enrichment/$Organism/$Strain/CDC_vs_core
  mkdir -p $OutDir
  AnnotTable=$(ls gene_pred/annotation/$Organism/$Strain/${Strain}_annotation_ncbi.tsv)

  CDC_set="contig_10_pilon contig_16_pilon contig_19_pilon contig_21_pilon"
  Core_set="contig_1_pilon contig_2_pilon contig_3_pilon contig_4_pilon contig_5_pilon contig_6_pilon contig_7_pilon"

  if [ $Strain == 1166 ]; then
    CDC_set="contig_14 contig_15 contig_18 contig_19 contig_20 contig_21"
    Core_set="contig_1 contig_2 contig_3 contig_4 contig_5 contig_6 contig_7 contig_8 contig_9 contig_10 contig_11 contig_12 contig_13 contig_16 contig_17 contig_22"
  elif [ $Strain == 650 ]; then
    CDC_set="contig_12 contig_14 contig_19 contig_20 contig_21 contig_22 contig_23"
    Core_set="contig_1 contig_2 contig_3 contig_4 contig_5 contig_6 contig_7 contig_8 contig_9 contig_10 contig_11 contig_13 contig_15 contig_16 contig_17 contig_18"
  fi

  mkdir -p $OutDir/tables
  ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium/analysis/gene_enrichment
  $ProgDir/IPR_prep_tables.py --interpro $AnnotTable --set1_name LS_contigs --set2_name core_contigs --contig_set1 $CDC_set --contig_set2 $Core_set --outdir $OutDir/tables
  $ProgDir/fisherstest.r --in_dir $OutDir/tables  --out_file $OutDir/results.tsv
  # Run fishertest.r commands
  cat $OutDir/results.tsv | cut -f1,3 | grep -e "e-" -e "0\.00" -e "0\.01" -e "0\.02" -e "0\.03" -e "0\.04" > $OutDir/significant_terms.txt
  SigIPR=$(cat $OutDir/significant_terms.txt | cut -f1)
  for IPR in $SigIPR; do
    cat $OutDir/tables/"$IPR"_fischertable.txt | grep "$IPR"
  done > $OutDir/${Strain}_IPR_sigterms.txt
done
```

Taken from: http://www.cs.tau.ac.il/~rshamir/ge/09/scribe/lec14a.pdf
'Since the signi􏰃cance test is performed for many groups, a multiple testing correction must be carried out in order to limit false positives. Both the Bonferroni and FDR methods are too stringent since there exist strong dependencies between groups (since they are often members of the same hierarchy). To get around these limitations, TANGO instead calculates the empirical p value distribution. For a given cluster Tj, TANGO samples many random sets of the same size & computes their p-values vs. each of the annotation sets Ai. Next, it also permutes gene IDs to eliminate dependency between annotation sets and target sets. This correction also applies for testing multiple clusters.'


```bash
cp gene_pred/interproscan/F.oxysporum_fsp_cepae/Fus2_canu_new/Fus2_canu_new_interproscan.tsv $OutDir/.
cp gene_pred/annotations/F.oxysporum_fsp_cepae/Fus2_canu_new/Fus2_canu_new_gene_annotations.tab $OutDir/.
```
