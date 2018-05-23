# Enrichment of annotation terms of Fp genes

## GO enrichment using topGO


GO enrichment of terms in Non-syn vs background:

```bash
  Prefix=non-syn_vs_background
  OutDir=analysis/enrichment/A.alternata_ssp_tenuissima/1166/$Prefix
  mkdir -p $OutDir
  # Append_InterProTSV=$OutDir/${Prefix}_interproscan.tsv
  # Fp_InterProTSV=$(ls gene_pred/interproscan/F.proliferatum/A8_ncbi/A8_ncbi_interproscan.tsv)
  # cat $Fp_InterProTSV | sed -e 's/^/Fp_/g' > $Append_InterProTSV
  # FoC_InterProTSV=$(ls gene_pred/interproscan/F.proliferatum/A8_ncbi/A8_ncbi_interproscan.tsv)
  # cat $FoC_InterProTSV | sed -e 's/^/FoC_/g' >> $Append_InterProTSV

  Append_InterProTSV=$(ls gene_pred/interproscan/A.alternata_ssp_tenuissima/1166/1166_interproscan.tsv)

  ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium/analysis/gene_enrichment
  $ProgDir/GO_prep_table.py --interpro $Append_InterProTSV > $OutDir/${Prefix}_GO_annots.tsv

  cat analysis/popgen/SNP_calling/arborescens_vs_1166/arborescens_vs_1166_filtered_no_indels.recode_nonsyn.vcf | grep -v '#' | cut -f8 | cut -f7 -d '|' | sort | uniq > analysis/popgen/SNP_calling/arborescens_vs_1166/arborescens_vs_1166_filtered_no_indels.recode_nonsyn_names.txt
  NonSynGenes=$(ls analysis/popgen/SNP_calling/arborescens_vs_1166/arborescens_vs_1166_filtered_no_indels.recode_nonsyn_names.txt)

  Set1Genes=$OutDir/non-syn_genes.txt
  Set2Genes=$OutDir/background_genes.txt
  cat $OutDir/${Prefix}_GO_annots.tsv | cut -f1 | grep -f $NonSynGenes | sed -e 's/$/\t0.001/g' > $Set1Genes
  cat $OutDir/${Prefix}_GO_annots.tsv | cut -f1 | grep -v -f $NonSynGenes | sed -e 's/$/\t1.00/g' > $Set2Genes
  AllGenes=$OutDir/${Prefix}_all_genes.txt
  cat $Set1Genes $Set2Genes > $AllGenes

  $ProgDir/GO_enrichment.r --all_genes $AllGenes --GO_annotations $OutDir/${Prefix}_GO_annots.tsv --out_dir $OutDir > $OutDir/output.txt
```

## Interproscan term enrichment using Fishers exact test in R

IPR enrichment of terms in Fp vs FoC:

```bash
  Prefix=non-syn_vs_background
  OutDir=analysis/enrichment/A.alternata_ssp_tenuissima/1166/$Prefix
  mkdir -p $OutDir/tables

  AnnotTab=$(ls gene_pred/annotation/A.alternata_ssp_tenuissima/1166/1166_annotation_ncbi.tsv)
  cat $Fp_AnnotTable | head -n1 > $Appended_AnnotTable
  cat $Fp_AnnotTable | tail -n+2 | sed -e "s/^/Fp_/g" | sed 's/contig/Fp_contig/g' >> $Appended_AnnotTable
  cat $FoC_AnnotTable | tail -n+2 | sed -e "s/^/FoC_/g" | sed 's/contig/FoC_contig/g' >> $Appended_AnnotTable

  Set1Contigs=$(cat $Appended_AnnotTable | tail -n+2 | cut -f2 | grep 'Fp_' | sort | uniq | tr -d '\n' | sed 's/Fp/ Fp/g')
  Set2Contigs=$(cat $Appended_AnnotTable | tail -n+2 | cut -f2 | grep 'FoC_' | sort | uniq | tr -d '\n' | sed 's/FoC/ FoC/g')

  # Fp vs FoC

  ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium/analysis/gene_enrichment
  $ProgDir/IPR_prep_tables.py --interpro $Appended_AnnotTable --set1_name Fp_genes --set2_name FoC_genes --contig_set1 $Set1Contigs  --contig_set2 $Set2Contigs --outdir $OutDir/tables
  $ProgDir/fisherstest.r --in_dir $OutDir/tables  --out_file $OutDir/results.tsv
  # Run fishertest.r commands
  cat $OutDir/results.tsv | cut -f1,2 | grep -e "e-" -e "0\.00" -e "0\.01" -e "0\.02" -e "0\.03" -e "0\.04"  > $OutDir/significant_terms.txt
  SigIPR=$(cat $OutDir/significant_terms.txt | cut -f1)
  for IPR in $SigIPR; do
    # echo $IPR
    cat $OutDir/tables/"$IPR"_fischertable.txt | grep "$IPR"
  done
```

Taken from: http://www.cs.tau.ac.il/~rshamir/ge/09/scribe/lec14a.pdf
'Since the signi􏰃cance test is performed for many groups, a multiple testing correction must be carried out in order to limit false positives. Both the Bonferroni and FDR methods are too stringent since there exist strong dependencies between groups (since they are often members of the same hierarchy). To get around these limitations, TANGO instead calculates the empirical p value distribution. For a given cluster Tj, TANGO samples many random sets of the same size & computes their p-values vs. each of the annotation sets Ai. Next, it also permutes gene IDs to eliminate dependency between annotation sets and target sets. This correction also applies for testing multiple clusters.'
