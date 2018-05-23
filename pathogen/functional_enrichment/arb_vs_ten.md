
Approach
For each IPR annotation in all provided A.ten and A.arb genome count the number




```bash
  OutDir=analysis/enrichment/arb_vs_ten
  mkdir -p $OutDir
  ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium/analysis/gene_enrichment
  $ProgDir/IPR_prep_tables.py --interpro $AnnotTable --set1_name LS_contigs --set2_name core_contigs --contig_set1 $CDC_set --contig_set2 $Core_set --outdir $OutDir/tables
  $ProgDir/fisherstest.r --in_dir $OutDir/tables  --out_file $OutDir/results.tsv
  # Run fishertest.r commands
```

```bash

AtenIsolates='675 1082 1164 24350'
AtenIPR=$(ls gene_pred/interproscan/A.*/*/*_interproscan.tsv | grep -e '675' -e '1082' -e '1164' -e '24350')
AarbIsolates='648 97.0013 97.0016'
AarbIPR=$(ls gene_pred/interproscan/A.*/*/*_interproscan.tsv | grep -e '648' -e '97.0013' -e '97.0016')
OutDir=analysis/enrichment/arb_vs_ten
mkdir -p $OutDir
ProgDir=/home/armita/git_repos/emr_repos/scripts/alternaria/pathogen/functional_enrichment
$ProgDir/Alt_species_enrichment.py --isolate_list_A $AtenIsolates --isolate_list_B $AarbIsolates --interpro_list_A $AtenIPR --interpro_list_B $AarbIPR > $OutDir/IPR_copies_by_isolate.txt
```

Output was downloaded to my local machine. It was then analysed using Rstudio. The commands are listed below:

```R
library(readr)
IPR_copies_by_isolate <- read_delim("~/Downloads/Aalt/enrichment/IPR_copies_by_isolate.txt",
     "\t", escape_double = FALSE, trim_ws = TRUE)
samp.with.rownames <- data.frame(samp[,-1], row.names=samp[,1])
df1 <- data.frame(IPR_copies_by_isolate)
rownames(df1) <- df1[,1]
df1$IPR_ID <- NULL
$IPR_description <- NULL
source("https://bioconductor.org/biocLite.R")
biocLite("HybridMTest")
library(HybridMTest)
grplbl <- c('Aten', 'Aten', 'Aten', 'Aten', 'Aarb', 'Aarb', 'Aarb')
df_test_results <- data.frame(row.oneway.anova(df1, grplbl))
df2 <- data.frame(IPR_copies_by_isolate)
rownames(df2) <- df2[,1]
df2$IPR_ID <- NULL
df_sigIPR <- subset(df_test_results, pval < 0.05)
df_merged_sigIPR <- merge(df_sigIPR, df2, by.x = 0, by.y = 0, by.all = x)
```


21 domains were identified as enriched domains.

These were:

```bash
OutDir=analysis/enrichment/arb_vs_ten
mkdir -p $OutDir
printf "IPR000120
IPR001663
IPR002067
IPR002401
IPR002403
IPR002523
IPR002925
IPR003034
IPR004803
IPR006598
IPR011051
IPR013078
IPR014710
IPR014851
IPR015881
IPR019786
IPR021858
IPR023631
IPR024079
IPR029033
IPR031348" > $OutDir/enriched_domains.txt

AnnotTab=$(ls gene_pred/annotation/A.*/*/*_annotation_ncbi.tsv | grep '1166')
cat $AnnotTab | grep -f $OutDir/enriched_domains.txt > $OutDir/enriched_domains_annot.tsv
cat $OutDir/enriched_domains_annot.tsv | grep ''
ls $PWD/$OutDir/enriched_domains_annot.tsv
cat $OutDir/enriched_domains_annot.tsv | cut -f20,21 | sort | uniq > $OutDir/enriched_domains_orthogroups.tsv
ls $PWD/$OutDir/enriched_domains_orthogroups.tsv
```
