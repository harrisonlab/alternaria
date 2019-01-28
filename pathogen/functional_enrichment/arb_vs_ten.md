
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

AtenIsolates='648 1082 1164 24350'
AtenIPR=$(ls gene_pred/interproscan/A.*/*/*_interproscan.tsv | grep -e '648' -e '1082' -e '1164' -e '24350')
AarbIsolates='675 97.0013 97.0016'
AarbIPR=$(ls gene_pred/interproscan/A.*/*/*_interproscan.tsv | grep -e '675' -e '97.0013' -e '97.0016')
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
# samp.with.rownames <- data.frame(samp[,-1], row.names=samp[,1])
df1 <- data.frame(IPR_copies_by_isolate)
rownames(df1) <- df1[,1]
df1$IPR_ID <- NULL
df1$IPR_description <- NULL
# source("https://bioconductor.org/biocLite.R")
# biocLite("HybridMTest")
library(HybridMTest)
library(genefilter)
grplbl <- c('Aten', 'Aten', 'Aten', 'Aten', 'Aarb', 'Aarb', 'Aarb')

# t.result <- apply(df1[,2:7], 1, function (x) t.test(x[1:4],x[4:7],paired=FALSE))
# counts$p_value <- unlist(lapply(t.result, function(x) x$p.value))

# df_test_results <- data.frame(row.oneway.anova(df1, grplbl))
df_test_results <- data.frame(rowttests(as.matrix(df1), factor(grplbl)))
df2 <- data.frame(IPR_copies_by_isolate)
rownames(df2) <- df2[,1]
df2$IPR_ID <- NULL
# df_sigIPR <- subset(df_test_results, pval < 0.05)
df_sigIPR <- subset(df_test_results, p.value < 0.05)
df_merged_sigIPR <- merge(df_sigIPR, df2, by.x = 0, by.y = 0, by.all = x)
df_merged_sigIPR$padj <- round(p.adjust(df_merged_sigIPR$p.value, method = "hochberg"), 3)
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


# Enrichment of transposon families

```r
TPSI_summary <- read.delim("~/Downloads/Aalt/enrichment/transposons/TPSI_summary.txt")
df1 <- subset(TPSI_summary, Organism == 'A.alternata_ssp._arborescens' | Organism == 'A.alternata_ssp._tenuissima')
df1$Strain <- NULL
df1$ISC1316 <- NULL
df1$Crypton <- NULL
library("ggpubr")
ggdensity(df1$DDE_1,
          main = "",
          xlab = "") + facet_wrap(df1$Organism)
ggqqplot(df1$DDE_1) #+ facet_wrap(df1$Organism)
shapiro.test(df1$DDE_1)
t <- data.frame(t(sapply(df1[-1], function(x)
      unlist(t.test(x~df1$Organism)[c("estimate","p.value","statistic","conf.int")]))))
round(p.adjust(t$p.value, method = "hochberg"), 3)
round(t$statistic, 2)

library(reshape2)
library(ggplot2)
library(scales)
library("ggpubr")
df2 <- melt(df1, id.vars=1)

# df2$combo <- paste(df2$variable, df2$Organism, sep="_")
# meds <- c(by(df2$value, df2$combo, median))

compare_means(value ~ Organism, data = df2,
              group.by = "variable", method="t.test")


# p <- ggplot(df2, aes(x=Organism, y=value)) +
#   scale_y_continuous(name = "", breaks= pretty_breaks()) +
#   geom_boxplot()
# p <- p + theme(axis.text.x=element_text(angle = 0, hjust = 0))
# p <- p + scale_x_discrete(name ="", labels=c("arb", "ten"))
# p <- p + facet_wrap(~df2$variable, nrow = 4, ncol = 3, strip.position = "top", scales="free_y")
# # p + geom_text(data=data.frame(), aes(x=names(meds), y=meds, label=1:20), col='red', size=10)
# # p + geom_text(size = 3,aes(x=Organism, y=1,
# #                  label=c("a", "b", "a", "b", "", "", "",
# #                    "", "", "",
# #                    "", "", "", "", "", "", "",
# #                      "", "", "")),
# #                                    data=df3)
# p <- p + stat_compare_means(method="t.test")
# outfile='transposon_enrichment.pdf'
# ggsave(outfile , plot = p, device = 'pdf', width =10, height = 10, units = "cm", limitsize = FALSE)
#

# --- Summary SE function ---
## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
    library(plyr)

    # New version of length which can handle NA's: if na.rm==T, don't count them
    length2 <- function (x, na.rm=FALSE) {
        if (na.rm) sum(!is.na(x))
        else       length(x)
    }

    # This does the summary. For each group's data frame, return a vector with
    # N, mean, and sd
    datac <- ddply(data, groupvars, .drop=.drop,
      .fun = function(xx, col) {
        c(N    = length2(xx[[col]], na.rm=na.rm),
          mean = mean   (xx[[col]], na.rm=na.rm),
          sd   = sd     (xx[[col]], na.rm=na.rm)
        )
      },
      measurevar
    )

    # Rename the "mean" column
    datac <- rename(datac, c("mean" = measurevar))

    datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean

    # Confidence interval multiplier for standard error
    # Calculate t-statistic for confidence interval:
    # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
    ciMult <- qt(conf.interval/2 + .5, datac$N-1)
    datac$ci <- datac$se * ciMult

    return(datac)
}
# --- /Summary SE function ends ---


df3 <- summarySE(df2, measurevar="value", groupvars=c("Organism", 'variable'))
df3$label_pos <- ((df3$value + df3$se) * 1.3)
# library(dplyr)
# df4 <- arrange(df3,desc(df3$'variable'),)
df3 <- df3[order(df3$variable),]

p<-ggplot(data=subset(df3), aes(x=Organism, y=value))
p <- p + geom_bar(stat="identity")
p <- p + theme(axis.text.x=element_text(angle = -45, hjust = 0))
p <- p + scale_y_continuous(name = "", breaks= pretty_breaks())
p <- p + scale_x_discrete(name ="", labels=c("arb", "ten"))
# p <- p + ylab('Number of lesions') + xlab('')
p <- p + geom_errorbar(aes(ymin=value-se, ymax=value+se),
                  width=.2,                    # Width of the error bars
                  position=position_dodge(.9))
# p <- p + geom_text(aes(x=name, y=1.2 * max(no..leisions), label=df4$sig),
#                                    data=)
 p <- p + geom_text(size = 3,aes(x=Organism, y=df3$label_pos,
                  label=c("a", "b", "a", "b", "", "", "",
                    "", "", "",
                    "", "", "", "", "", "", "",
                      "", "", "")),
                                    data=df3)
p <- p + facet_wrap("variable", nrow = 4, ncol = 3, strip.position = "top", scales="free_y")
outfile='transposon_enrichment.pdf'
ggsave(outfile , plot = p, device = 'pdf', width =15, height = 15, units = "cm", limitsize = FALSE)




df1 <- subset(TPSI_summary, Organism == 'A.alternata_ssp._arborescens' | Organism == 'A.alternata_ssp._tenuissima' | Organism == 'A.alternata_ssp._tenuissima apple pathotype')
df1$Strain <- NULL
df1$ISC1316 <- NULL
df1$Crypton <- NULL
ggdensity(df1$DDE_1,
          main = "",
          xlab = "") + facet_wrap(df1$Organism)
ggqqplot(df1$DDE_1) #+ facet_wrap(df1$Organism)
shapiro.test(df1$DDE_1)
# t <- data.frame(t(sapply(df1[-1], function(x)
#       unlist(t.test(x~df1$Organism)[c("estimate","p.value","statistic","conf.int")]))))
# round(p.adjust(t$p.value, method = "hochberg"), 3)
# round(t$statistic, 2)
df2 <- melt(df1, id.vars=1)
compare_means(value ~ Organism, data = df2,
              group.by = "variable", method="anova")

df3 <- summarySE(df2, measurevar="value", groupvars=c("Organism", 'variable'))
df3$label_pos <- ((df3$value + df3$se) * 1.3)
# library(dplyr)
# df4 <- arrange(df3,desc(df3$'variable'),)
df3 <- df3[order(df3$variable),]


p<-ggplot(data=subset(df3), aes(x=Organism, y=value))
p <- p + geom_bar(stat="identity")
# p <- p + geom_boxplot()
p <- p + theme(axis.text.x=element_text(angle = -45, hjust = 0))
p <- p + scale_y_continuous(name = "", breaks= pretty_breaks())
p <- p + scale_x_discrete(name ="", labels=c("arb", "ten"))
# p <- p + ylab('Number of lesions') + xlab('')
p <- p + geom_errorbar(aes(ymin=value-se, ymax=value+se),
                  width=.2,                    # Width of the error bars
                  position=position_dodge(.9))
p <- p + facet_wrap("variable", nrow = 4, ncol = 3, strip.position = "top", scales="free_y")
# p <- p + stat_compare_means(method="anova")
p <- p + stat_compare_means(method="anova", hide.ns='True')

```

This showed that gypsy elements and DDE elements were expanded in the arborescens
genomes.


Genes in 2kb of these elements were identified:

```bash
for TpsiGff in $(ls repeat_masked/A.alternata_ssp._arborescens/675/ncbi_edits_repmask/675_contigs_unmasked.fa.TPSI.allHits.chains.bestPerLocus.gff3); do
  Strain=$(echo $TpsiGff | rev | cut -f3 -d '/' | rev)
  Organism=$(echo $TpsiGff | rev | cut -f4 -d '/' | rev)
  echo "$Organism - $Strain"
OutDir=analysis/enrichment/arb_vs_ten/transposon/$Organism/$Strain
mkdir -p $OutDir
GeneGff=$(ls gene_pred/final/$Organism/$Strain/final/final_genes_appended_renamed.gff3)
AnnotTab=$(ls gene_pred/annotation/$Organism/$Strain/${Strain}_annotation_ncbi.tsv)
for Element in gypsy DDE_1; do
echo "$Element"
cat $TpsiGff | grep "${Element}" > $OutDir/${Element}_tpsi.gff
cat $OutDir/${Element}_tpsi.gff | wc -l
ProgDir="/home/armita/git_repos/emr_repos/tools/pathogen/mimp_finder"
$ProgDir/gffexpander.pl +- 2000 $OutDir/${Element}_tpsi.gff > $OutDir/"$Strain"_TPSI_${Element}_exp.gff

bedtools intersect -u -a $GeneGff -b $OutDir/"$Strain"_TPSI_${Element}_exp.gff > $OutDir/${Strain}_genes_in_2kb_${Element}_incl_transposons.gff

bedtools intersect -v -a $OutDir/${Strain}_genes_in_2kb_${Element}_incl_transposons.gff -b $OutDir/${Element}_tpsi.gff > $OutDir/${Strain}_genes_in_2kb_${Element}.gff

cat $OutDir/${Strain}_genes_in_2kb_${Element}.gff | grep -w 'gene' | cut -f9 | cut -f2 -d '=' | tr -d ';' > $OutDir/${Strain}_genes_in_2kb_${Element}_headers.txt
cat $OutDir/${Strain}_genes_in_2kb_${Element}_headers.txt | wc -l
# cat $AnnotTab | grep -w -f $OutDir/${Strain}_genes_in_2kb_${Element}_headers.txt > $OutDir/${Strain}_genes_in_2kb_${Element}_table.tsv
ProgDir=/home/armita/git_repos/emr_repos/scripts/alternaria/pathogen/annotation
$ProgDir/annot_tab_extract_genes.py --annot $AnnotTab --genes $OutDir/${Strain}_genes_in_2kb_${Element}_headers.txt > $OutDir/${Strain}_genes_in_2kb_${Element}_table.tsv
cat $OutDir/${Strain}_genes_in_2kb_${Element}_table.tsv | wc -l
done
done
```

# Testing for enrichment by orthogroups

```bash
OutDir=analysis/enrichment/arb_vs_ten/orthogroups
mkdir -p $OutDir
Orthogroups=$(ls analysis/orthology/orthomcl/At_Aa_Ag_all_isolates/formatted/Results_May31/Orthogroups.txt)
All="Aa_1 Aa_2 Aa_3 Ag_1 At_1 At_2 At_3 At_4 At_5 At_6 At_7 At_8"
Species="A.alternata_ssp._arborescens A.alternata_ssp._arborescens A.alternata_ssp._arborescens A.gaisen_pear_pathotype A.alternata_ssp._tenuissima A.alternata_ssp._tenuissima A.alternata_ssp._tenuissima A.alternata_ssp._tenuissima A.alternata_ssp._tenuissima_apple_pathotype A.alternata_ssp._tenuissima_apple_pathotype A.alternata_ssp._tenuissima_apple_pathotype A.alternata_ssp._tenuissima_apple_pathotype"
# Set1="Aa_1 Aa_2 Aa_3"
# Set2="At_1 At_2 At_3 At_4"
ProgDir=/home/armita/git_repos/emr_repos/scripts/alternaria/pathogen/functional_enrichment
$ProgDir/summarise_orthogroups_alt.py --orthogroups $Orthogroups --orthomclIDs $All --species $Species --prefix $OutDir/summarised_orthogroups
# ProgDir=/home/armita/git_repos/emr_repos/scripts/verticillium_pathogenomics
# $ProgDir/summarise_orthogroups.py --orthogroups $Orthogroups --set1 $Grp1 --set2 $Grp2 --all $All --effectors $Effectors
# > $OutDir/effector_orthogroups_summary.txt

ls $PWD/$OutDir/summarised_orthogroups.tsv

cat $OutDir/summarised_orthogroups_expansion_status.tsv | grep 'arborescens' | cut -f1 > $OutDir/arborescens_expanded_orthogroups.txt
cat $OutDir/arborescens_expanded_orthogroups.txt | wc -l
cat $OutDir/summarised_orthogroups_expansion_status.tsv | grep 'tenuissima'  | cut -f1 > $OutDir/tenuissima_expanded_orthogroups.txt
cat $OutDir/tenuissima_expanded_orthogroups.txt | wc -l
```
Investigate arborescens orthogroups
```bash
OutDir=analysis/enrichment/arb_vs_ten/orthogroups
AnnotTab=$(ls /home/groups/harrisonlab/project_files/alternaria/gene_pred/annotation/*/*/*_annotation_ncbi.tsv | grep '675')
cat $AnnotTab | grep -wf $OutDir/arborescens_expanded_orthogroups.txt > $OutDir/675_annotab_arborescens_expanded.tsv
cat $OutDir/675_annotab_arborescens_expanded.tsv | wc -l
cat $OutDir/675_annotab_arborescens_expanded.tsv | grep 'EffP' | wc -l
cat $OutDir/675_annotab_arborescens_expanded.tsv | grep 'At_1(0):At_2(0):At_3(0):At_4(0)' | wc -l

cat $AnnotTab | grep 'HET' | wc -l
cat $OutDir/675_annotab_arborescens_expanded.tsv | grep 'HET' | wc -l
#cat $AnnotTab | grep 'HET' | grep -v 'IPR010730' | less -S

cat $OutDir/675_annotab_arborescens_expanded.tsv  | grep 'CAZY' | wc -l

cat $OutDir/675_annotab_arborescens_expanded.tsv | grep -e 'Transcription factor' -e 'Zinc finger'  | grep 'At_1(0):At_2(0):At_3(0):At_4(0)' | wc -l

cat $AnnotTab | cut -f1 > $OutDir/arborescens_expanded_orthogroups_genes.txt

GeneGff=$(ls /home/groups/harrisonlab/project_files/alternaria/gene_pred/final/A.alternata_ssp._arborescens/675/final/final_genes_appended_renamed.gff3)
cat $GeneGff | grep -w 'mRNA' | grep -w -f $OutDir/arborescens_expanded_orthogroups_genes.txt > $OutDir/arborescens_expanded_orthogroup_genes.gff
#
# ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
# $ProgDir/extract_gff_for_sigP_hits.pl $OutFileHeaders $Gff effectorP ID > $EffectorP_Gff

```

Investigate tenuissima orthogroups
```bash
OutDir=analysis/enrichment/arb_vs_ten/orthogroups
AnnotTab=$(ls /home/groups/harrisonlab/project_files/alternaria/gene_pred/annotation/*/*/*_annotation_ncbi.tsv | grep '648')
cat $AnnotTab | grep -wf $OutDir/tenuissima_expanded_orthogroups.txt > $OutDir/648_annotab_tenuissima_expanded.tsv
cat $OutDir/648_annotab_tenuissima_expanded.tsv | wc -l
cat $OutDir/648_annotab_tenuissima_expanded.tsv | grep 'EffP' | wc -l
cat $OutDir/648_annotab_tenuissima_expanded.tsv | grep 'Aa_1(0):Aa_2(0):Aa_3(0)' | wc -l

cat $AnnotTab | grep 'HET' | wc -l
cat $OutDir/648_annotab_tenuissima_expanded.tsv | grep 'HET' | wc -l
#cat $AnnotTab | grep 'HET' | grep -v 'IPR010730' | less -S

cat $OutDir/648_annotab_tenuissima_expanded.tsv  | grep 'CAZY' | wc -l

cat $OutDir/648_annotab_tenuissima_expanded.tsv | grep -e 'Transcription factor' -e 'Zinc finger'  | grep 'Aa_1(0):Aa_2(0):Aa_3(0)' | wc -l
```





The file was downloaded to my local computer and used for analysis in R

```r
library(scales)
summarised_orthogroups <- read.delim("~/Downloads/Aalt/enrichment/orthogroups/summarised_orthogroups.tsv")
df1 <- subset(summarised_orthogroups, Species == 'A.alternata_ssp._arborescens' | Species == 'A.alternata_ssp._tenuissima')
#df2 <- df1[, colSums(df1 != 0) > 0]
# Remove orthogroups with identical numbers
df2 <- df1[vapply(df1, function(x) length(unique(x)) > 1, logical(1L))]
# Convert columns to numeric
df2[3:ncol(df2)] <- lapply(df2[3:ncol(df2)], as.numeric)
# Drop the OrthomclID column
df2$OrthomclID <- NULL
# Take those columns where all samples in one species are greater and identical
# arb_greater <- function(x) {
#   min(x[df2$Species == 'A.alternata_ssp._arborescens'])
#   max(x[df2$Species == 'A.alternata_ssp._tenuissima'])
#   }


df3 <- df2[vapply(df2, function(x) length(unique(x[df2$Species == 'A.alternata_ssp._arborescens'])) == 1 && length(unique(x[df2$Species == 'A.alternata_ssp._tenuissima'])) == 1, logical(1L))]
df3$Species <- df2$Species

df4 <- df2[vapply(df2, function(x) length(unique(x[df2$Species == 'A.alternata_ssp._arborescens'])) > 1 | length(unique(x[df2$Species == 'A.alternata_ssp._tenuissima'])) > 1, logical(1L))]
# df4 <- cbind(df2$Species, df4)

# df4[3:ncol(df2)] <- lapply(df2[3:ncol(df2)], as.numeric)

t <- data.frame(t(sapply(df4, function(x)
      unlist(t.test(x~df2$Species)[c("estimate","p.value","statistic","conf.int")]))))
df5 <- data.frame(colnames(df4))
df5$pval <- round(t$p.value, 3)
df5$bh_pval <- round(p.adjust(t$p.value, method = "hochberg"), 3)
df5$tstat <- round(t$statistic, 2)
```


# Numbers of genes, secreted proteins and effectors were tested:

```r
gene_numbers <- read.delim("~/Downloads/Aalt/enrichment/gene_numbers/gene_numbers.txt")
View(gene_numbers)
gene_numbers$Strain <- as.factor(gene_numbers$Strain)


library(reshape2)
library(ggplot2)
library(scales)
df1 <- gene_numbers[c("Organism", "Strain", "Total.genes", "Secreted.genes", "Secreted...effectorP", "Secreted.CAZYmes", "Secondary.metabolite.clusters")]
colnames(df1) <- c("Organism", "Strain", "Total genes", "Secreted genes", "Secreted effectorP", "Secreted CAZYmes", "Secondary metabolite clusters")
df1$Organism <- factor(df1$Organism, levels = c("pear pathotype", "apple pathotype", "tenuissima clade", "arborescens clade"))
df2 <- melt(df1, id.vars=c(1,2))
p <- ggplot(df2, aes(x=Organism, y=value)) +
  scale_y_continuous(name = "", breaks= pretty_breaks()) +
  geom_boxplot()
p <- p + theme(axis.text.x=element_text(angle = -45, hjust = 0))
p <- p + scale_x_discrete(name ="")
p <- p + facet_wrap(~df2$variable, nrow = 3, ncol = 2, strip.position = "top", scales="free_y")
p <- p + theme(plot.margin=unit(c(0,1,0,0),"cm"))
outfile='gene_enrichment.pdf'
ggsave(outfile , plot = p, device = 'pdf', width =15, height = 15, units = "cm", limitsize = FALSE)
outfile='gene_enrichment.tiff'
ggsave(outfile , plot = p, device = 'tiff', width =15, height = 15, units = "cm", limitsize = FALSE)

df3 <- subset(gene_numbers, Organism != "pear pathotype")

fit <- aov(Total.genes ~ Organism, data = df3)
summary(fit)
TukeyHSD(fit)

fit <- aov(Secreted.genes ~ Organism, data = df3)
summary(fit)
TukeyHSD(fit)

fit <- aov(Secreted...effectorP ~ Organism, data = df3)
summary(fit)
TukeyHSD(fit)

fit <- aov(Secreted.CAZYmes ~ Organism, data = df3)
summary(fit)
TukeyHSD(fit)

fit <- aov(Secondary.metabolite.clusters ~ Organism, data = df3)
summary(fit)
TukeyHSD(fit)
```
