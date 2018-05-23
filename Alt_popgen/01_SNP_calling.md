
# 1. Alignment of Pcac raw reads vs the 414 genome

Alignment of reads from a single run:

```bash
  Reference=$(ls repeat_masked/A.alternata_ssp_tenuissima/1166/filtered_contigs/1166_contigs_unmasked.fa)
  for StrainPath in $(ls -d ../../../../home/groups/harrisonlab/project_files/alternaria/qc_dna/paired/*/*); do
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/spades
    Strain=$(echo $StrainPath | rev | cut -f1 -d '/' | rev)
    Organism=$(echo $StrainPath | rev | cut -f2 -d '/' | rev)
    F_Read=$(ls $StrainPath/F/*_trim.fq.gz)
    R_Read=$(ls $StrainPath/R/*_trim.fq.gz)
    echo $F_Read
    echo $R_Read
    OutDir=analysis/genome_alignment/bowtie/$Organism/$Strain/vs_1166
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/genome_alignment
    qsub $ProgDir/bowtie/sub_bowtie.sh $Reference $F_Read $R_Read $OutDir
  done
```

# 2. Pre SNP calling cleanup


## 2.1 Rename input mapping files in each folder by prefixing with the strain ID

```bash
  for File in $(ls analysis/genome_alignment/bowtie/*/*/vs_1166/*_sorted.bam); do
    Strain=$(echo $File | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $File | rev | cut -f4 -d '/' | rev)
    echo "$Organism - $Strain"
    OutDir=analysis/popgen/$Organism/$Strain
    CurDir=$PWD
    mkdir -p $OutDir
    cd $OutDir
    cp -s $CurDir/$File "$Strain"_vs_1166_aligned_sorted.bam
    cd $CurDir
  done
```

## 2.2 Remove multimapping reads, discordant reads. PCR and optical duplicates, and add read group and sample name to each mapped read (preferably, the shortest ID possible)

Convention used:
qsub $ProgDir/sub_pre_snp_calling.sh <INPUT SAM FILE> <SAMPLE_ID>

```bash
for Sam in $(ls analysis/popgen/*/*/*_vs_1166_aligned_sorted.bam); do
Strain=$(echo $Sam | rev | cut -f2 -d '/' | rev)
Organism=$(echo $Sam | rev | cut -f3 -d '/' | rev)
echo "$Organism - $Strain"
ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/Pcac_popgen
qsub $ProgDir/sub_pre_snp_calling.sh $Sam $Strain
done
```

# 3. Run SNP calling

#Runs a SNP calling script from Maria in order to be able to draw up a phylogeny
To change in each analysis:

<!-- ```bash
input=/home/groups/harrisonlab/project_files/phytophthora_fragariae/analysis/genome_alignment/bowtie
reference=repeat_masked/P.fragariae/Bc16/filtered_contigs_repmask/95m_contigs_unmasked.fa

filename=$(basename "$reference")
output="${filename%.*}.dict"
``` -->

##Prepare genome reference indexes required by GATK

```bash
Reference=$(ls repeat_masked/A.alternata_ssp_tenuissima/1166/filtered_contigs/1166_contigs_unmasked.fa)
OutName=$(echo $Reference | sed 's/.fa/.dict/g')
OutDir=$(dirname $Reference)
mkdir -p $OutDir
ProgDir=/home/sobczm/bin/picard-tools-2.5.0
java -jar $ProgDir/picard.jar CreateSequenceDictionary R=$Reference O=$OutName
samtools faidx $Reference
```

###Copy index file to same folder as BAM alignments
<!--
```bash
Reference=$(ls repeat_masked/P.cactorum/414_v2/filtered_contigs_repmask/414_v2_contigs_unmasked.fa)
for AlignDir in $(ls -d analysis/popgen/P.*/*/); do
    Index="$Reference".dict
    Directory=analysis/genome_alignment/bowtie/*/$Strain/vs_414/
    cp $Index $AlignDir/.
done
``` -->

Move to the directory where the output of SNP calling should be placed. Then
Start SNP calling with GATK.
The submission script required need to be custom-prepared for each analysis,
depending on what samples are being analysed. See inside the submission script
below:

```bash
CurDir=$PWD
OutDir=analysis/popgen/SNP_calling
mkdir -p $OutDir
cd $OutDir
ProgDir=/home/armita/git_repos/emr_repos/scripts/alternaria/Alt_popgen
qsub $ProgDir/sub_SNP_calling_multithreaded.sh
cd $CurDir
```

## Filter SNPs based on this region being present in all isolates

Only retain biallelic high-quality SNPS with no missing data (for any individual) for genetic analyses below (in some cases, may allow some missing data in order to retain more SNPs, or first remove poorly sequenced individuals with too much missing data and then filter the SNPs).

```bash
# cp analysis/popgen/SNP_calling/414_contigs_unmasked_temp.vcf analysis/popgen/SNP_calling/414_contigs_unmasked.vcf
Vcf=$(ls analysis/popgen/SNP_calling/1166_contigs_unmasked.vcf)
ProgDir=/home/armita/git_repos/emr_repos/scripts/popgen/snp
# mq=40
# qual=30
# dp=10
# gq=30
# na=0.95
# removeindel=Y
# $VcfLib/vcffilter -f "QUAL > $qual & MQ > $mq"
# $vcftools/vcftools --vcf temp.vcf --max-missing $na --remove-indels --recode --out ${filename%.vcf}_filtered
qsub $ProgDir/sub_vcf_parser.sh $Vcf 40 30 10 30 1 Y
```

```bash
mv 1166_contigs_unmasked_filtered.vcf analysis/popgen/SNP_calling/1166_contigs_unmasked_filtered.vcf
```


## Remove sequencing errors from vcf files:

```bash
Vcf=$(ls analysis/popgen/SNP_calling/1166_contigs_unmasked_filtered.vcf)
OutDir=$(dirname $Vcf)
Errors=$OutDir/1166_error_SNPs.tsv
FilteredVcf=$OutDir/1166_contigs_unmasked_filtered_no_errors.vcf
ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/Pcac_popgen
$ProgDir/flag_error_SNPs.py --ploidy 'haploid' --inp_vcf $Vcf --ref_isolate 1166 --errors $Errors --filtered $FilteredVcf
echo "The number of probable errors from homozygous SNPs being called from reference illumina reads vs the reference assembly is:"
cat $Errors
# cat $Errors | wc -l
echo "These have been removed from the vcf file"
```

```
contig_10	1837658
contig_10	1837760
```

<!--
In some organisms, may want to thin (subsample) SNPs in high linkage diseqilibrium down to
1 SNP  per e.g. 10 kbp just for the population structure analyses.
```bash
VcfTools=/home/sobczm/bin/vcftools/bin
$VcfTools/vcftools --vcf $input_vcf --thin 10000 --recode --out ${input_vcf%.vcf}_thinned
```
-->

## Collect VCF stats

General VCF stats (remember that vcftools needs to have the PERL library exported)

```bash
  VcfTools=/home/sobczm/bin/vcftools/bin
  export PERL5LIB="$VcfTools:$PERL5LIB"
  # Vcf=$(ls analysis/popgen/SNP_calling/*_filtered_no_errors.vcf)
  # Stats=$(echo $Vcf | sed 's/.vcf/.stat/g')
  # perl $VcfTools/vcf-stats $Vcf > $Stats
  VcfFiltered=$(ls analysis/popgen/SNP_calling/1166_contigs_unmasked_filtered_no_errors.vcf)
  Stats=$(echo $VcfFiltered | sed 's/.vcf/.stat/g')
  perl $VcfTools/vcf-stats $VcfFiltered > $Stats
```

Calculate the index for percentage of shared SNP alleles between the individuals.

```bash
  for Vcf in $(ls analysis/popgen/SNP_calling/*filtered_no_errors.vcf); do
      ProgDir=/home/armita/git_repos/emr_repos/scripts/popgen/snp
      $ProgDir/similarity_percentage.py $Vcf
  done
```

# Visualise the output as heatmap and clustering dendrogram
```bash
for Log in $(ls analysis/popgen/SNP_calling/*distance.log); do
  ProgDir=/home/armita/git_repos/emr_repos/scripts/popgen/snp
  Rscript --vanilla $ProgDir/distance_matrix.R $Log
  mv Rplots.pdf analysis/popgen/SNP_calling/.
done
```


## Carry out PCA and plot the results

```bash
for Vcf in $(ls analysis/popgen/SNP_calling/*filtered_no_errors.vcf); do
    echo $Vcf
    ProgDir=/home/armita/git_repos/emr_repos/scripts/popgen/snp
    # Out=$(basename $Vcf)
    Out=analysis/popgen/SNP_calling
    echo $Out
    /home/deakig/R3.4/bin/Rscript --vanilla $ProgDir/pca.R $Vcf $Out/PCA.pdf
done
```

```
The total number of samples: 12
The total number of SNPs: 48000
```


## Calculate a NJ tree

These commands didnt work as P. idaei is too distant for sufficient sites to be shared
between isolates

based on all the SNPs. Outputs a basic display of the tree, plus a Newick file to be used for displaying the tree in FigTree and beautifying it.

Remove all missing data for nj tree construction

```bash
  for Vcf in $(ls analysis/popgen/SNP_calling/*_filtered_no_errors.vcf); do
    echo $Vcf
    Out=$(basename $Vcf .vcf)
    echo $Out
    VcfTools=/home/sobczm/bin/vcftools/bin
    $VcfTools/vcftools --vcf $Vcf --mac 1 --max-missing 1.0 --recode --out analysis/popgen/SNP_calling/"$Out"_no_missing
  done
```

```
After filtering, kept 12 out of 12 Individuals
Outputting VCF file...
After filtering, kept 48414 out of a possible 48414 Sites
```


```bash
for Vcf in $(ls analysis/popgen/SNP_calling/*_no_missing.recode.vcf); do
    echo $Vcf
    Ploidy=1
    ProgDir=/home/armita/git_repos/emr_repos/scripts/popgen/snp
    $ProgDir/nj_tree.sh $Vcf $Ploidy
    mv Rplots.pdf analysis/popgen/SNP_calling/NJ_tree.pdf
done
```



# Identify SNPs in gene models:

Create custom SnpEff genome database

```bash
SnpEff=/home/sobczm/bin/snpEff
nano $SnpEff/snpEff.config
```


Add the following lines to the section with databases:

```
#---
# EMR Databases
#----
# Fus2 genome
Fus2v1.0.genome : Fus2
# Bc16 genome
Bc16v1.0.genome: BC-16
# P414 genome
P414v1.0.genome: 414
# Alt1166 genome
1166v1.0.genome: 1166
# Alt650 genome
650v1.0.genome: 650
```

Collect input files

```bash
Reference=$(ls repeat_masked/A.alternata_ssp_tenuissima/1166/filtered_contigs/1166_contigs_unmasked.fa)
Gff=$(ls gene_pred/final/A.alternata_ssp_tenuissima/1166/final/final_genes_appended_renamed.gff3)
SnpEff=/home/sobczm/bin/snpEff
mkdir $SnpEff/data/1166v1.0
cp $Reference $SnpEff/data/1166v1.0/sequences.fa
cp $Gff $SnpEff/data/1166v1.0/genes.gff

#Build database using GFF3 annotation
java -jar $SnpEff/snpEff.jar build -gff3 -v 1166v1.0
```


## Annotate VCF files
```bash
CurDir=/data/scratch/armita/alternaria
cd $CurDir
for a in $(ls analysis/popgen/SNP_calling/*_filtered_no_errors.vcf); do
    echo $a
    filename=$(basename "$a")
    Prefix=${filename%.vcf}
    OutDir=$(ls -d analysis/popgen/SNP_calling)
    SnpEff=/home/sobczm/bin/snpEff
    java -Xmx4g -jar $SnpEff/snpEff.jar -v -ud 0 1166v1.0 $a > $OutDir/"$Prefix"_annotated.vcf
    mv snpEff_genes.txt $OutDir/snpEff_genes_"$Prefix".txt
    mv snpEff_summary.html $OutDir/snpEff_summary_"$Prefix".html
    # mv 414_v2_contigs_unmasked_filtered* $OutDir/.
    #-
    #Create subsamples of SNPs containing those in a given category
    #-
    #genic (includes 5', 3' UTRs)
    java -jar $SnpEff/SnpSift.jar filter "(ANN[*].EFFECT has 'missense_variant') || (ANN[*].EFFECT has 'nonsense_variant') || (ANN[*].EFFECT has 'synonymous_variant') || (ANN[*].EFFECT has 'intron_variant') || (ANN[*].EFFECT has '5_prime_UTR_variant') || (ANN[*].EFFECT has '3_prime_UTR_variant')" $OutDir/"$Prefix"_annotated.vcf > $OutDir/"$Prefix"_gene.vcf
    #coding
    java -jar $SnpEff/SnpSift.jar filter "(ANN[0].EFFECT has 'missense_variant') || (ANN[0].EFFECT has 'nonsense_variant') || (ANN[0].EFFECT has 'synonymous_variant')" $OutDir/"$Prefix"_annotated.vcf > $OutDir/"$Prefix"_coding.vcf
    #non-synonymous
    java -jar $SnpEff/SnpSift.jar filter "(ANN[0].EFFECT has 'missense_variant') || (ANN[0].EFFECT has 'nonsense_variant')" $OutDir/"$Prefix"_annotated.vcf > $OutDir/"$Prefix"_nonsyn.vcf
    #synonymous
    java -jar $SnpEff/SnpSift.jar filter "(ANN[0].EFFECT has 'synonymous_variant')" $OutDir/"$Prefix"_annotated.vcf > $OutDir/"$Prefix"_syn.vcf
    #Four-fold degenrate sites (output file suffix: 4fd)
    ProgDir=/home/sobczm/bin/popgen/summary_stats
    python $ProgDir/parse_snpeff_synonymous.py $OutDir/"$Prefix"_syn.vcf
    AllSnps=$(cat $OutDir/"$Prefix"_annotated.vcf | grep -v '#' | wc -l)
    GeneSnps=$(cat $OutDir/"$Prefix"_gene.vcf | grep -v '#' | wc -l)
    CdsSnps=$(cat $OutDir/"$Prefix"_coding.vcf | grep -v '#' | wc -l)
    NonsynSnps=$(cat $OutDir/"$Prefix"_nonsyn.vcf | grep -v '#' | wc -l)
    SynSnps=$(cat $OutDir/"$Prefix"_syn.vcf | grep -v '#' | wc -l)
    printf "$AllSnps\t$GeneSnps\t$CdsSnps\t$NonsynSnps\t$SynSnps\n"
done
```

```
48414	32779	29777	6784	22993
```

```bash
Vcf=$(ls analysis/popgen/SNP_calling/1166_contigs_unmasked_filtered_no_errors_annotated.vcf)
for i in $(seq 10 21); do
Strain=$(cat $Vcf| grep -v '##' | grep '#' | cut -f$i)
SNPs=$(cat $Vcf| grep -v '##' | cut -f$i | grep '^1:' | wc -l)
echo "$Strain - $SNPs"
done
```
Total SNPs
```
1166 - 0
635 - 4301
743 - 4299
1177 - 6151

648 - 6165
1082 - 6652
1164 - 5170
24350 - 6530

675 - 20016
97.0013 - 20030
97.0016 - 20126

650 - 23278
```

Gene SNPs
```bash
Vcf=$(ls analysis/popgen/SNP_calling/1166_contigs_unmasked_filtered_no_errors_gene.vcf)
for i in $(seq 10 21); do
Strain=$(cat $Vcf| grep -v '##' | grep '#' | cut -f$i)
SNPs=$(cat $Vcf| grep -v '##' | cut -f$i | grep '^1:' | wc -l)
Genes=$(cat $Vcf| grep -v '##' | cut -f8,$i | grep -P "\t1:" | cut -f4 -d '|' | sort | uniq | wc -l)
printf "$Strain\t$SNPs\t$Genes\n"
done
```

```
1082	4753	2185
1164	3675	1722
1166	0	0
1177	4349	1958
24350	4650	2108
635	3072	1438
648	4381	2038
650	16010	4938
675	13555	4726
743	3070	1438
97.0013	13583	4718
97.0016	13614	4712
```

Non-syn SNPs
```bash
AnnotTab=$(ls gene_pred/annotation/A.alternata_ssp_tenuissima/1166/1166_annotation_ncbi.tsv)
OutDir=$(dirname $AnnotTab)

# SigPHeaders=$(ls gene_pred/braker_signalp-4.1/A.alternata_ssp_tenuissima/1166/1166_final_sp_no_trans_mem_headers.txt)
SigPHeaders=$OutDir/1166_secreted.txt
cat $AnnotTab | cut -f1,8 | grep 'Yes' | cut -f1 > $SigPHeaders
EffPHeaders=$OutDir/1166_effectorP.txt
cat $AnnotTab | cut -f1,8,9 | grep "Yes.Yes" | cut -f1 > $EffPHeaders
CAZymeHeaders=$OutDir/1166_CAZymes.txt
cat $AnnotTab | cut -f1,8,10 | grep 'CAZY' | grep 'Yes' | cut -f1 > $CAZymeHeaders
SecMetHeaders=$OutDir/1166_SecMet.txt
cat $AnnotTab | grep 'AS_' | cut -f1 > $SecMetHeaders


Vcf=$(ls analysis/popgen/SNP_calling/1166_contigs_unmasked_filtered_no_errors_nonsyn.vcf)
for i in $(seq 10 21); do
Strain=$(cat $Vcf| grep -v '##' | grep '#' | cut -f$i)
SNPs=$(cat $Vcf| grep -v '##' | cut -f$i | grep '^1:' | wc -l)
Genes=$(cat $Vcf| grep -v '##' | cut -f8,$i | grep -P "\t1:" | cut -f7 -d '|' | sort | uniq | wc -l)

Secreted=$(cat $Vcf | grep -v '##' | cut -f8,$i | grep -P "\t1:" | cut -f7 -d '|' | sort | uniq | grep -f $SigPHeaders | wc -l)
CAZYmes=$(cat $Vcf | grep -v '##' | cut -f8,$i | grep -P "\t1:" | cut -f7 -d '|' | sort | uniq | grep -f $CAZymeHeaders | wc -l)
Effectors=$(cat $Vcf | grep -v '##' | cut -f8,$i | grep -P "\t1:" | cut -f7 -d '|' | sort | uniq | grep -f $EffPHeaders | wc -l)
Toxins=$(cat $Vcf | grep -v '##' | cut -f8,$i | grep -P "\t1:" | cut -f7 -d '|' | sort | uniq | grep -f $SecMetHeaders | wc -l)
printf "$Strain\t$SNPs\t$Genes\t$Secreted\t$CAZYmes\t$Effectors\t$Toxins\n"
done
```

```
1082	953	718	44	21	2	30
1164	702	566	35	14	3	25
1166	0	0	0	0	0	0
1177	856	650	48	19	2	41
24350	926	705	44	22	3	45
635	596	478	22	10	2	21
648	856	666	37	15	1	35
650	3131	1880	133	48	7	95
675	2555	1723	121	54	11	89
743	595	478	22	10	2	21
97.0013	2569	1753	128	61	12	89
97.0016	2602	1752	125	59	12	92
```

## SNP distribution

```bash
Vcf=$(ls analysis/popgen/SNP_calling/1166_contigs_unmasked_filtered_no_errors_annotated.vcf)
cat $Vcf | grep -v '#' | cut -f1 | uniq -c
```

```
5563 contig_1
4544 contig_2
4543 contig_3
4900 contig_4
4345 contig_5
4778 contig_6
3800 contig_7
3670 contig_8
3246 contig_9
2786 contig_10
3328 contig_11
2154 contig_12
 329 contig_13
 208 contig_16
 206 contig_17
  14 contig_22
```

```bash
Vcf=$(ls analysis/popgen/SNP_calling/1166_contigs_unmasked_filtered_no_errors_syn.vcf)
cat $Vcf | grep -v '#' | cut -f1 | uniq -c
```

```
2817 contig_1
2228 contig_2
1840 contig_3
2330 contig_4
2073 contig_5
2155 contig_6
1811 contig_7
1777 contig_8
1441 contig_9
1450 contig_10
1603 contig_11
 978 contig_12
 229 contig_13
 110 contig_16
 137 contig_17
  14 contig_22
```


```bash
Vcf=$(ls analysis/popgen/SNP_calling/1166_contigs_unmasked_filtered_no_errors_nonsyn.vcf)
cat $Vcf | grep -v '#' | cut -f1 | uniq -c
```
```
840 contig_1
598 contig_2
647 contig_3
733 contig_4
697 contig_5
704 contig_6
526 contig_7
481 contig_8
404 contig_9
325 contig_10
423 contig_11
293 contig_12
 42 contig_13
 42 contig_16
 29 contig_17
```


# 4.2 arborescens isolates in comparison to 1166

```bash
  Prefix=arborescens_vs_1166
  OutDir=analysis/popgen/SNP_calling/$Prefix
  mkdir -p $OutDir

  Vcf=$(ls analysis/popgen/SNP_calling/*_filtered_no_errors.vcf)
  ExcludeList="1082 1164 1166 1177 24350 635 648 650 743"
  VcfLib=/home/sobczm/bin/vcflib/bin
  $VcfLib/vcfremovesamples $Vcf $ExcludeList > $OutDir/$Prefix.vcf

  VcfTools=/home/sobczm/bin/vcftools/bin
  $VcfTools/vcftools --vcf $OutDir/$Prefix.vcf --max-missing 0.95 --remove-indels --mac 1 --recode --out $OutDir/"$Prefix"_filtered_no_indels

  for Vcf in $(ls $OutDir/"$Prefix"_filtered_no_indels.recode.vcf); do
      echo $Vcf
      # ProgDir=/home/armita/git_repos/emr_repos/scripts/popgen/summary_stats
      # $ProgDir/annotate_snps_genome.sh $Vcf P414v1.0

      filename=$(basename "$Vcf")
      Prefix=$(echo $filename | sed 's/.vcf//g')
      SnpEff=/home/sobczm/bin/snpEff
      java -Xmx4g -jar $SnpEff/snpEff.jar -v -ud 0 1166v1.0 $Vcf > $OutDir/"$Prefix"_annotated.vcf
      mv snpEff_genes.txt $OutDir/snpEff_genes_$Prefix.txt
      mv snpEff_summary.html $OutDir/snpEff_summary_$Prefix.html

      #Create subsamples of SNPs containing those in a given category

      #genic (includes 5', 3' UTRs)
      java -jar $SnpEff/SnpSift.jar filter "(ANN[*].EFFECT has 'missense_variant') || (ANN[*].EFFECT has 'nonsense_variant') || (ANN[*].EFFECT has 'synonymous_variant') || (ANN[*].EFFECT has 'intron_variant') || (ANN[*].EFFECT has '5_prime_UTR_variant') || (ANN[*].EFFECT has '3_prime_UTR_variant')" $OutDir/"$Prefix"_annotated.vcf > $OutDir/"$Prefix"_gene.vcf
      #coding
      java -jar $SnpEff/SnpSift.jar filter "(ANN[0].EFFECT has 'missense_variant') || (ANN[0].EFFECT has 'nonsense_variant') || (ANN[0].EFFECT has 'synonymous_variant')" $OutDir/${filename%.vcf}_annotated.vcf > $OutDir/"$Prefix"_coding.vcf
      #non-synonymous
      java -jar $SnpEff/SnpSift.jar filter "(ANN[0].EFFECT has 'missense_variant') || (ANN[0].EFFECT has 'nonsense_variant')" $OutDir/"$Prefix"_annotated.vcf > $OutDir/"$Prefix"_nonsyn.vcf
      #synonymous
      java -jar $SnpEff/SnpSift.jar filter "(ANN[0].EFFECT has 'synonymous_variant')" $OutDir/"$Prefix"_annotated.vcf > $OutDir/"$Prefix"_syn.vcf
      #Four-fold degenrate sites (output file suffix: 4fd)
      ProgDir=/home/sobczm/bin/popgen/summary_stats
      python $ProgDir/parse_snpeff_synonymous.py $OutDir/"$Prefix"_syn.vcf

      # Identify SNP frequency within each group:
      $VcfTools/vcftools --vcf $OutDir/"$Prefix"_syn.vcf --freq --out $OutDir/"$Prefix"_syn
      $VcfTools/vcftools --vcf $OutDir/"$Prefix"_nonsyn.vcf --freq --out $OutDir/"$Prefix"_nonsyn
  done
```

```bash
cat analysis/popgen/SNP_calling/arborescens_vs_1166/arborescens_vs_1166_filtered_no_indels.recode_nonsyn.vcf | grep -v '#' | wc -l
cat analysis/popgen/SNP_calling/arborescens_vs_1166/arborescens_vs_1166_filtered_no_indels.recode_nonsyn.vcf  | grep -v '#' | cut -f8 | cut -f7 -d '|' | sort | uniq | wc -l
cat analysis/popgen/SNP_calling/arborescens_vs_1166/arborescens_vs_1166_filtered_no_indels.recode_nonsyn.vcf | grep -v '#' | cut -f8 | cut -f7 -d '|' | sort | uniq > tmp.txt
AnnotTab=$(ls gene_pred/annotation/A.alternata_ssp_tenuissima/1166/1166_annotation_ncbi.tsv)
OutDir=$(dirname $AnnotTab)
for Gene in $(cat tmp.txt); do
cat $AnnotTab | grep -P "^$Gene"
done > gene_pred/annotation/A.alternata_ssp_tenuissima/1166/A.arborescens_non-syn_SNP_genes.tsv
```


```bash
echo "Number of genes per contig"
cat $AnnotTab | cut -f2 | uniq -c
echo "Number of non-syn SNPs per contig"
cat gene_pred/annotation/A.alternata_ssp_tenuissima/1166/A.arborescens_non-syn_SNP_genes.tsv | cut -f2 | sort | uniq -c | sort -nr
```

```
1569 contig_1
1419 contig_2
1164 contig_3
1113 contig_4
1060 contig_5
 991 contig_6
 961 contig_7
 964 contig_8
 933 contig_9
 874 contig_10
 825 contig_11
 507 contig_12
 283 contig_13
 189 contig_14
 169 contig_15
 165 contig_16
 169 contig_17
 110 contig_18
  71 contig_19
  35 contig_20
  50 contig_21
  42 contig_22
```

```
91 contig_6
70 contig_1
68 contig_4
64 contig_3
60 contig_11
54 contig_2
45 contig_7
45 contig_5
43 contig_8
39 contig_9
34 contig_10
23 contig_12
 6 contig_16
 6 contig_13
 5 contig_17
 ```

 ```bash

cat gene_pred/annotation/A.alternata_ssp_tenuissima/1166/A.arborescens_non-syn_SNP_genes.tsv | grep 'contig_6'
 ```

### Number of secreted genes

```bash
NonSynTable=$(ls gene_pred/annotation/A.alternata_ssp_tenuissima/1166/A.arborescens_non-syn_SNP_genes.tsv)
echo "Number of secreted genes"
cat $AnnotTab | cut -f8 | grep 'Yes' | wc -l
echo "Number of non-syn secreted genes"
cat $NonSynTable | cut -f8 | grep 'Yes' | wc -l
echo "Number of EffP genes"
cat $AnnotTab | cut -f8,9 | grep 'Yes.Yes' | wc -l
echo "Number of non-syn EffP genes"
cat $NonSynTable | cut -f8,9 | grep 'Yes.Yes' | wc -l
echo "Number of CAZY genes"
cat $AnnotTab | cut -f8,10 | grep 'Yes' | grep 'CAZY' | wc -l
echo "Number of non-syn CAZY genes"
cat $NonSynTable | cut -f8,10 | grep 'Yes' | grep 'CAZY' | wc -l
echo "Number of Antismash genes"
cat $AnnotTab | grep 'AS_' | wc -l
echo "Number of non-syn antismash genes"
cat $NonSynTable | grep 'AS_' | wc -l
```

```bash
cat $NonSynTable | grep 'AS_' | cut -f12 | sort | uniq -c | sort -nr
```
```
     12 Cluster_13:AS_13:nrps
      4 Cluster_18:AS_18:t1pks
      3 Cluster_7:AS_7:t1pks
      3 Cluster_5:AS_5:nrps
      3 Cluster_15:AS_15:other
      2 Cluster_8:AS_8:terpene
      2 Cluster_6:AS_6:t1pks
      2 Cluster_10:AS_10:terpene
      1 Cluster_26:AS_26:t1pks-nrps
      1 Cluster_20:AS_20:other
      1 Cluster_19:AS_19:t1pks
      1 Cluster_16:AS_16:other
      1 Cluster_12:AS_12:t1pks-nrps
```
