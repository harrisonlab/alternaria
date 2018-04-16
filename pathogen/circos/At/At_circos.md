 <!-- ```bash
  OutDir=analysis/circos/At_circos
  mkdir -p $OutDir
  ProgDir=/home/armita/git_repos/emr_repos/scripts/alternaria/pathogen/circos
  Assembly=$(ls repeat_masked/*/*/filtered_contigs/*_contigs_softmasked_repeatmasker_TPSI_appended.fa | grep '1166')
  $ProgDir/fasta2circos.py --genome $Assembly --contig_prefix "" > $OutDir/Assembly.txt
  for GffFile in $(ls gene_pred/final/*/*/final/final_genes_appended.gff3 | grep '1166'); do
    echo $GffFile
    Chr=$(echo $GffFile | rev |cut -f1 -d'/' | rev | cut -f6 -d '_')
    $ProgDir/gff2circos_scatterplot.py --gff $GffFile --feature Chr"$Chr"_gene_homolog > $OutDir/FoL_chr"$Chr"_genes.txt
  done
  circos -conf /home/armita/git_repos/emr_repos/scripts/alternaria/pathogen/circos/At/At_circos.conf -outputdir ./$OutDir
```


```bash
  ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
  Assembly=$(ls repeat_masked/*/*/filtered_contigs/*_contigs_softmasked_repeatmasker_TPSI_appended.fa | grep '1166')
  $ProgDir/fasta2circos.py --genome $Assembly --contig_prefix "" > tmp5/Assembly.txt
  for GffFile in $(ls analysis/blast_homology/*/Fus2_pacbio_merged_richards/*_chr_*_gene_single_copy.aa_hits.gff); do
    echo $GffFile
    Chr=$(echo $GffFile | rev |cut -f1 -d'/' | rev | cut -f6 -d '_')
    $ProgDir/gff2circos_scatterplot.py --gff $GffFile --feature Chr"$Chr"_gene_homolog > tmp5/FoL_chr"$Chr"_genes.txt
  done
  circos -conf /home/armita/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos/Fus2/Fus2_circos.conf -outputdir ./tmp5
``` -->

```bash
qlogin
cd /data/scratch/armita/alternaria
  OutDir=analysis/circos/At_circos
for ReadsBam in $(ls analysis/genome_alignment/bowtie/*/*/vs_1166/1166_contigs_unmasked.fa_aligned_sorted.bam | grep -e '1166' -e '648' -e '650' -e '97.0016' | grep -w -v -e '675' -e '97.0013' -e '97.0016' -e '650' -e '1082' -e '1164' -e '1166'); do
  Strain=$(echo $ReadsBam | rev | cut -f3 -d '/' | rev)
  Organism=$(echo $ReadsBam | rev | cut -f4 -d '/' | rev)
  # AlignDir=$(dirname $ReadsBam)
  echo "$Organism - $Strain"
  bedtools coverage -a $OutDir/100kb_windows.gff -b $ReadsBam > $OutDir/"$Strain"_coverage.bed
  # Convert coverage bed files into circos format
  ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
  $ProgDir/coverage_bed2circos.py --bed $OutDir/"$Strain"_coverage.bed > $OutDir/"$Strain"_coverage_scatterplot.txt
done
```

```bash
  # ProgDir=/home/armita/git_repos/emr_repos/scripts/alternaria/pathogen/circos
  OutDir=analysis/circos/At_circos
  mkdir -p $OutDir
  ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
  Assembly=$(ls repeat_masked/*/*/filtered_contigs/*_contigs_softmasked_repeatmasker_TPSI_appended.fa | grep '1166')
  # Convert the Fus2 genome into circos format
  $ProgDir/fasta2circos.py --genome $Assembly --contig_prefix "" > $OutDir/Assembly.txt

  # Make 100kb windows for plots
  $ProgDir/fasta2gff_windows.py --genome $Assembly > $OutDir/100kb_windows.gff

  # Identify GC content in 100kb windows
  $ProgDir/gc_content2circos.py --genome $Assembly --gff $OutDir/100kb_windows.gff > $OutDir/GC_scatterplot.txt

  ls analysis/genome_alignment/bowtie/*/*/vs_1166/*_depth_10kb.tsv

  # Plot location of Fus2 mimps and secreted effectorP genes as a scatterplot
  GffCAZY=$(ls gene_pred/CAZY/*/1166/*_CAZY_secreted.gff)
  ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
  $ProgDir/gff2circos_scatterplot.py --gff $GffCAZY --feature gene --value 1 > $OutDir/CAZY_plot.txt
  GffEffP=$(ls analysis/effectorP/*/1166/*_EffectorP_secreted.gff)
  ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
  $ProgDir/gff2circos_scatterplot.py --gff $GffEffP --feature gene --value 1 > $OutDir/effectorP_plot.txt
  GffAntiSmash=$(ls gene_pred/secondary_metabolites/antismash/*/1166/*_secondary_metabolite_regions.gff_secmet_clusters.gff)
  ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
  $ProgDir/gff2circos_scatterplot.py --gff $GffAntiSmash --feature indole indole-nrps nrps nrps-t1pks other t1pks t1pks-nrps t3pks terpene --value 1 > $OutDir/antismash_plot.txt

  BlastHits=$(ls analysis/blast_homology/*/1166/*_CDC_genes.fa_homologs.gff)
  # GffSIX=$OutDir/CDC_genes.gff
  # cat $BlastHits | grep -v -e 'MIMP' -e 'C5' -e 'CRX' > $GffSix
  ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
  $ProgDir/gff2circos_scatterplot.py --gff $BlastHits --feature toxin_homolog --value 1 > $OutDir/CDC_genes_plot.txt

  circos -conf /home/armita/git_repos/emr_repos/scripts/alternaria/pathogen/circos/At/At_circos.conf -outputdir ./$OutDir
```
