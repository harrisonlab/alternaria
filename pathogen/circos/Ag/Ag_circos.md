
```bash
qlogin
cd /data/scratch/armita/alternaria
OutDir=analysis/circos/Ag_circos
for ReadsBam in $(ls analysis/genome_alignment/bowtie/*/*/vs_650/650_contigs_unmasked.fa_aligned_sorted.bam | grep -w -e '1166' -e '648' -e '650' -e '97.0016' | grep -e '648' -e '1166'); do
 Strain=$(echo $ReadsBam | rev | cut -f3 -d '/' | rev)
 Organism=$(echo $ReadsBam | rev | cut -f4 -d '/' | rev)
 echo "$Organism - $Strain"
 bedtools coverage -a $OutDir/100kb_windows.gff -b $ReadsBam > $OutDir/"$Strain"_coverage.bed
 # Convert coverage bed files into circos format
 ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
 $ProgDir/coverage_bed2circos.py --bed $OutDir/"$Strain"_coverage.bed > $OutDir/"$Strain"_coverage_scatterplot.txt
done
```

```bash
 # ProgDir=/home/armita/git_repos/emr_repos/scripts/alternaria/pathogen/circos
 OutDir=analysis/circos/Ag_circos
 mkdir -p $OutDir
 ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
 Assembly=$(ls repeat_masked/*/*/filtered_contigs/*_contigs_softmasked_repeatmasker_TPSI_appended.fa | grep '650')
 # Convert the Fus2 genome into circos format
 $ProgDir/fasta2circos.py --genome $Assembly --contig_prefix "" > $OutDir/Assembly.txt

 # Make 100kb windows for plots
 $ProgDir/fasta2gff_windows.py --genome $Assembly > $OutDir/100kb_windows.gff

 # Identify GC content in 100kb windows
 $ProgDir/gc_content2circos.py --genome $Assembly --gff $OutDir/100kb_windows.gff > $OutDir/GC_scatterplot.txt

 ls analysis/genome_alignment/bowtie/*/*/vs_650/*_depth_10kb.tsv

 # Plot location of Fus2 mimps and secreted effectorP genes as a scatterplot
 GffCAZY=$(ls gene_pred/CAZY/*/650/*_CAZY_secreted.gff)
 ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
 $ProgDir/gff2circos_scatterplot.py --gff $GffCAZY --feature gene --value 1 > $OutDir/CAZY_plot.txt
 GffEffP=$(ls analysis/effectorP/*/650/*_EffectorP_secreted.gff)
 ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
 $ProgDir/gff2circos_scatterplot.py --gff $GffEffP --feature gene --value 1 > $OutDir/effectorP_plot.txt
  GffAntiSmash=$(ls gene_pred/secondary_metabolites/antismash/*/650/*_secondary_metabolite_regions.gff_secmet_clusters.gff)
 ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
 $ProgDir/gff2circos_scatterplot.py --gff $GffAntiSmash --feature indole indole-nrps nrps nrps-t1pks other t1pks t1pks-nrps t3pks terpene --value 1 > $OutDir/antismash_plot.txt

 BlastHits=$(ls analysis/blast_homology/*/650/*_CDC_genes.fa_homologs_AKT.gff)
 # GffSIX=$OutDir/CDC_genes.gff
 # cat $BlastHits | grep -v -e 'MIMP' -e 'C5' -e 'CRX' > $GffSix
 ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
 $ProgDir/gff2circos_scatterplot.py --gff $BlastHits --feature toxin_homolog --value 1 > $OutDir/CDC_genes_plot.txt

 circos -conf /home/armita/git_repos/emr_repos/scripts/alternaria/pathogen/circos/Ag/Ag_circos.conf -outputdir ./$OutDir
```
