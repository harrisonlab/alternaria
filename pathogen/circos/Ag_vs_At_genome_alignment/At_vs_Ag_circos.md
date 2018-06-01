```bash
  OutDir=analysis/circos/At_vs_Ag_genome_alignment_circos
  mkdir -p $OutDir

  Ag_genome=$(ls repeat_masked/*/*/filtered_contigs/*_contigs_softmasked_repeatmasker_TPSI_appended.fa | grep '650')
  ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
  $ProgDir/fasta2circos.py --genome $Ag_genome --contig_prefix "Ag_" > $OutDir/Ag_genome.txt

  At_genome=$(ls repeat_masked/*/*/filtered_contigs/*_contigs_softmasked_repeatmasker_TPSI_appended.fa | grep '1166')
  ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
  $ProgDir/fasta2circos.py --genome $At_genome --contig_prefix "At_" > $OutDir/At_genome.txt

  cat $OutDir/At_genome.txt > $OutDir/At_Ag_genome.txt
  tac $OutDir/Ag_genome.txt >> $OutDir/At_Ag_genome.txt

  # cat $OutDir/At_Ag_genome.txt | grep -v 'DS231' | grep -v -e 'chr23' -e 'chr24' -e 'chr25' -e 'chr26' -e 'chr27' -e 'chr28' -e 'chr29' -e 'chr30' -e 'chr31' -e 'chr32' -e 'chr33' -e 'chr34' > $OutDir/At_Ag_genome_edited.txt
```

```bash
OutDir=analysis/circos/At_vs_Ag_genome_alignment_circos
Coords=$(ls analysis/genome_alignment/mummer/A.alternata_ssp_tenuissima/1166/1166_vs_650/1166_vs_650_coords.tsv)
ProgDir=/home/armita/git_repos/emr_repos/scripts/alternaria/pathogen/genome_alignment
$ProgDir/nucmer_coords2circos.py --inp_coords $Coords --queery_id At --ref_id Ag > $OutDir/At_vs_Ag_links.txt
cat $OutDir/At_vs_Ag_links.txt > $OutDir/At_vs_Ag_links_edited.txt
```



A file showing contig orientations was made:
```bash
  cat $OutDir/At_Ag_genome.txt | cut -f3 -d ' ' | sed "s/$/; /g" | tr -d '\n' > $OutDir/At_contig_order.txt
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/circos
  $ProgDir/find_contig_orientation.py --links_file $OutDir/At_vs_Ag_links_edited.txt > $OutDir/At_vs_Ag_contig_orientation.txt
```
<!--
The number of bp in syntenous contigs was identified using:

```bash
  cat $OutDir/At_vs_Ag_contig_orientation.txt | tail -n3 | grep -v 'orientation' | sed 's/, /\n/g' > $OutDir/At_vs_Ag_syntenous_contigs.txt
  cat $OutDir/At_Ag_genome.txt | grep -f $OutDir/At_vs_Ag_syntenous_contigs.txt | cut -f6 -d ' ' | awk '{s+=$1} END {print s}'
``` -->

Contig order was selected by taking the first line of that file and then also
taking the reversed order of FoL contigs using the command:

```bash
cat $OutDir/At_vs_Ag_contig_orientation.txt | grep -A1 'Order of all seen contigs' | tail -n1 | sed "s/, /\n/g" > tmp.txt
cat $OutDir/At_vs_Ag_contig_orientation.txt | grep -A1 'Order of all seen contigs' | tail -n1
cat $OutDir/At_Ag_genome.txt | grep 'At' | grep -w -v -f tmp.txt | cut -f3 -d ' '| tr -d '\n' | sed 's/At/, At/g'
cat $OutDir/At_Ag_genome.txt | grep 'Ag' | cut -f3 -d ' ' | tr -d '\n' | sed 's/Ag/, Ag/g'
```

```bash
ProgDir=/home/armita/git_repos/emr_repos/scripts/alternaria/pathogen/circos
circos -conf $ProgDir/Ag_vs_At_genome_alignment/At_vs_Ag_circos.conf -outputdir $OutDir
mv $OutDir/circos.png $OutDir/At_vs_Ag_circos.png
mv $OutDir/circos.svg $OutDir/At_vs_Ag_circos.svg
```
<!--
```bash
OutDir=analysis/circos/At_vs_Ag_genome_alignment_circos
cat $OutDir/At_Ag_genome_edited2.txt | grep -v '4287' > $OutDir/At_Ag_genome_final.txt
mkdir -p $OutDir/by_FoC_chr
for Num in $(seq 1 22); do
  Chr="contig_"$Num"_pilon"
  echo "$Chr"
  OrthologyTxt=analysis/orthology/orthomcl/FoC_vs_Fo_vs_FoL_publication_ncbi/FoC_vs_Fo_vs_FoL_publication_ncbi_orthogroups.txt
  ProgDir=~/git_repos/emr_repossh s/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
  $ProgDir/orthology2ribbons_internal.py \
  --chr1 $Chr \
  --orthology $OrthologyTxt \
  --name1 Fus2 \
  --gff1 gene_pred/final_genes/F.oxysporum_fsp_cepae/Fus2_canu_new/final/final_genes_appended.gff3 \
  | sort | uniq \
  > $OutDir/At_vs_Ag_links_edited.txt
  ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
  circos -conf $ProgDir/Fus2/Fus2_FoL/At_vs_Ag_circos.conf -outputdir $OutDir
  mv $OutDir/circos.png $OutDir/by_FoC_chr/At_vs_Ag_LS_"$Chr"_circos.png
  mv $OutDir/circos.svg $OutDir/by_FoC_chr/At_vs_Ag_LS_"$Chr"_circos.svg
done
``` -->
<!--
The frequency of gene duplications within and between chromosomes was investigated:

```bash
OutDir=analysis/circos/At_vs_Ag_genome_alignment_circos
for Num in $(seq 1 22); do
Chr="contig_"$Num"_pilon"
ChrList="$ChrList $Chr"
done
echo "$ChrList"
OrthologyTxt=analysis/orthology/orthomcl/FoC_vs_Fo_vs_FoL_publication_ncbi/FoC_vs_Fo_vs_FoL_publication_ncbi_orthogroups.txt
ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
$ProgDir/orthology2ribbons_internal.py \
--chr1 $ChrList \
--orthology $OrthologyTxt \
--name1 Fus2 \
--gff1 gene_pred/final_genes/F.oxysporum_fsp_cepae/Fus2_canu_new/final/final_genes_appended.gff3 \
| sort | uniq > $OutDir/Fus2_all_Fus2_links.txt
cat $OutDir/Fus2_all_Fus2_links.txt | cut -f1,4 | sort | uniq -c > $OutDir/Fus2_all_Fus2_links_occurence.txt

``` -->
