```bash
  OutDir=analysis/circos/Ag_vs_As_circos
  mkdir -p $OutDir

  Ag_genome=$(ls repeat_masked/*/*/filtered_contigs/*_contigs_softmasked_repeatmasker_TPSI_appended.fa | grep '650')
  ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
  $ProgDir/fasta2circos.py --genome $Ag_genome --contig_prefix "Ag_" > $OutDir/Ag_genome.txt

  As_genome=$(ls /home/groups/harrisonlab/project_files/alternaria/assembly/misc_publications/Alternaria_solani_altNL03003/genome.ctg.fa)
  ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
  $ProgDir/fasta2circos.py --genome $As_genome --contig_prefix "As_" > $OutDir/As_genome.txt

  cat $OutDir/Ag_genome.txt > $OutDir/Ag_As_genome.txt
  tac $OutDir/As_genome.txt >> $OutDir/Ag_As_genome.txt
```

Telomere locations on contigs:

```bash
cat analysis/telomere/A.gaisen/650/telomere_hits_circos.txt | sed 's/contig/Ag_contig/g' | sort -k3 -n -t'_' > $OutDir/Ag_vs_As_telomere_hits.txt
cat /home/groups/harrisonlab/project_files/alternaria/analysis/telomere/A.solani/altNL03003/telomere_hits_circos.txt  | sed 's/CP/As_CP/g' | sort -k3 -n -t'_' >> $OutDir/Ag_vs_As_telomere_hits.txt
```

```bash
OutDir=analysis/circos/Ag_vs_As_circos
Coords=$(ls analysis/genome_alignment/mummer/A.gaisen/650/650_vs_A.solani/650_vs_A.solani_coords.tsv)
ProgDir=/home/armita/git_repos/emr_repos/scripts/alternaria/pathogen/genome_alignment
$ProgDir/nucmer_coords2circos.py --inp_coords $Coords --queery_id Ag --ref_id As > $OutDir/Ag_vs_As_links.txt
cat $OutDir/Ag_vs_As_links.txt > $OutDir/Ag_vs_As_links_edited.txt
```


A file showing contig orientations was made:
```bash
  cat $OutDir/Ag_As_genome.txt | cut -f3 -d ' ' | sed "s/$/; /g" | tr -d '\n' > $OutDir/Ag_contig_order.txt
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/circos
  $ProgDir/find_contig_orientation.py --links_file $OutDir/Ag_vs_As_links_edited.txt > $OutDir/Ag_vs_As_contig_orientation.txt
```


Contig order was selected by taking the first line of that file and then also
taking the reversed order of contigs using the command:

```bash
cat $OutDir/Ag_vs_As_contig_orientation.txt | grep -A1 'Order of all seen contigs' | tail -n1 | sed "s/, /\n/g" > tmp.txt
cat $OutDir/Ag_As_genome.txt | grep 'As' | grep -w -v -f tmp.txt | cut -f3 -d ' '| tr -d '\n' | sed 's/As/, As/g'
# cat $OutDir/Ag_As_genome.txt | grep 'Ag' | cut -f3 -d ' ' | tr -d '\n' | sed 's/Ag/, Ag/g' >> tmp.txt

echo "Order of unseen Ag contigs and remaining As contigs"
cat $OutDir/Ag_As_genome.txt | grep -w -v -f tmp.txt | cut -f3 -d ' '| tr -d '\n' | sed 's/Ag/, Ag/g' | sed 's/As/, As/g'

```

```bash
ProgDir=/home/armita/git_repos/emr_repos/scripts/alternaria/pathogen/circos
circos -conf $ProgDir/Ag_vs_As/Ag_vs_As_circos.conf -outputdir $OutDir
mv $OutDir/circos.png $OutDir/Ag_vs_As_circos.png
mv $OutDir/circos.svg $OutDir/Ag_vs_As_circos.svg
ls $PWD/$OutDir/Ag_vs_As_circos.png
```
<!--
```bash
OutDir=analysis/circos/Ag_vs_As_circos
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
  > $OutDir/Ag_vs_As_links_edited.txt
  ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
  circos -conf $ProgDir/Fus2/Fus2_FoL/Ag_vs_As_circos.conf -outputdir $OutDir
  mv $OutDir/circos.png $OutDir/by_FoC_chr/Ag_vs_As_LS_"$Chr"_circos.png
  mv $OutDir/circos.svg $OutDir/by_FoC_chr/Ag_vs_As_LS_"$Chr"_circos.svg
done
``` -->
<!--
The frequency of gene duplications within and between chromosomes was investigated:

```bash
OutDir=analysis/circos/Ag_vs_As_circos
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
