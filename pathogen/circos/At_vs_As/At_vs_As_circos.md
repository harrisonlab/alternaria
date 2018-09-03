
```bash
  OutDir=analysis/circos/At_vs_As_circos
  mkdir -p $OutDir

  At_genome=$(ls repeat_masked/*/*/filtered_contigs/*_contigs_softmasked_repeatmasker_TPSI_appended.fa | grep '1166')
  ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
  $ProgDir/fasta2circos.py --genome $At_genome --contig_prefix "At_" > $OutDir/At_genome.txt

  As_genome=$(ls /home/groups/harrisonlab/project_files/alternaria/assembly/misc_publications/Alternaria_solani_altNL03003/genome.ctg.fa)
  ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
  $ProgDir/fasta2circos.py --genome $As_genome --contig_prefix "As_" > $OutDir/As_genome.txt

  cat $OutDir/At_genome.txt > $OutDir/At_As_genome.txt
  tac $OutDir/As_genome.txt >> $OutDir/At_As_genome.txt

```

Telomere locations on contigs:

```bash
cat analysis/telomere/A.alternata_ssp_tenuissima/1166/telomere_hits_circos.txt | sed 's/contig/At_contig/g' | sort -k3 -n -t'_' > $OutDir/At_vs_As_telomere_hits.txt
cat /home/groups/harrisonlab/project_files/alternaria/analysis/telomere/A.solani/altNL03003/telomere_hits_circos.txt  | sed 's/CP/As_CP/g' | sort -k3 -n -t'_' >> $OutDir/At_vs_As_telomere_hits.txt

```

```bash
OutDir=analysis/circos/At_vs_As_circos
Coords=$(ls analysis/genome_alignment/mummer/A.alternata_ssp_tenuissima/1166/1166_vs_A.solani/1166_vs_A.solani_coords.tsv)
ProgDir=/home/armita/git_repos/emr_repos/scripts/alternaria/pathogen/genome_alignment
$ProgDir/nucmer_coords2circos.py --inp_coords $Coords --queery_id At --ref_id As > $OutDir/At_vs_As_links.txt
cat $OutDir/At_vs_As_links.txt > $OutDir/At_vs_As_links_edited.txt
```


A file showing contig orientations was made:
```bash
  cat $OutDir/At_As_genome.txt | cut -f3 -d ' ' | sed "s/$/; /g" | tr -d '\n' > $OutDir/At_contig_order.txt
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/circos
  $ProgDir/find_contig_orientation.py --links_file $OutDir/At_vs_As_links_edited.txt > $OutDir/At_vs_As_contig_orientation.txt
```


Contig order was selected by taking the first line of that file and then also
taking the reversed order of contigs using the command:

```bash
cat $OutDir/At_vs_As_contig_orientation.txt | grep -A1 'Order of all seen contigs' | tail -n1 | sed "s/, /\n/g" > tmp.txt
cat $OutDir/At_As_genome.txt | grep 'As' | grep -w -v -f tmp.txt | cut -f3 -d ' '| tr -d '\n' | sed 's/As/, As/g'
# cat $OutDir/At_As_genome.txt | grep 'Ag' | cut -f3 -d ' ' | tr -d '\n' | sed 's/Ag/, Ag/g' >> tmp.txt

echo "Order of unseen At contigs and remaining As contigs"
cat $OutDir/At_As_genome.txt | grep -w -v -f tmp.txt | cut -f3 -d ' '| tr -d '\n' | sed 's/At/, At/g' | sed 's/As/, As/g'
```

```bash
ProgDir=/home/armita/git_repos/emr_repos/scripts/alternaria/pathogen/circos
circos -conf $ProgDir/At_vs_As/At_vs_As_circos.conf -outputdir $OutDir
mv $OutDir/circos.png $OutDir/At_vs_As_circos.png
mv $OutDir/circos.svg $OutDir/At_vs_As_circos.svg
ls $PWD/$OutDir/At_vs_As_circos.png
```
