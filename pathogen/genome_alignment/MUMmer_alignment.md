
# 5 Promer alignment of Assemblies

## 5.1 against Fus2 genome

MUMmer was run to align assemblies against the reference genome.

```bash
Reference=$(ls repeat_masked/*/*/filtered_contigs/1166_contigs_unmasked.fa)
for Query in $(ls repeat_masked/*/*/filtered_contigs/650_contigs_unmasked.fa); do
Strain=$(echo $Query | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Query | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
Prefix="$Strain"_vs_1166
OutDir=analysis/genome_alignment/mummer/$Organism/$Strain/$Prefix
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/genome_alignment/MUMmer
qsub $ProgDir/sub_nucmer.sh $Reference $Query $Prefix $OutDir
done
```

```bash
Reference=$(ls repeat_masked/*/*/filtered_contigs/650_contigs_hardmasked_repeatmasker_TPSI_appended.fa)
for Query in $(ls repeat_masked/*/*/filtered_contigs/1166_contigs_hardmasked_repeatmasker_TPSI_appended.fa); do
Strain=$(echo $Query | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Query | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
Prefix="$Strain"_vs_650
OutDir=analysis/genome_alignment/mummer/$Organism/$Strain/$Prefix
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/genome_alignment/MUMmer
qsub $ProgDir/sub_nucmer.sh $Reference $Query $Prefix $OutDir
done
```
