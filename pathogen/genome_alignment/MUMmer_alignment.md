
# 5 Promer alignment of Assemblies

## 5.1 against 1166 genome

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

## 5.1 against 650 genome

```bash
Reference=$(ls repeat_masked/*/*/filtered_contigs/650_contigs_unmasked.fa)
for Query in $(ls repeat_masked/*/*/filtered_contigs/1166_contigs_unmasked.fa); do
Strain=$(echo $Query | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Query | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
Prefix="$Strain"_vs_650
OutDir=analysis/genome_alignment/mummer/$Organism/$Strain/$Prefix
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/genome_alignment/MUMmer
qsub $ProgDir/sub_nucmer.sh $Reference $Query $Prefix $OutDir
done
```

## 5.1 against the reference A. solani genome

```bash
Reference=$(ls ../../../../home/groups/harrisonlab/project_files/alternaria/assembly/misc_publications/Alternaria_solani_altNL03003/genome.ctg.fa)
for Query in $(ls repeat_masked/*/*/filtered_contigs/*_contigs_unmasked.fa | grep '650'); do
Strain=$(echo $Query | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Query | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
# Prefix="$Strain"_vs_A.solani
Prefix="$Strain"_vs_A.solani_2
OutDir=analysis/genome_alignment/mummer/$Organism/$Strain/$Prefix
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/genome_alignment/MUMmer
qsub $ProgDir/sub_nucmer.sh $Reference $Query $Prefix $OutDir
# $ProgDir/sub_nucmer.sh $Reference $Query $Prefix $OutDir
done
```
