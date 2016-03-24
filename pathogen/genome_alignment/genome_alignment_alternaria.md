# Genome alignment

Genome alignment was performed using the pipeline recomended by Eva Stukenbrock
at the Max Planck Institute for Evolutionary Biology, Kiel.

## 1) Generation of a phylogenetic tree for isolates

A phylogenetic tree was provided to TBA

```bash
TreeString='"((648 1082 1164 24350 (635 743 1166 1177)) (675 97.0013 97.0016) (650))"'
```

## 2) Generate multiple alignment between genomes

The TBA package http://www.bx.psu.edu/miller_lab/ was used to generate multiple
alignments.

### 2.1) generating a series of pair-wise alignments

generating a series of pair-wise alignments to “seed” the multiple alignment process

```bash
  ProjDir=/home/groups/harrisonlab/project_files/alternaria
  WorkDir=analysis/genome_alignment/tmp4
  mkdir -p $WorkDir
  cd $WorkDir

  for Assembly in $(ls $ProjDir/repeat_masked/*/*/filtered_contigs_repmask/*_contigs_softmasked.fa); do
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
    Strain=$(echo $Assembly| rev | cut -f3 -d '/' | rev | sed 's/\./-/g')
    # The commands to define strain also replace any '.' with a '-' character.
    cat $Assembly | sed "s/>/>Alt_$Strain:/g" > Alt_$Strain.fa
    get_standard_headers Alt_$Strain.fa > Alt_"$Strain"_headers.txt
    ProgDir=~/git_repos/emr_repos/tools/pathogen/lineage_specific_regions
    $ProgDir/parse_tba_headers.py --inp_fasta Alt_$Strain.fa --new_headers Alt_"$Strain"_headers.txt > Alt_$Strain
  done

  # Reference=$ProjDir/assembly/external_groups/A.brassicicola/EGS42–002/Alternaria_brassicicola_masked_assembly.fasta
  # cat $Reference | sed "s/>/>Alt_$Strain:/g" > Alt_$Strain.fa

  all_bz - \
    "((Alt_648 Alt_1082 Alt_1164 Alt_24350 (Alt_635 Alt_743 Alt_1166 Alt_1177)) (Alt_675 Alt_97-0013 Alt_97-0016) (Alt_650))" >& all_bz.log

  cat all_bz.log | grep 'blastzWrapper' > commands_part1.log
  while read Commands; do
    Nodes="blacklace02.blacklace|blacklace02.blacklace|blacklace03.blacklace|blacklace04.blacklace|blacklace05.blacklace|blacklace06.blacklace|blacklace07.blacklace|blacklace08.blacklace"
    qsub -S /bin/bash -b y -pe smp 1 -l virtual_free=1G -N blastzWrapper -l h="$Nodes" -cwd "$Commands"
  done < commands_part1.log

  cat all_bz.log | grep 'single_cov2' > commands_part2.log
  while read Commands; do
    qsub -S /bin/bash -b y -pe smp 1 -l virtual_free=1G -N single_cov2 -l h="$Nodes" -cwd "$Commands"
  done < commands_part2.log

```

### 2.2) generating the multiple alignment

```bash
# # tba - "(((648 1082 1164 24350 (635 743 1166 1177)) (675 97.0013 97.0016) 650))" *.*.maf tba.maf >& tba.log
# for Strain in 648 1082 1164 24350 635 743 1166 1177 675 97.0013 97.0016 650; do
# cat $Strain | sed 's/>/>A/g' > A$Strain
# done
# # maf_list=$(ls *.maf | grep -v -e '97.00' -e 'tba' -e '1082' -e '1164' -e '743' -e '24350' -e '1166' -e '1177' -e '675' | tr -d '\n' | sed 's/maf/maf /g')
# maf_list=$(ls *.*.maf | tr -d '\n' | sed 's/maf/maf /g')
# for Strain in 648 1082 1164 24350 635 743 1166 1177 675 97.0013 97.0016 650; do
# for File in $(ls *$Strain*); do
# echo $File;
# sed -i -E 's/^s /s A/g' $File;
# NewName=$(echo $File | sed -r "s/$Strain/A$Strain/g")
# mv $File $NewName
# done
# done
# maf_list=$(ls A*.maf)
tba E=A1177 "((Alt_648 Alt_1082 Alt_1164 Alt_24350 (Alt_635 Alt_743 Alt_1166 Alt_1177)) (Alt_675 Alt_97-0013 Alt_97-0016) (Alt_650))" A*.maf tba.maf >& tba.log
```

### 2.3) “projecting” the alignment onto a reference sequence

## 3) Process alignment file

The Alignment file was processed using Maffilter.
http://biopp.univ-montp2.fr/forge/maffilter
This also allowed calling of SNPs, extraction of CDS and preparation for
phylogenetic analysis.
Importantly, it allows identification of contigs that do not align to other
strains, facilitating identification of CDCs.

```bash
  mkdir ../maffilter
  cd ../maffilter
  # maffilter input.file=mydata.maf.gz input.file.compression=gzip input.format=Maf output.log=mydata.maffilter.log
  maffilter input.file=../tmp2/tba.maf input.file.compression=none input.format=Maf output.log=mydata.maffilter.log maf.filter=" \
  SequenceStatistics(                     \
      statistics=( BlockLength,                    \
          AlnScore,                       \
          BlockCounts,                       \
          SiteStatistics(                 \
                  species=A648 A1082 A1164 A24350 A635 A743 A1166 A1177 A675 A97.0013 A97.0016 A650      \
          ),                   \
        ),                   \
      ref_species=species1,               \
      file=data.statistics.csv       \
      )"

# Output(                                 \
#     file=data.filtered.maf,          \
#     compression=none,                   \
#     mask=yes) \
# "
# Subset(
#     species=(A635, A743, A1166, A1177), \
#     strict=yes, \
#     remove_duplicates=yes), \
# "Output(                                 \
#         file=data.filtered.maf,          \
#         compression=none,                   \
#         mask=yes)"
```







## Progressive Mauve was used to align genomes:

```bash
ProjDir=/home/groups/harrisonlab/project_files/alternaria
WorkDir=$ProjDir/analysis/genome_alignment/mauve
mkdir -p $WorkDir
# Reference=$ProjDir/assembly/external_groups/A.brassicicola/EGS42–002/Alternaria_brassicicola_masked_assembly.fasta
# cp $Reference $WorkDir/A.brasicicola.fa
Reference=$(ls $ProjDir/repeat_masked/*/*/filtered_contigs_repmask/*_contigs_softmasked.fa)
# cd $WorkDir

# Use move_contigs to order genomes based on reference for each phylogroup
for Assembly in $(ls $ProjDir/repeat_masked/*/*/filtered_contigs_repmask/*_contigs_softmasked.fa | grep -v '1177'); do
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
Strain=$(echo $Assembly| rev | cut -f3 -d '/' | rev)
echo "$Organism - $Strain"
MauveDir=~/prog/mauve/mauve_snapshot_2015-02-13
OutDir=$WorkDir/"$Strain"_contigs
ProgDir=~/git_repos/emr_repos/tools/seq_tools/genome_alignment/mauve
rm -r $OutDir
qsub $ProgDir/mauve_order_contigs.sh $MauveDir $Reference  $Assembly $OutDir
# qsub $ProgDir/mauve_order_contigs.sh $MauveDir $WorkDir/A.brasicicola.fa  $Assembly $OutDir
done


# Generate alignment of genome sequences using progressive mauve - remember to use linux-64
./progressiveMauve --output=/home/hulinm/local/src/mauve_2.3.1/all.xmfa  /home/hulinm/pseudomonas_data/pseudomonas/assembly/ordered_contigs/fastafiles/5244.fasta
```
