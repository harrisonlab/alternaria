
Contents of the alternaria genomes database were downloaded:
```bash
assembly/Alt_database/Alternaria_alternata_ATCC11680/genome.ctg.fa
```

Additional genomes were downloaded from other publications:

```bash
Species='Alternaria_brassicicola'
Strain='Abra3'
OutDir=assembly/misc_publications/${Species}_${Strain}
mkdir -p $OutDir
wget -P $OutDir ftp://ftp.ncbi.nlm.nih.gov/sra/wgs_aux/PH/FN/PHFN01/PHFN01.1.fsa_nt.gz
gunzip $OutDir/*.gz
cp -s $PWD/$OutDir/*.fsa_nt $OutDir/genome.ctg.fa

Species='Alternaria_brassicicola'
Strain='ATCC96836'
OutDir=assembly/misc_publications/${Species}_${Strain}
mkdir -p $OutDir
wget -P $OutDir ftp://ftp.ncbi.nlm.nih.gov/sra/wgs_aux/AC/IW/ACIW01/ACIW01.1.fsa_nt.gz
gunzip $OutDir/*.gz
cp -s $PWD/$OutDir/*.fsa_nt $OutDir/genome.ctg.fa

Species='Alternaria_alternata'
Strain='B2a'
OutDir=assembly/misc_publications/${Species}_${Strain}
mkdir -p $OutDir
wget -P $OutDir ftp://ftp.ncbi.nlm.nih.gov/sra/wgs_aux/LS/HC/LSHC01/LSHC01.1.fsa_nt.gz
gunzip $OutDir/*.gz
cp -s $PWD/$OutDir/*.fsa_nt $OutDir/genome.ctg.fa

Species='Alternaria_alternata'
Strain='Z7'
OutDir=assembly/misc_publications/${Species}_${Strain}
mkdir -p $OutDir
wget -P $OutDir ftp://ftp.ncbi.nlm.nih.gov/sra/wgs_aux/LP/VP/LPVP01/LPVP01.1.fsa_nt.gz
gunzip $OutDir/*.gz
cp -s $PWD/$OutDir/*.fsa_nt $OutDir/genome.ctg.fa

Species='Alternaria_alternata'
Strain='ATCC34957'
OutDir=assembly/misc_publications/${Species}_${Strain}
mkdir -p $OutDir
wget -P $OutDir ftp://ftp.ncbi.nlm.nih.gov/sra/wgs_aux/LM/XP/LMXP01/LMXP01.1.fsa_nt.gz
gunzip $OutDir/*.gz
cp -s $PWD/$OutDir/*.fsa_nt $OutDir/genome.ctg.fa


Species='Alternaria_arborescens'
Strain='EGS39-128'
OutDir=assembly/misc_publications/${Species}_${Strain}
mkdir -p $OutDir
wget -P $OutDir ftp://ftp.ncbi.nlm.nih.gov/sra/wgs_aux/AI/IC/AIIC01/AIIC01.1.fsa_nt.gz
gunzip $OutDir/*.gz
cp -s $PWD/$OutDir/*.fsa_nt $OutDir/genome.ctg.fa

Species='Alternaria_solani'
Strain='altNL03003'
OutDir=assembly/misc_publications/${Species}_${Strain}
mkdir -p $OutDir
wget -P $OutDir ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/952/155/GCA_002952155.1_ASM295215v1/GCA_002952155.1_ASM295215v1_genomic.fna.gz
gunzip $OutDir/*.gz
cp -s $PWD/$OutDir/*.fna $OutDir/genome.ctg.fa

# Species='Ulocladium_consortalis'
# Strain='JCM1940'
# OutDir=assembly/misc_publications/${Species}_${Strain}
# mkdir -p $OutDir
# wget -P $OutDir ftp://ftp.ncbi.nlm.nih.gov/sra/wgs_aux/BC/GG/BCGG01/BCGG01.1.fsa_nt.gz
# gunzip $OutDir/*.gz
# cp -s $PWD/$OutDir/*.fsa_nt $OutDir/genome.ctg.fa
```



Quast and BUSCO

```bash
for Assembly in $(ls assembly/*/*/genome.ctg.fa | grep -e 'Alt_database' -e 'misc_publications' | grep 'misc_publications' | grep 'solani'); do
Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev | cut -f3 -d '_')
Organism=$(echo $Assembly | rev | cut -f2 -d '/' | rev | cut -f1,2 -d '_')
OutDir=$(dirname $Assembly)
echo "$Organism - $Strain"
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
qsub $ProgDir/sub_quast.sh $Assembly $OutDir
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/ascomycota_odb9)
# OutDir=$(dirname $Assembly)
OutDir=gene_pred/busco/$Organism/$Strain/assembly
qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
done
```


```bash
  for File in $(ls gene_pred/busco/*/*/assembly/*/short_summary_*.txt); do
  Strain=$(echo $File| rev | cut -d '/' -f4 | rev)
  Organism=$(echo $File | rev | cut -d '/' -f5 | rev)
  Complete=$(cat $File | grep "(C)" | cut -f2)
  Single=$(cat $File | grep "(S)" | cut -f2)
  Fragmented=$(cat $File | grep "(F)" | cut -f2)
  Missing=$(cat $File | grep "(M)" | cut -f2)
  Total=$(cat $File | grep "Total" | cut -f2)
  echo -e "$Organism\t$Strain\t$Complete\t$Single\t$Fragmented\t$Missing\t$Total"
  done
```

```
A.alternata_ssp._arborescens	675	1299	1296	6	10	1315
A.alternata_ssp._arborescens	97.0013	1299	1297	6	10	1315
A.alternata_ssp._arborescens	97.0016	1301	1298	3	11	1315
A.alternata_ssp._gaisen	650	1298	1296	4	13	1315
A.alternata_ssp._tenuissima	1082	1299	1296	4	12	1315
A.alternata_ssp._tenuissima	1164	1302	1299	3	10	1315
A.alternata_ssp_tenuissima	1166	1302	1301	3	10	1315
A.alternata_ssp._tenuissima	1166	1302	1299	5	8	1315
A.alternata_ssp._tenuissima	1177	1303	1299	3	9	1315
A.alternata_ssp._tenuissima	24350	1304	1302	2	9	1315
A.alternata_ssp._tenuissima	635	1303	1299	3	9	1315
A.alternata_ssp._tenuissima	648	1304	1302	3	8	1315
A.alternata_ssp._tenuissima	743	1303	1300	4	8	1315
A.gaisen	650	1260	1257	5	50	1315

Alternaria_alternata    B2a     1286    1284    6       23      1315
Alternaria_alternata    BMP0270 1302    1300    4       9       1315
Alternaria_alternata    Z7      1302    1298    4       9       1315
Alternaria_arborescens  BMP0308 1013    1011    203     99      1315
Alternaria_arborescens  EGS39-128       1301    1299    4       10      1315
Alternaria_brassicicola ATCC96836       1127    1125    106     82      1315
Alternaria_brassicicola Abra3   1299    1297    5       11      1315
Alternaria_capsici      BMP0180 1290    1289    9       16      1315
Alternaria_carthami     BMP1963 1261    1259    25      29      1315
Alternaria_citriarbusti BMP2343 1264    1262    25      26      1315
Alternaria_crassa       BMP0172 1286    1284    13      16      1315
Alternaria_dauci        BMP0167 1011    1007    133     171     1315
Alternaria_destruens    BMP0317 439     438     480     396     1315
Alternaria_fragariae    BMP3062 1268    1266    24      23      1315
Alternaria_gaisen       BMP2338 1052    1051    157     106     1315
Alternaria_limoniasperae        BMP2335 1261    1259    26      28      1315
Alternaria_longipes     BMP0313 1299    1297    6       10      1315
Alternaria_macrospora   BMP1949 1212    1211    57      46      1315
Alternaria_mali BMP3063 1153    1152    87      75      1315
Alternaria_mali BMP3064 1233    1231    41      41      1315
Alternaria_porri        BMP0178 859     858     244     212     1315
Alternaria_solani       BMP0185 1178    1177    94      43      1315
Alternaria_tagetica     BMP0179 1294    1291    9       12      1315
Alternaria_tangelonis   BMP2327 1221    1219    48      46      1315
Alternaria_tenuissima   BMP0304 1302    1300    4       9       1315
Alternaria_tomatophila  BMP2032 1217    1215    68      30      1315
Alternaria_turkisafria  BMP3436 1211    1209    50      54      1315
```

Summarise contig IDs by isolate
```bash
for File in $(ls assembly/Alt_database/*/genome.ctg.fa); do
  Species=$(echo $File | cut -f3 -d '/' | cut -f1,2 -d '_' | sed 's/lternaria_/. /g')
  Strain=$(echo $File | cut -f3 -d '/' | cut -f3 -d '_' )
  Abv=$(cat $File | grep '>' | head -n1 | tr -d '>' | sed "s/CTG.*//g")
  printf "$Species\t$Strain\t$Abv\t$Species $Strain\n"
done
```

```
A. alternata	ATCC11680	ATN	A. alternata ATCC11680
A. alternata	ATCC66891	AAT	A. alternata ATCC66891
A. alternata	BMP0270	AA2	A. alternata BMP0270
A. arborescens	BMP0308	AAB	A. arborescens BMP0308
A. brassicicola	ATCC96836	ABR	A. brassicicola ATCC96836
A. capsici	BMP0180	ACS	A. capsici BMP0180
A. carthami	BMP1963	ACM	A. carthami BMP1963
A. citriarbusti	BMP2343	ACT	A. citriarbusti BMP2343
A. crassa	BMP0172	ACR	A. crassa BMP0172
A. dauci	BMP0167	ADC	A. dauci BMP0167
A. destruens	BMP0317	ADT	A. destruens BMP0317
A. fragariae	BMP3062	AFG	A. fragariae BMP3062
A. gaisen	BMP2338	AGS	A. gaisen BMP2338
A. limoniasperae	BMP2335	ATK	A. limoniasperae BMP2335
A. longipes	BMP0313	ALG	A. longipes BMP0313
A. macrospora	BMP1949	AMR	A. macrospora BMP1949
A. mali	BMP3063	AML	A. mali BMP3063
A. mali	BMP3064	AM2	A. mali BMP3064
A. porri	BMP0178	APR	A. porri BMP0178
A. solani	BMP0185	ASL	A. solani BMP0185
A. tagetica	BMP0179	ATT	A. tagetica BMP0179
A. tangelonis	BMP2327	ALM	A. tangelonis BMP2327
A. tenuissima	BMP0304	AT2	A. tenuissima BMP0304
A. tomatophila	BMP2032	ATM	A. tomatophila BMP2032
A. turkisafria	BMP3436	ATG	A. turkisafria BMP3436
```


## Mating type genes:

```bash
  # mkdir -p analysis/blast/MAT_genes
  for Subject in $(ls assembly/*/*/genome.ctg.fa); do
    Strain=$(echo $Subject | rev | cut -f2 -d '/' | rev | cut -f3 -d '_')
    Organism=$(echo $Subject | rev | cut -f2 -d '/' | rev | cut -f1,2 -d '_')
    OutDir=analysis/blast_homology/$Organism/$Strain
    ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
    Query=analysis/blast_homology/MAT_genes/MAT_genes.fasta
    qsub $ProgDir/blast_pipe.sh $Query dna $Subject $OutDir
  done
```

```bash
qlogin -pe smp 12
cd /home/groups/harrisonlab/project_files/alternaria
QueryFasta=$(ls /data/scratch/armita/alternaria/analysis/blast/MAT_genes/MAT_genes.fasta)
dbType="nucl"
# for dbFasta in $(ls assembly/Alt_database/*/genome.ctg.fa); do
for dbFasta in $(ls assembly/misc_publications/*/genome.ctg.fa); do
Strain=$(echo $dbFasta | rev | cut -f2 -d '/' | rev | cut -f3 -d '_')
Organism=$(echo $dbFasta | rev | cut -f2 -d '/' | rev | cut -f1,2 -d '_')
echo "$Organism - $Strain"
Prefix="${Strain}_MAT"
Eval="1e-30"
OutDir=analysis/blast_homology/$Organism/$Strain
mkdir -p $OutDir
makeblastdb -in $dbFasta -input_type fasta -dbtype $dbType -title $Prefix.db -parse_seqids -out $OutDir/$Prefix.db
blastn -num_threads 12 -db $OutDir/$Prefix.db -query $QueryFasta -outfmt 6 -num_alignments 1 -out $OutDir/${Prefix}_hits.txt -evalue $Eval
done
```

```bash
for File in $(ls analysis/blast_homology/*/*/*_MAT_hits.txt); do
  Organism=$(echo $File | rev | cut -f3 -d '/' | rev)
  Strain=$(echo $File | rev | cut -f2 -d '/' | rev)
  for Hit in $(cat $File | cut -f1 | cut -f1 -d '_'); do
    printf "$Organism\t$Strain\t$Hit\n"
  done
done
```

```
A.alternata_ssp._arborescens	675	MAT1-2-1
A.alternata_ssp._arborescens	97.0013	MAT1-1-1
A.alternata_ssp._arborescens	97.0016	MAT1-2-1
A.alternata_ssp._gaisen	650	MAT1-1-1
A.alternata_ssp._tenuissima	1082	MAT1-2-1
A.alternata_ssp._tenuissima	1164	MAT1-2-1
A.alternata_ssp._tenuissima	1166	MAT1-2-1
A.alternata_ssp._tenuissima	1177	MAT1-1-1
A.alternata_ssp._tenuissima	24350	MAT1-1-1
A.alternata_ssp._tenuissima	635	MAT1-2-1
A.alternata_ssp._tenuissima	648	MAT1-1-1
A.alternata_ssp._tenuissima	743	MAT1-2-1
Alternaria_alternata	ATCC11680	MAT1-1-1
Alternaria_alternata	ATCC34957	MAT1-2-1
Alternaria_alternata	ATCC66891	MAT1-1-1
Alternaria_alternata	B2a	MAT1-1-1
Alternaria_alternata	BMP0270	MAT1-1-1
Alternaria_alternata	Z7	MAT1-1-1
Alternaria_arborescens	BMP0308	MAT1-2-1
Alternaria_arborescens	EGS39-128	MAT1-2-1
Alternaria_brassicicola	ATCC96836	MAT1-2-1
Alternaria_brassicicola	Abra3	MAT1-1-1
Alternaria_capsici	BMP0180	MAT1-2-1
Alternaria_carthami	BMP1963	MAT1-1-1
Alternaria_citriarbusti	BMP2343	MAT1-1-1
Alternaria_crassa	BMP0172	MAT1-1-1
Alternaria_dauci	BMP0167	MAT1-2-1
Alternaria_destruens	BMP0317	MAT1-1-1
Alternaria_destruens	BMP0317	MAT1-1-1
Alternaria_fragariae	BMP3062	MAT1-1-1
Alternaria_gaisen	BMP2338	MAT1-1-1
Alternaria_limoniasperae	BMP2335	MAT1-2-1
Alternaria_longipes	BMP0313	MAT1-1-1
Alternaria_macrospora	BMP1949	MAT1-1-1
Alternaria_mali	BMP3063	MAT1-1-1
Alternaria_mali	BMP3064	MAT1-1-1
Alternaria_porri	BMP0178	MAT1-1-1
Alternaria_solani	BMP0185	MAT1-1-1
Alternaria_tagetica	BMP0179	MAT1-1-1
Alternaria_tangelonis	BMP2327	MAT1-2-1
Alternaria_tenuissima	BMP0304	MAT1-2-1
Alternaria_tomatophila	BMP2032	MAT1-2-1
Alternaria_turkisafria	BMP3436	MAT1-1-1
```


## Presence and genes with homology to Alternaria toxins

```bash
for Fasta in $(ls assembly/*/*/genome.ctg.fa); do
Strain=$(echo $Fasta | rev | cut -f2 -d '/' | rev | cut -f3 -d '_')
Organism=$(echo $Fasta | rev | cut -f2 -d '/' | rev | cut -f1,2 -d '_')
cp $Fasta $(dirname $Fasta)"/${Organism}_${Strain}_contigs.fa"
done
for File in $(ls assembly/*/*/genes.gff3); do
Strain=$(echo $File | rev | cut -f2 -d '/' | rev | cut -f3 -d '_')
Organism=$(echo $File | rev | cut -f2 -d '/' | rev | cut -f1,2 -d '_')
cp $File $(dirname $File)"/${Organism}_${Strain}_contigs.gff3"
done
```

```bash
qlogin -pe smp 12
cd /home/groups/harrisonlab/project_files/alternaria
QueryFasta=$(ls analysis/blast_homology/CDC_genes/A.alternata_CDC_genes.fa)
dbType="nucl"
for dbFasta in $(ls assembly/*/*/genome.ctg.fa); do
Strain=$(echo $dbFasta | rev | cut -f2 -d '/' | rev | cut -f3 -d '_')
Organism=$(echo $dbFasta | rev | cut -f2 -d '/' | rev | cut -f1,2 -d '_')
echo "$Organism - $Strain"
Prefix="${Strain}_toxgenes"
Eval="1e-30"
OutDir=analysis/blast_homology/$Organism/$Strain
mkdir -p $OutDir
makeblastdb -in $dbFasta -input_type fasta -dbtype $dbType -title $Prefix.db -parse_seqids -out $OutDir/$Prefix.db
blastn -num_threads 4 -db $OutDir/$Prefix.db -query $QueryFasta -outfmt 6 -num_alignments 1 -out $OutDir/${Prefix}_hits.txt -evalue $Eval
done
```

```bash
for File in $(ls analysis/blast_homology/*/*/*_toxgenes_hits.txt); do
  Organism=$(echo $File | rev | cut -f3 -d '/' | rev)
  Strain=$(echo $File | rev | cut -f2 -d '/' | rev)
  for Hit in $(cat $File | cut -f1 | cut -f1 -d ' '); do
    printf "$Organism\t$Strain\t$Hit\n"
  done
done > analysis/blast_homology/CDC_genes/ref_genome_summarised_hits.tsv
```

```bash
for File in $(ls analysis/blast_homology/*/*/*_toxgenes_hits.txt); do
  Organism=$(echo $File | rev | cut -f3 -d '/' | rev)
  Strain=$(echo $File | rev | cut -f2 -d '/' | rev)
  cat $File | awk '{OFS="\t"; print $2, "toxin_homolog", "tblastx", $9, $10, $11, "+", ".", "ID=\"" $1 "_blast_hit1\";"}' | awk '$4 > $5 { temp = $4; $4 = $5; $5 = temp; $7 = "-" } 1' OFS='\t' | sed 's/gene_blast/blast/g' > ${File%.txt}.gff
done
```

```bash
for File in $(ls analysis/blast_homology/Alternaria_*/*/*_toxgenes_hits.txt); do echo $File;
cat $File | sort -nr -k3 | head -n5;
done
```

## Telomere identification

Telomeres were identified within the A. solani genome


```bash
  for Assembly in $(ls assembly/misc_publications/Alternaria_solani_altNL03003/genome.ctg.fa); do
  Organism="A.solani"
  Strain="altNL03003"
  echo "$Organism - $Strain"
  OutDir=analysis/telomere/$Organism/$Strain
  mkdir -p $OutDir
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/telomeres
  $ProgDir/annotate_telomeres.py --fasta $Assembly --out $OutDir/telomere_hits
  done
  # cat $OutDir/telomere_hits.txt | sort -nr -k7 | less
```
