assembly/Alt_database/Alternaria_alternata_ATCC11680/genome.ctg.fa



Quast and BUSCO

```bash
for Assembly in $(ls assembly/Alt_database/*/genome.ctg.fa); do
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
Alternaria_alternata	ATCC11680	1300	1298	5	10	1315
Alternaria_alternata	ATCC66891	1299	1297	7	9	1315
Alternaria_alternata	BMP0270	1302	1300	4	9	1315
Alternaria_arborescens	BMP0308	1013	1011	203	99	1315
Alternaria_brassicicola	ATCC96836	1132	1125	111	72	1315
Alternaria_capsici	BMP0180	1290	1289	9	16	1315
Alternaria_carthami	BMP1963	1261	1259	25	29	1315
Alternaria_citriarbusti	BMP2343	1264	1262	25	26	1315
Alternaria_crassa	BMP0172	1286	1284	13	16	1315
Alternaria_dauci	BMP0167	1011	1007	133	171	1315
Alternaria_destruens	BMP0317	439	438	480	396	1315
Alternaria_fragariae	BMP3062	1268	1266	24	23	1315
Alternaria_gaisen	BMP2338	1052	1051	157	106	1315
Alternaria_limoniasperae	BMP2335	1261	1259	26	28	1315
Alternaria_longipes	BMP0313	1299	1297	6	10	1315
Alternaria_macrospora	BMP1949	1212	1211	57	46	1315
Alternaria_mali	BMP3063	1153	1152	87	75	1315
Alternaria_mali	BMP3064	1233	1231	41	41	1315
Alternaria_porri	BMP0178	859	858	244	212	1315
Alternaria_solani	BMP0185	1178	1177	94	43	1315
Alternaria_tagetica	BMP0179	1294	1291	9	12	1315
Alternaria_tangelonis	BMP2327	1221	1219	48	46	1315
Alternaria_tenuissima	BMP0304	1302	1300	4	9	1315
Alternaria_tomatophila	BMP2032	1217	1215	68	30	1315
Alternaria_turkisafria	BMP3436	1211	1209	50	54	1315
```
