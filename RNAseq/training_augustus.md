These commands were used to train augustus using RNAseq data.


# 1) QC

Perform qc of RNAseq timecourse data
```bash
  cd /home/groups/harrisonlab/project_files/alternaria/
  for FilePath in $(ls -d raw_rna/paired/A.*/*); do
    echo $FilePath
    FileF=$(ls $FilePath/F/*.gz)
    FileR=$(ls $FilePath/R/*.gz)
    IlluminaAdapters=/home/armita/git_repos/emr_repos/tools/seq_tools/ncbi_adapters.fa
    qsub /home/armita/git_repos/emr_repos/tools/seq_tools/rna_qc/rna_qc_fastq-mcf.sh $FileF $FileR $IlluminaAdapters RNA
  done
```

# 2) Align reads vs. host genomes
Alignments of RNAseq reads were made against the Fus2 Genome using tophat:

## 2.1) Alignment

```bash
  for FilePath in $(ls -d qc_rna/paired/A.*/*); do
    Treatment=$(echo $FilePath | rev | cut -f2 -d '/' | rev)
    Genome=repeat_masked/A.alternata_ssp._gaisen/650/A.alternata_ssp._gaisen_650_67_repmask/650_contigs_unmasked.fa
    OutDir=gene_pred/training_augustus/A.alternata_ssp._gaisen/650/$Treatment
    FileF=$(ls $FilePath/F/*_trim.fq.gz)
    FileR=$(ls $FilePath/R/*_trim.fq.gz)
    qsub /home/armita/git_repos/emr_repos/tools/seq_tools/RNAseq/tophat_alignment.sh $Genome $FileF $FileR $
  done
  for FilePath in $(ls -d qc_rna/paired/A.*/*); do
    Treatment=$(echo $FilePath | rev | cut -f2 -d '/' | rev)
    Genome=repeat_masked/A.alternata_ssp._tenuissima/1166/A.alternata_ssp._tenuissima_1166_43_repmask/1166_contigs_unmasked.fa
    OutDir=gene_pred/training_augustus/A.alternata_ssp._tenuissima/1166/$Treatment
    FileF=$(ls $FilePath/F/*_trim.fq.gz)
    FileR=$(ls $FilePath/R/*_trim.fq.gz)
    qsub /home/armita/git_repos/emr_repos/tools/seq_tools/RNAseq/tophat_alignment.sh $Genome $FileF $FileR $OutDir
  done
```

<!-- ## 2.2) Summarising alignments
Results were summarised using the following commands:

```bash
  for File in $(ls bowtie_alignment.sh.e*); do
    echo $File;
    Logfile=$(echo $File | sed "s/.sh.e/.sh.o/g");
    head -n 2 $Logfile;
    cat $File;
    echo "";  
  done
```

The text output of these alignments is stored in summised_output.txt -->

# 3) Move to timecourse directory

Data was copied from the alignment directory to a working directory for this
timecourse experiment.

```bash
  WorkDir=gene_pred/training_augustus
  mkdir -p $WorkDir
  cp -r alignment/A.* $WorkDir/.
```

# 4) Assemble transcripts

Cufflinks was used to assemble transcripts from reads aligned to the genome.

```bash
  for RnaSample in $(ls -d gene_pred/training_augustus/*/* | grep -v 'merged' | grep -v 'quantified'); do
  echo $RnaSample;
  Alignment=$(ls gene_pred/training_augustus/*/*/$RnaSample/accepted_hits.bam)
  echo $Alignment
  cufflinks -o gene_pred/training_augustus/*/*/$RnaSample/cufflinks -p 16 --max-intron-length 4000 $Alignment
  done
  # for Media in $(ls timecourse/v2_genes/F.oxysporum_fsp_cepae | grep -v 'merged' | grep -v 'quantified'); do
  # echo $Media;
  # Alignment=$(ls "timecourse/v2_genes/F.oxysporum_fsp_cepae/$Media/accepted_hits.bam")
  # cufflinks -o timecourse/v2_genes/Fus2/$Media/cufflinks -p 16 --max-intron-length 4000 $Alignment
  # done
```

# 5) Merge assembled transcripts

```bash
  ls timecourse/v2_genes/Fus2/*/cufflinks/transcripts.gtf | sort -g -k4 -t '/' > transcript_list.txt
  cuffmerge -o timecourse/v2_genes/Fus2/merged --num-threads 16 transcript_list.txt
  rm transcript_list.txt
```

# 6) quantify expression

```bash
for TimePoint in $(ls timecourse/v2_genes/Fus2 | grep -v 'merged' | grep -v 'quantified' | sort -g -k4 -t '/'); do
  echo $TimePoint;
  Alignment=$(ls timecourse/v2_genes/Fus2/$TimePoint/*_alignment.bam)
  echo $Alignment
  cuffquant timecourse/v2_genes/Fus2/merged/merged.gtf -o timecourse/v2_genes/Fus2/quantified/$TimePoint -p 16 $Alignment
done
```

```bash
  Alignments=$(ls timecourse/v2_genes/Fus2/quantified/*/abundances.cxb | sort -n -k5 -t '/')
  Labels=$(ls timecourse/v2_genes/Fus2/quantified/* | sort -n -k5 -t '/')
  cuffdiff --time-series -o timecourse/v2_genes/Fus2/quantified -p 16 timecourse/v2_genes/Fus2/merged/merged.gtf $Alignments
```
