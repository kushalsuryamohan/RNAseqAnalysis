# Pipeline for analyzing standard Illumina RNA-seq data

The pipeline follows the steps outlined below to analyze standard Illumina RNA-seq data: 
1) Pre-Alignment QC
2) Alignment
3) Provision for visualizing alignments using IGV
4) Alignment QC
5) Expression and Differential Expression Analysis
6) Reference-Free Abundance Estimation using Kallisto
7) Isoform Discovery and Alternative Expression
  a) Reference Guided Transcript Assembly
  b) Differential Splicing
  c) Splicing Visualization


##1) Pre-Alignment QC
Run FASTQC on RNA-seq reads to assess quality of data
If data are in folder $RNA_seq/data
```
cd $RNA_seq/data
fastqc *.fastq.gz
```
Alternatively, use http://www.bioinformatics.babraham.ac.uk/projects/fastqc/ to generate QC reports for each set of read data

## 2) Alignment

```
#Create directory for alignment data
mkdir alignments
cd alignments
for f in `cat ${files}`; do STAR --genomeDir ${GENOMEDIR} \
--readFilesIn fastq/$f\_R1.fastq fastq/$f\_R2.fastq \
--runThreadN 24 --outSAMtype BAM SortedByCoordinate --outFileNamePrefix aligned/$f.; done

```


## 3) 


