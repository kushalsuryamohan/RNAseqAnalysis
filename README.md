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
hisat2 [options]* -x <ht2-idx> {-1 <m1> -2 <m2> | -U <r> | --sra-acc <SRA accession number>} [-S <sam>]

where,
'-p 8' tells HISAT2 to use eight CPUs for bowtie alignments.
'--rna-strandness RF' specifies strandness of RNAseq library. We will specify RF since the TruSeq strand-specific library was used to make these libraries. See here for options.
'--rg-id $ID' specifies a read group ID that is a unique identifier.
'--rg SM:$SAMPLE_NAME' specifies a read group sample name. This together with rg-id will allow you to determine which reads came from which sample in the merged bam later on.
'--rg LB:$LIBRARY_NAME' specifies a read group library name. This together with rg-id will allow you to determine which reads came from which library in the merged bam later on.
'--rg PL:ILLUMINA' specifies a read group sequencing platform.
'--rg PU:$PLATFORM_UNIT' specifies a read group sequencing platform unit. Typically this consists of FLOWCELL-BARCODE.LANE
'--dta' Reports alignments tailored for transcript assemblers.
'-x /path/to/hisat2/index' The HISAT2 index filename prefix (minus the trailing .X.ht2) built earlier including splice sites and exons.
'-1 /path/to/read1.fastq.gz' The read 1 FASTQ file, optionally gzip(.gz) or bzip2(.bz2) compressed.
'-2 /path/to/read2.fastq.gz' The read 2 FASTQ file, optionally gzip(.gz) or bzip2(.bz2) compressed.
'-S /path/to/output.sam' The output SAM format text file of alignments.

```
#Create directory for alignment data
mkdir alignments

```


## 3) 


