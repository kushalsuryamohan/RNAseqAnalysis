#!/bin/sh
set -e
set -u
set -o pipefail

module load apps/cufflinks

function usage(){
printf "\nThis script is to assemble the transcriptome using a reference genome and includes the following steps: \
1. Build genome reference index using STAR
2. Alignment/Mapping of each sample's RNA-seq reads to the supplied reference genome using STAR;  \
2. Assemble transcripts for each sample using Cufflinks; \
3. Merge transcripts from individual sample to generate multi-tissue 'transcriptome' using Cuffmerge;
4. Visualize transcriptome
5. Identify differentially expressed transcripts using Cuffdiff

$(basename "$0") [-h] [-OPTIONS] -d <working dir> -g <path to genome file> -I <text file of paths to each sample's RNA-seq dataset> \
-n <number of CPU cores> -star <Path to STAR binary> 
where:
    -h  show this help text
    -d  PATH to the working directory for transcriptome assembly
    -I  path to text file containing list of paths to RNA-seq data
    -n  number of CPU cores
    -g  input genome reference fasta file
    -star Path to STAR binary
    -S  the steps to processed (default=1234): 1: Build Index; 2. Align reads; 
        3. Assemble transcripts; 4. Merge transcript models; 5. Visualization; 6. Differential expression analysis;
        \n\n"
}


######### Default parameters #########  
WORKDIR=''
GENOME=''

ThreadN=24
ProcStep=1234 


######### Parse input #########
while getopts 'h:d:I:n:g:star:S:' option; do
  case "$option" in
    h) usage
       exit
       ;;
    d) WORKDIR=${OPTARG}
       ;;
    g) GENOME=${OPTARG}
       ;;
    I) ILMN_PATH=${OPTARG}
       ;;
    star) STAR=${OPTARG}
       ;;
    n) ThreadN=${OPTARG}
       ;;
    S) ProcStep=${OPTARG}
       ;;
    *) usage
       exit 1
       ;;
  esac
done    
shift $((OPTIND-1))


if [ -z "${WORKDIR}" ] 
then
    echo -e "ERROR: no working directory provided"
    exit 1
fi

if [ -z "${GENOME}" ] 
then
    echo -e "ERROR: no genome sequence provided"
    exit 1
fi

if [ -z "${STAR}" ] 
then
    echo -e "ERROR: Path to STAR not specified"
    exit 1 
fi

if [ -z "${ILMN_PATH}" ]
then
  echo -e "ERROR: Illumina RNA-seq data filepath text file not provided. Please provide path to a text file containing paths to each sample's RNA-seq data"
  exit 1



cd ${WORKDIR}

## 1. Build reference genome index
# Usage: STAR  --runMode genomeGenerate --runThreadN <# cpus> --genomeDir <genome output directory> --genomeFastaFiles <input Genome FASTA file>
bsub -e build_index.err -o build_index.out -J BUILD_INDEX -R "rusage[mem=128]" "${STAR} --runMode genomeGenerate --runThreadN ${ThreadN} --genomeDir ./ --genomeFastaFiles ${GENOME}"

# paried-end data:
#Usage: STAR --genomeDir mm9-starIndex/  --runThreadN 24 --readFilesIn read1.fastq read2.fastq --outFileNamePrefix Experiment1Star --genomeLoad LoadAndKeep
for path in `cat ${ILMN_PATH}`
do
  array=(${path//// }) #split line on '/' to get filename
  array_len=$((${#array[@]} - 1)) #length of string array
  tmp=${array[${array_len}]}
  filename=(${tmp//./ })
  bsub -e build_index.err -o build_index.out -J BUILD_INDEX -R "rusage[mem=128]" "${STAR} --genomeDir ./  --runThreadN 48 --readFilesIn ${path}/*_R1.fastq.gz ${path}/*_R2.fastq.gz --outFileNamePrefix ${filename[0]}"

mv /ProjectName/assemblies/sample${SAMPLE_ID}/transcripts.gtf /ProjectName/assemblies/sample${SAMPLE_ID}_transcripts.gtf
EOF
bsub cufflinks_${SAMPLE_ID}.sh
done
