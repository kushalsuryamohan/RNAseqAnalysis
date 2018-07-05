
#!/bin/sh
set -e 
set -u
set -o pipefail

module load apps/gcc/4.9.0
module load apps/python/2.7.13
module load apps/java/8u112
module load apps/parallel/20180322 
module load apps/biopython/1.56
module load apps/bioperl/1.6.924

function usage(){
printf "\nDescription: This script is to assemble the transcriptome and produce alignments of transcripts\
The script can be split into 2 major steps: 1: Read Alignment, and 2: Reference-guided Transcript assembly 

$(basename "$0") [-h] [-OPTIONS] -p <Project name> -r <reference genome> \
-i <bowtie index> -I <path to Illumina reads> 
where:
   -h  show this help text
   
 \n\n"  
}

######### Parse input #########
while getopts 'h:d:s:P:N:I:B:p:n:S:' option; do
  case "$option" in
    h) usage
       exit 1
       ;;
    d) WORKDIR=${OPTARG}
       ;;
    s) GenomeSZ=${OPTARG}
       ;;
    *) usage
       exit 1
       ;;
  esac
done
shift $((OPTIND-1))

if [ -z "${WORKDIR}" ]; 
then
    echo -e "Error: no working directory provided"
    exit 1
fi

if [ -z "${Prefix}" ];
then
    echo -e "Error: no output prefix provided"
    exit 1
fi


for SAMPLE_ID in {1..N}
do
cat > /ProjectName/alignments/scripts/tophat_${SAMPLE_ID}.sh <<EOF
$TOPHAT_BINARY -G $GENE_REFERENCE -p $P -o /ProjectName/alignments/sample${SAMPLE_ID} $BOWTIE_INDEX /ProjectName/data/sample${SAMPLE_ID}_R1.fastq /ProjectName/data/sample${SAMPLE_ID}_R2.fastq
mv /ProjectName/alignments/sample${SAMPLE_ID}/accepted_hits.bam /ProjectName/alignments/sample${SAMPLE_ID}.bam
EOF

bsub /ProjectName/alignments/scripts/tophat_${SAMPLE_ID}.sh

done
