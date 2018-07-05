
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


for SAMPLE_ID in {1..N}
do
cat > /ProjectName/alignments/scripts/tophat_${SAMPLE_ID}.sh <<EOF
$TOPHAT_BINARY -G $GENE_REFERENCE -p $P -o /ProjectName/alignments/sample${SAMPLE_ID} $BOWTIE_INDEX /ProjectName/data/sample${SAMPLE_ID}_1.fastq /ProjectName/data/sample${SAMPLE_ID}_2.fastq
mv /ProjectName/alignments/sample${SAMPLE_ID}/accepted_hits.bam /ProjectName/alignments/sample${SAMPLE_ID}.bam
EOF

bsub /ProjectName/alignments/scripts/tophat_${SAMPLE_ID}.sh

done
