#!/bin/sh

CUFFLINKS_BINARY=/path/to/cufflinks/binary/cufflinks
P=4 #use 4 threads/cores

for SAMPLE_ID in {1..N}
do
cat > /ProjectName/assemblies/scripts/cufflinks_${SAMPLE_ID}.sh <<EOF 
#!/bin/sh

$CUFFLINKS_BINARY -q -p $P -o /ProjectName/assemblies/sample${SAMPLE_ID} /ProjectName/alignments/sample${SAMPLE_ID}.bam

mv /ProjectName/assemblies/sample${SAMPLE_ID}/transcripts.gtf /ProjectName/assemblies/sample${SAMPLE_ID}_transcripts.gtf
EOF
qsub -l mf=20G,h_vmem=5G -m e -M myemail@email.com -pe local $P cufflinks_${SAMPLE_ID}.sh
done
