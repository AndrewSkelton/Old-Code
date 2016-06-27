#!/bin/bash
mkdir ACAN_upstream
for i in *.bam
do
  filename="${i%.*}"
  samtools view -h -o ./ACAN_upstream/$filename'.sam' $i 15:89338134-89347877
  samtools view -bS ./ACAN_upstream/$filename'.sam' > ./ACAN_upstream/$filename'_ACAN_Upstream.bam'
  samtools index ./ACAN_upstream/$filename'_ACAN_Upstream.bam'
done
