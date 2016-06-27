#!/bin/bash
mkdir -p ../cuffquant/abundances
for i in *.bam
do
  file=$(basename "$i")
  filename="${file%.*}"
  qsub ./cuff_quant.sh $file $filename
done
