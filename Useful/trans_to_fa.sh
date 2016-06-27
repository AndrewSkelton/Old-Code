#!/bin/bash

grep $1 RNA_Seq_Comparison_NOF_OA_Cell-Lines/merged.gtf > ./tmp.gtf

gffread -w $1'.fa' -g /opt/databases/genomes/Ensembl_2/genome.fa ./tmp.gtf 

rm tmp.gtf
