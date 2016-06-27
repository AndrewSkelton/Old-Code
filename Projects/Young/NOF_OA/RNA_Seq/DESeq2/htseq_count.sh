#!/bin/bash

#-------------------------------------------------------------------------------------------#
#  Author      : Andrew J Skelton                                                           |
#  Language    : Bash                                                                       |
#  Study       : NOF and OA RNA Seq Study                                                   |
#  Data Owner  : Newcastle University - Prof. David Young                                   |
#  Description : Use htseq_count to count features at the gene level utilising GRCh38       |
#-------------------------------------------------------------------------------------------#



##'Set script variables
##'-----------------------------------------------------------------------------------------#
IN_DIR="/home/andrew/Young/Pipeline_Analysis/tophat/tophat_grch38/alignment/sorted"
REF_GTF="/data/genomes/GRCh38_Ensembl_Homo_Sapiens/hs_grch38.gtf"
COUNT_OUT="/home/andrew/Young/Pipeline_Analysis/htseq_count/grch38"
##'-----------------------------------------------------------------------------------------#



##'Loop through filenames and set a child process of htseq_count running
##'-----------------------------------------------------------------------------------------#
declare -a arr=("A" "B")
for i in "${arr[@]}"
do
  for j in {1..8}
  do
    OUT_PREFIX='Flowcell_'${i}'_'${j}
    htseq-count -f bam \
                -s no \
                -m intersection-nonempty \
                ${IN_DIR}/${OUT_PREFIX}'.sorted.bam' \
                ${REF_GTF} > ${COUNT_OUT}/${OUT_PREFIX}'.count' &
  done
done
##'-----------------------------------------------------------------------------------------#
