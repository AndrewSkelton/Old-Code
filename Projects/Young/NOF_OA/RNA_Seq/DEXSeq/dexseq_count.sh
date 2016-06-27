#!/bin/bash

#-------------------------------------------------------------------------------------------#
#  Author      : Andrew J Skelton                                                           |
#  Language    : Bash                                                                       |
#  Study       : NOF and OA RNA Seq Study                                                   |
#  Data Owner  : Newcastle University - Prof. David Young                                   |
#  Description : Use DEXSeq_Count to count exon level features sample by sample using       |
#                GRCh38 annotation                                                          |
#-------------------------------------------------------------------------------------------#



##'Set script variables
##'-----------------------------------------------------------------------------------------#
IN_DIR="/home/andrew/Young/Pipeline_Analysis/tophat/tophat_grch38/alignment/sorted"
REF_GFF="/data/genomes/GRCh38_Ensembl_Homo_Sapiens/hs_grch38_dexseq_ref/hs_grch38_dexseq.gff"
SCRIPT_IN="/data/genomes/GRCh38_Ensembl_Homo_Sapiens/hs_grch38_dexseq_ref/dexseq_count.py"
COUNT_OUT="/home/andrew/Young/Pipeline_Analysis/dexseq/grch38"
##'-----------------------------------------------------------------------------------------#



##'Loop through filenames and set a child process of DEXSeq_Count running
##'-----------------------------------------------------------------------------------------#
declare -a arr=("A" "B")
for i in "${arr[@]}"
do
  for j in {1..8}
  do
    OUT_PREFIX='Flowcell_'${i}'_'${j}
    python ${SCRIPT_IN} -p yes \
                        -s no \
                        -f bam \
                        -r name \
                        ${REF_GFF} \
                        ${IN_DIR}/${OUT_PREFIX}'.sorted.bam' \
                        ${COUNT_OUT}/${OUT_PREFIX}'.count' &
  done
done
##'-----------------------------------------------------------------------------------------#
