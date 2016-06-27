#!/bin/bash

declare -a arr=("A" "B")
RAW_DATA_PATH="/home/andrew/Young/Pipeline_Analysis/tophat/tophat_grch38/alignment"
GTF_ANNO="/home/andrew/Young/Pipeline_Analysis/cuffmerge/grch38/merged.gtf"
FA_ANNO="/data/genomes/GRCh38_Ensembl_Homo_Sapiens/hs_grch38.fa"
CUFFQUANT_OUT_DIR="/home/andrew/Young/Pipeline_Analysis/cuffquant/grch38"
SCRIPT_PATH="/home/andrew/Young/Pipeline_Analysis/scripts"
mkdir -p ${CUFFQUANT_OUT_DIR}/abundances

for i in "${arr[@]}"
do
  for j in {1..8}
  do
    OUT_PREFIX='Flowcell_'${i}'_'${j}
    mkdir -p ${CUFFQUANT_OUT_DIR}/${OUT_PREFIX}
    ${SCRIPT_PATH}/cuffquant_sulaco.sh \
                                      ${CUFFQUANT_OUT_DIR}/${OUT_PREFIX} \
                                      ${CUFFQUANT_OUT_DIR} \
                                      ${OUT_PREFIX} \
                                      ${FA_ANNO} \
                                      ${GTF_ANNO} \
                                      ${RAW_DATA_PATH}/${OUT_PREFIX}'.bam' &
  done
done

