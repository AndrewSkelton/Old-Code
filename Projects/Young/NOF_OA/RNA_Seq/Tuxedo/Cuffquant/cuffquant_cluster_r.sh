#!/bin/bash

declare -a arr=("A" "B")
RAW_DATA_PATH="/home/nas151/WORKING_DATA/tophat/tophat_grch38/alignment"
GTF_ANNO="/home/nas151/WORKING_DATA/cuffmerge/grch38_No_Novel/merged.gtf"
FA_ANNO="/opt/databases/genomes/GRCh38_Ensembl_Homo_Sapiens/hs_grch38.fa"
CUFFQUANT_OUT_DIR="/home/nas151/WORKING_DATA/cuffquant/grch38_No_Novel"

mkdir -p ${CUFFQUANT_OUT_DIR}/abundances

for i in "${arr[@]}"
do
  for j in {1..8}
  do
    OUT_PREFIX='Flowcell_'${i}'_'${j}
    mkdir -p ${CUFFQUANT_OUT_DIR}/${OUT_PREFIX}
    qsub cuffquant.sh \
                    ${CUFFQUANT_OUT_DIR}/${OUT_PREFIX} \
                    ${CUFFQUANT_OUT_DIR} \
                    ${OUT_PREFIX} \
                    ${FA_ANNO} \
                    ${GTF_ANNO} \
                    ${RAW_DATA_PATH}/${OUT_PREFIX}'.bam' &
  done
done
