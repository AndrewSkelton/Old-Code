#!/bin/bash

declare -a arr=("A" "B")
RAW_DATA_PATH="/home/andrew/Young/Pipeline_Analysis/raw_data"
KALLISTO_INDEX="/data/genomes/GRCh38_Ensembl_Homo_Sapiens/hs_grch38_kallisto_index.idx"
KALLISTO_OUT_DIR="/home/andrew/Young/Pipeline_Analysis/kallisto/grch38/"
SCRIPT_PATH="/home/andrew/Young/Pipeline_Analysis/scripts/Sulaco"
#SALMON_GENE_MAP="/data/genomes/GRCh38_Ensembl_Homo_Sapiens/hs_grch38_salmon_gene_map.txt"

mkdir -p ${KALLISTO_OUT_DIR}/alignment

for i in "${arr[@]}"
do
  for j in {1..8}
  do
    FORWARD_IN='Flowcell_'${i}'_'${j}'_Forward.fastq.gz'
    REVERSE_IN='Flowcell_'${i}'_'${j}'_Reverse.fastq.gz'
    OUT_PREFIX='Flowcell_'${i}'_'${j}

    mkdir -p ${KALLISTO_OUT_DIR}/${OUT_PREFIX}

    sh kallisto_r.sh \
              ${KALLISTO_INDEX} \
              ${RAW_DATA_PATH}/${FORWARD_IN} \
              ${RAW_DATA_PATH}/${REVERSE_IN} \
              ${KALLISTO_OUT_DIR}/${OUT_PREFIX} \
              ${KALLISTO_OUT_DIR}/alignment \
              ${OUT_PREFIX} &

  done
done
