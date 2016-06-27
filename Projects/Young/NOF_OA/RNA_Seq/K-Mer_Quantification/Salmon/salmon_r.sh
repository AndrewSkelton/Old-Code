#!/bin/bash

declare -a arr=("A" "B")
RAW_DATA_PATH="/home/andrew/Young/Pipeline_Analysis/raw_data"
SALMON_INDEX="/data/genomes/GRCh38_Ensembl_Homo_Sapiens/hs_grch38_salmon_idx"
SALMON_OUT_DIR="/home/andrew/Young/Pipeline_Analysis/salmon/grch38_gene"
SCRIPT_PATH="/home/andrew/Young/Pipeline_Analysis/scripts/"
SALMON_GENE_MAP="/data/genomes/GRCh38_Ensembl_Homo_Sapiens/hs_grch38_salmon_gene_map.txt"

mkdir -p ${SALMON_OUT_DIR}/alignment

for i in "${arr[@]}"
do
  for j in {1..8}
  do
    FORWARD_IN='Flowcell_'${i}'_'${j}'_Forward.fastq.gz'
    REVERSE_IN='Flowcell_'${i}'_'${j}'_Reverse.fastq.gz'
    OUT_PREFIX='Flowcell_'${i}'_'${j}

    mkdir -p ${SALMON_OUT_DIR}/${OUT_PREFIX}

    ${SCRIPT_PATH}/salmon.sh \
              		     ${RAW_DATA_PATH}/${FORWARD_IN} \
              		     ${RAW_DATA_PATH}/${REVERSE_IN} \
              		     ${SALMON_INDEX} \
              		     ${SALMON_OUT_DIR}/${OUT_PREFIX} \
			     ${SALMON_GENE_MAP} &

  done
done

