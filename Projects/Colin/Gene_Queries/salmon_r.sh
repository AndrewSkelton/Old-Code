#!/bin/bash

#-------------------------------------------------------------------------------------------#
#  Author      : Andrew J Skelton                                                           |
#  Language    : Bash                                                                       |
#  Study       : Quantification of SRX544934 from ENA                                       |
#  Data Owner  : ENA - At the request of Colin Shepard                                      |
#  Description : Launches the Salmon Quantification process													        |
#-------------------------------------------------------------------------------------------#



##'Set Variables
##'-----------------------------------------------------------------------------------------#
declare -a arr=("SRR364315" "SRR364316")
RAW_DATA_PATH="/home/andrew/Colin/SRX544934/raw_data"
SALMON_INDEX="/data/genomes/GRCh38_Ensembl_Homo_Sapiens/hs_grch38_salmon_idx"
SALMON_OUT_DIR="/home/andrew/Colin/SRX544934/Salmon"
SCRIPT_PATH="/home/andrew/Colin/SRX544934/Scripts"
SALMON_GENE_MAP="/data/genomes/GRCh38_Ensembl_Homo_Sapiens/hs_grch38_salmon_gene_map.txt"
##'-----------------------------------------------------------------------------------------#



##'Launch Salmon Qunatification
##'-----------------------------------------------------------------------------------------#
mkdir -p ${SALMON_OUT_DIR}/alignment_GRCh38

for i in "${arr[@]}"
do
  FORWARD_IN=${i}'_1.fastq.gz'
  REVERSE_IN=${i}'_2.fastq.gz'
  OUT_PREFIX=${i}

  mkdir -p ${SALMON_OUT_DIR}/${OUT_PREFIX}

  ${SCRIPT_PATH}/salmon.sh \
            		     ${RAW_DATA_PATH}/${FORWARD_IN} \
            		     ${RAW_DATA_PATH}/${REVERSE_IN} \
            		     ${SALMON_INDEX} \
            		     ${SALMON_OUT_DIR}/${OUT_PREFIX} \
  	                 ${SALMON_GENE_MAP} &
done
##'-----------------------------------------------------------------------------------------#
