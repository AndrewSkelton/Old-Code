#!/bin/bash

RAW_DATA_PATH="/home/andrew/Young/Pipeline_Analysis/cuffquant/grch38/abundances"
GTF_ANNO="/home/andrew/Young/Pipeline_Analysis/cuffmerge/grch38/merged.gtf"
FA_ANNO="/data/genomes/GRCh38_Ensembl_Homo_Sapiens/hs_grch38.fa"
CUFFDIFF_OUT_DIR="/home/andrew/Young/Pipeline_Analysis/cuffdiff/grch38"
mkdir -p ${CUFFDIFF_OUT_DIR}

cuffdiff -o ${CUFFDIFF_OUT_DIR} \
         -L NOF,OA \
         -p 20 \
         -b ${FA_ANNO} \
         -u \
         --max-bundle-frags 500000000 \
         --library-norm-method geometric \
         ${GTF_ANNO} \
         ${RAW_DATA_PATH}/Flowcell_A_2.cxb,${RAW_DATA_PATH}/Flowcell_A_5.cxb,${RAW_DATA_PATH}/Flowcell_A_7.cxb,${RAW_DATA_PATH}/Flowcell_B_2.cxb,${RAW_DATA_PATH}/Flowcell_B_5.cxb,${RAW_DATA_PATH}/Flowcell_B_6.cxb \
         ${RAW_DATA_PATH}/Flowcell_A_1.cxb,${RAW_DATA_PATH}/Flowcell_A_3.cxb,${RAW_DATA_PATH}/Flowcell_A_4.cxb,${RAW_DATA_PATH}/Flowcell_A_6.cxb,${RAW_DATA_PATH}/Flowcell_A_8.cxb,${RAW_DATA_PATH}/Flowcell_B_1.cxb,${RAW_DATA_PATH}/Flowcell_B_3.cxb,${RAW_DATA_PATH}/Flowcell_B_4.cxb,${RAW_DATA_PATH}/Flowcell_B_7.cxb,${RAW_DATA_PATH}/Flowcell_B_8.cxb

