#!/bin/bash

#1 - output directory
#2 - GTF
#3 - transcriptome index
#4 - std dev
#5 - inner distance
#6 - bt2 index
#7 - forward reads
#8 - reverse reads

declare -a arr=("A" "B")
RAW_DATA_PATH="/home/andrew/Young/Pipeline_Analysis/raw_data"
METRICS="/home/andrew/Young/Pipeline_Analysis/raw_data/Experiment_Metrics.log"

GTF_ANNO="/data/genomes/GRCh38_Ensembl_Homo_Sapiens/hs_grch38.gtf"
TRANSCRIPTOME_INDEX="/data/genomes/GRCh38_Ensembl_Homo_Sapiens/hs_grch38_transcriptome_idx/hs_grch38_transcriptome"
BT2_INDEX="/data/genomes/GRCh38_Ensembl_Homo_Sapiens/hs_grch38_bt2_idx/hs_grch38_bt2_idx"
FA_ANNO="/data/genomes/GRCh38_Ensembl_Homo_Sapiens/hs_grch38.fa"

TOPHAT_OUT_DIR="/home/andrew/Young/Pipeline_Analysis/tophat/tophat_grch38"
#CUFFLINKS_OUT_DIR="/home/nas151/WORKING_DATA/NOF_OA_Analysis_Master_Branch/cufflinks/cufflinks_grch38"

mkdir -p ${TOPHAT_OUT_DIR}/alignment \
${TOPHAT_OUT_DIR}/junctions \
${TOPHAT_OUT_DIR}/summary \
${TOPHAT_OUT_DIR}/insertions \
${TOPHAT_OUT_DIR}/deletions \
${TOPHAT_OUT_DIR}/prep \
${TOPHAT_OUT_DIR}/unmapped \
${TOPHAT_OUT_DIR}/alignment \
#${CUFFLINKS_OUT_DIR}/transcripts \
#${CUFFLINKS_OUT_DIR}/isoform_tracking \
#${CUFFLINKS_OUT_DIR}/gene_tracking \
#${CUFFLINKS_OUT_DIR}/skipped \

for i in "${arr[@]}"
do
  for j in {1..8}
  do
    FORWARD_IN='Flowcell_'${i}'_'${j}'_Forward.fastq.gz'
    REVERSE_IN='Flowcell_'${i}'_'${j}'_Reverse.fastq.gz'
    OUT_PREFIX='Flowcell_'${i}'_'${j}

    TOPHAT_JOB_NAME='Flowcell_'${i}'_'${j}'_TOPHAT_GRCH378'
    CUFFLINKS_JOB_NAME='Flowcell_'${i}'_'${j}'_CUFFLINKS_GRCH378'

    SAMPLE_METRICS=`grep ${OUT_PREFIX} ${METRICS}`
    INNER_DIST=`grep ${OUT_PREFIX} ${METRICS} | cut -f 2 -`
    STD_DEV=`grep ${OUT_PREFIX} ${METRICS} | cut -f 3 -`

    mkdir -p ${TOPHAT_OUT_DIR}/${OUT_PREFIX}
    #mkdir -p ${CUFFLINKS_OUT_DIR}/${OUT_PREFIX}

    sh ./tophat_sulaco.sh \
                              ${TOPHAT_OUT_DIR}/${OUT_PREFIX} \
                              ${GTF_ANNO} \
                              ${TRANSCRIPTOME_INDEX} \
                              ${STD_DEV} \
                              ${INNER_DIST} \
                              ${BT2_INDEX} \
                              ${RAW_DATA_PATH}/${FORWARD_IN} \
                              ${RAW_DATA_PATH}/${REVERSE_IN} \
                              ${TOPHAT_OUT_DIR} \
                              ${OUT_PREFIX} &

  done
done

