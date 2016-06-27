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
RAW_DATA_PATH="/home/nas151/WORKING_DATA/NOF_OA_Analysis_Master_Branch/raw_data"
METRICS="/home/nas151/WORKING_DATA/NOF_OA_Analysis_Master_Branch/bt2/metrics/Experiment_Metrics.log"

GTF_ANNO="/opt/databases/genomes/GRCh38_Ensembl_Homo_Sapiens/hs_grch38.gtf"
TRANSCRIPTOME_INDEX="/opt/databases/genomes/GRCh38_Ensembl_Homo_Sapiens/hs_grch38_transcriptome_idx/hs_grch38_transcriptome"
BT2_INDEX="/opt/databases/genomes/GRCh38_Ensembl_Homo_Sapiens/hs_grch38_bt2_idx/hs_grch38_bt2_idx"
FA_ANNO="/opt/databases/genomes/GRCh38_Ensembl_Homo_Sapiens/hs_grch38.fa"

TOPHAT_OUT_DIR="/sharedlustre/users/nas151/tophat/tophat_grch38"
CUFFLINKS_OUT_DIR="/sharedlustre/users/nas151/cufflinks/cufflinks_grch38_no_novel"

mkdir -p ${CUFFLINKS_OUT_DIR}/transcripts \
         ${CUFFLINKS_OUT_DIR}/isoform_tracking \
         ${CUFFLINKS_OUT_DIR}/gene_tracking \
         ${CUFFLINKS_OUT_DIR}/skipped \

for i in "${arr[@]}"
do
  for j in {1..8}
  do
    CUFFLINKS_JOB_NAME='Flowcell_'${i}'_'${j}'_CUFFLINKS_GRCH38'
    OUT_PREFIX='Flowcell_'${i}'_'${j}
    mkdir -p ${CUFFLINKS_OUT_DIR}/${OUT_PREFIX}
    
    qsub -N ${CUFFLINKS_JOB_NAME} \
            cufflinks_grch38_cluster_r.sh \
                                         ${GTF_ANNO} \
                                         ${FA_ANNO} \
                                         ${OUT_PREFIX} \
                                         ${CUFFLINKS_OUT_DIR}/${OUT_PREFIX} \
                                         ${TOPHAT_OUT_DIR}/alignment/${OUT_PREFIX}'.bam' \
                                         ${CUFFLINKS_OUT_DIR} \
                                         ${OUT_PREFIX}

  done
done
