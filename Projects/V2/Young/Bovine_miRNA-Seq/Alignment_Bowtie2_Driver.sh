#!/bin/bash

#-------------------------------------------------------------------------------------------#
#  Author      : Andrew J Skelton                                                           |
#  Language    : Bash                                                                       |
#  Study       : miRNA Sequencing of Pooled Bovine Cartilage Samples                        |
#  Data Owner  : Newcastle University - Prof. David Young, (Cardiff Collaboration           |
#  Description : Bowtie2 Driver - Sulaco                                                    |
#-------------------------------------------------------------------------------------------#



##'Set script variables - Human Hairpin
##'-----------------------------------------------------------------------------------------#
RAW_DATA_PATH="/home/andrew/Paulina/Resequencing/trimmed"
OUT_DIR="/home/andrew/Paulina/Resequencing/April/Alignment/hairpin/human"
BT2_INDEX="/home/andrew/Paulina/Resequencing/April/bowtie2Idx/DNA/hairpin/human/human_hairpin_dna_bt2idx"
##'-----------------------------------------------------------------------------------------#

##'Run Bowtie2 - Human Hairpin
##'-----------------------------------------------------------------------------------------#
for i in ${RAW_DATA_PATH}/*.gz
do
    filename=$(basename "$i")
    extension="${filename##*.}"
    filename="${filename%.*}"
    filename="${filename%.*}"

    SAMPLE_IN=${i}
    OUT_PREFIX=${filename}

    sh bt2_sulaco.sh \
                    ${i} \
                    ${OUT_PREFIX} \
                    ${BT2_INDEX} \
                    ${OUT_DIR} &
done
##'-----------------------------------------------------------------------------------------#

##'Wait for Return of all child processes
##'-----------------------------------------------------------------------------------------#
wait
##'-----------------------------------------------------------------------------------------#



##'Set script variables - Human Mature
##'-----------------------------------------------------------------------------------------#
OUT_DIR="/home/andrew/Paulina/Resequencing/April/Alignment/mature/human"
BT2_INDEX="/home/andrew/Paulina/Resequencing/April/bowtie2Idx/DNA/mature/human/human_mature_dna_bt2idx"
##'-----------------------------------------------------------------------------------------#

##'Run Bowtie2 - Human Mature
##'-----------------------------------------------------------------------------------------#
for i in ${RAW_DATA_PATH}/*.gz
do
    filename=$(basename "$i")
    extension="${filename##*.}"
    filename="${filename%.*}"
    filename="${filename%.*}"

    SAMPLE_IN=${i}
    OUT_PREFIX=${filename}

    sh bt2_sulaco.sh \
                    ${i} \
                    ${OUT_PREFIX} \
                    ${BT2_INDEX} \
                    ${OUT_DIR} &
done
##'-----------------------------------------------------------------------------------------#

##'Wait for Return of all child processes
##'-----------------------------------------------------------------------------------------#
wait
##'-----------------------------------------------------------------------------------------#



##'Set script variables - Cow Hairpin
##'-----------------------------------------------------------------------------------------#
OUT_DIR="/home/andrew/Paulina/Resequencing/April/Alignment/hairpin/cow"
BT2_INDEX="/home/andrew/Paulina/Resequencing/April/bowtie2Idx/DNA/hairpin/cow/cow_hairpin_dna_bt2idx"
##'-----------------------------------------------------------------------------------------#

##'Run Bowtie2 - Cow Hairpin
##'-----------------------------------------------------------------------------------------#
for i in ${RAW_DATA_PATH}/*.gz
do
    filename=$(basename "$i")
    extension="${filename##*.}"
    filename="${filename%.*}"
    filename="${filename%.*}"

    SAMPLE_IN=${i}
    OUT_PREFIX=${filename}

    sh bt2_sulaco.sh \
                    ${i} \
                    ${OUT_PREFIX} \
                    ${BT2_INDEX} \
                    ${OUT_DIR} &
done
##'-----------------------------------------------------------------------------------------#

##'Wait for Return of all child processes
##'-----------------------------------------------------------------------------------------#
wait
##'-----------------------------------------------------------------------------------------#



##'Set script variables - Cow Mature
##'-----------------------------------------------------------------------------------------#
OUT_DIR="/home/andrew/Paulina/Resequencing/April/Alignment/mature/cow"
BT2_INDEX="/home/andrew/Paulina/Resequencing/April/bowtie2Idx/DNA/mature/cow/cow_mature_dna_bt2idx"
##'-----------------------------------------------------------------------------------------#

##'Run Bowtie2 - Cow Mature
##'-----------------------------------------------------------------------------------------#
for i in ${RAW_DATA_PATH}/*.gz
do
    filename=$(basename "$i")
    extension="${filename##*.}"
    filename="${filename%.*}"
    filename="${filename%.*}"

    SAMPLE_IN=${i}
    OUT_PREFIX=${filename}

    sh bt2_sulaco.sh \
                    ${i} \
                    ${OUT_PREFIX} \
                    ${BT2_INDEX} \
                    ${OUT_DIR} &
done
##'-----------------------------------------------------------------------------------------#

##'Wait for Return of all child processes
##'-----------------------------------------------------------------------------------------#
wait
##'-----------------------------------------------------------------------------------------#
