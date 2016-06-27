#!/bin/bash

#-------------------------------------------------------------------------------------------#
#  Author      : Andrew J Skelton                                                           |
#  Language    : Bash                                                                       |
#  Study       : miRNA Sequencing of Pooled Bovine Cartilage Samples                        |
#  Data Owner  : Newcastle University - Prof. David Young, (Cardiff Collaboration           |
#  Description : Bowtie2 Runner - Sulaco                                                    |
#-------------------------------------------------------------------------------------------#



##'Run Bowtie2 with input parameters
##'-----------------------------------------------------------------------------------------#
bowtie2 -p 2 \
        --very-sensitive \
        -x ${3} \
        -U ${1} \
        -S ${4}/${2}'.sam' 2> ${4}/${2}.log

samtools view -bSh ${4}/${2}'.sam' > ${4}/${2}'.bam'
rm ${4}/${2}'.sam'

samtools sort ${4}/${2}'.bam' ${4}/${2}'.sorted'
rm ${4}/${2}'.bam'

samtools index ${4}/${2}'.sorted.bam'
mv ${2}/${2}'.sorted.bam*' ${5}
##'-----------------------------------------------------------------------------------------#
