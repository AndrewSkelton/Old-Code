#!/bin/bash
#$ -cwd -V
#$ -pe smp 5
#$ -l h_vmem=20G
#$ -e ~/log/
#$ -o ~/log
source ~/.bash_profile

#1 - Reference GTF
#2 - Reference Fasta
#3 - Label
#4 - output dir
#5 - accepted hits bam file

cufflinks -p 5 \
          -G ${1} \
          -b ${2} \
          -u \
          --max-bundle-frags 500000000 \
          -L ${3} \
          -v \
          -o ${4} \
          ${5}

mv ${4}/transcripts.gtf ${6}/transcripts/${7}.gtf
mv ${4}/isoforms.fpkm_tracking ${6}/isoform_tracking/${7}.fpkm_tracking
mv ${4}/genes.fpkm_tracking ${6}/gene_tracking/${7}.fpkm_tracking
mv ${4}/skipped.gtf ${6}/skipped/${7}.gtf
