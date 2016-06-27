#!/bin/bash
#$ -cwd -V
#$ -pe smp 20
#$ -l h_vmem=100G
#$ -e ~/log/
#$ -o ~/log
source ~/.bash_profile

cufflinks -p 20 -v -G /opt/databases/genomes/Ensembl_2/genome.gtf \
		  -o ./Cufflinks/$1 ./Tophat/accepted_hits/$1'.bam'

mv ./Cufflinks/$1/transcripts.gtf ./Cufflinks/transcripts/$1'.gtf'
