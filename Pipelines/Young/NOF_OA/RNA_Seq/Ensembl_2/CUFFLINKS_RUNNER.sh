#!/bin/bash
#$ -cwd -V
#$ -pe smp 10
#$ -l h_vmem=50G
#$ -e ~/log/
#$ -o ~/log
source ~/.bash_profile

cufflinks -p 10 -q --max-bundle-frags 30000000 --GTF-guide /opt/databases/genomes/Ensembl_2/genome.gtf \
		  -o ./Cufflinks/$1 ./Tophat/accepted_hits/$1'.bam'

mv ./Cufflinks/$1/transcripts.gtf ./Cufflinks/transcripts/$1'.gtf'
