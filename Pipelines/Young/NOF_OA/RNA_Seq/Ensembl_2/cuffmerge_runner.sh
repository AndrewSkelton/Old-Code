#!/bin/bash
#$ -cwd -V
#$ -pe smp 1
#$ -l h_vmem=80G
#$ -e ~/log/
#$ -o ~/log
source ~/.bash_profile

find ./Cufflinks/transcripts -iname '*.gtf' > assemblies.txt
cuffmerge -s /opt/databases/genomes/Ensembl_2/genome.fa -g /opt/databases/genomes/Ensembl_2/genome.gtf -o ./Cuffmerge ./assemblies.txt
