#!/bin/bash
#$ -cwd -V
#$ -pe smp 20
#$ -l h_vmem=100G
#$ -e ~/log/
#$ -o ~/log
source ~/.bash_profile

cuffmerge -p 20 -s /opt/databases/genomes/Ensembl_2/genome.fa \
          -g /opt/databases/genomes/Ensembl_2/genome.gtf \
          -o ./Cuffmerge ./assemblies.txt
