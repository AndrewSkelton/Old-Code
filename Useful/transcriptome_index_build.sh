#!/bin/bash
#$ -cwd -V
#$ -pe smp 1
#$ -l h_vmem=70G
#$ -e ~/log/
#$ -o ~/log
source ~/.bash_profile
tophat -G /opt/databases/genomes/Ensembl/Homo_sapiens.GRCh37.Ensembl.gtf --transcriptome-index=transcriptome_index/Ensembl /opt/databases/genomes/Ensembl/Homo_sapiens.GRCh37.Ensembl
