#!/bin/bash
#$ -cwd -V
#$ -pe smp 20
#$ -l h_vmem=120G
#$ -e ~/log/
#$ -o ~/log
source ~/.bash_profile

cuffquant -p 20 -v -o ../cuffquant/$2 \
          ../merged.gtf $1
cp ../cuffquant/$2/abundances.cxb ../cuffquant/abundances/$2'.cxb'
