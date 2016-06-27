#!/bin/bash
#$ -cwd -V
#$ -pe smp 5
#$ -l h_vmem=20G
#$ -e ~/log/
#$ -o ~/log
source ~/.bash_profile

cuffquant -o ${1} \
          -p 5 \
          -b ${4} \
          -u \
          --max-bundle-frags 500000000 \
          ${5} \
          ${6}

mv ${1}/abundances.cxb ${2}/abundances/${3}.cxb
