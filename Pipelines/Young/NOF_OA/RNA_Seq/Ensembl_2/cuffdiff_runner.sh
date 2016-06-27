#!/bin/bash
#$ -cwd -V
#$ -pe smp 20
#$ -l h_vmem=80G
#$ -e ~/log/
#$ -o ~/log
source ~/.bash_profile

cuffdiff -q -p 20 -u -b /data/genomes/grch37/grch37.fa -o ./Cuffdiff ./Cuffmerge/merged.gtf \
Tophat/accepted_hits/Flowcell_A_2.bam,Tophat/accepted_hits/Flowcell_A_5.bam,Tophat/accepted_hits/Flowcell_A_7.bam,Tophat/accepted_hits/Flowcell_B_2.bam,Tophat/accepted_hits/Flowcell_B_5.bam,Tophat/accepted_hits/Flowcell_B_6.bam \
Tophat/accepted_hits/Flowcell_A_1.bam,Tophat/accepted_hits/Flowcell_A_3.bam,Tophat/accepted_hits/Flowcell_A_4.bam,Tophat/accepted_hits/Flowcell_A_6.bam,Tophat/accepted_hits/Flowcell_A_8.bam,\
Tophat/accepted_hits/Flowcell_B_1.bam,Tophat/accepted_hits/Flowcell_B_3.bam,Tophat/accepted_hits/Flowcell_B_4.bam,Tophat/accepted_hits/Flowcell_B_7.bam,Tophat/accepted_hits/Flowcell_B_8.bam &> cuffdiff_Merged.log
