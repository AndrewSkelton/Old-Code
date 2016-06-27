#!/bin/bash
#$ -cwd -V
#$ -pe smp 20
#$ -l h_vmem=120G
#$ -e ~/log/
#$ -o ~/log
source ~/.bash_profile

cuffdiff -v -p 20 -o ./Cuffdiff --max-bundle-frags 10000000 -L NOF,OA,HCH,MSCBM,HOB ../merged.gtf \
        Flowcell_A_2_NOF_filtered_noXY.cxb,Flowcell_A_5_NOF_filtered_noXY.cxb,Flowcell_A_7_NOF_filtered_noXY.cxb,\
        Flowcell_B_2_NOF_filtered_noXY.cxb,Flowcell_B_5_NOF_filtered_noXY.cxb,Flowcell_B_6_NOF_filtered_noXY.cxb \
        Flowcell_B_7_OA_filtered_noXY.cxb,Flowcell_B_8_OA_filtered_noXY.cxb,Flowcell_B_3_OA_filtered_noXY.cxb,\
        Flowcell_B_4_OA_filtered_noXY.cxb,Flowcell_A_8_OA_filtered_noXY.cxb,Flowcell_B_1_OA_filtered_noXY.cxb,\
        Flowcell_A_6_OA_filtered_noXY.cxb,Flowcell_A_3_OA_filtered_noXY.cxb,Flowcell_A_4_OA_filtered_noXY.cxb,\
        Flowcell_A_1_OA_filtered_noXY.cxb \
        hch_rep1_filtered_noXY.cxb,hch_rep2_filtered_noXY.cxb \
        mscbm_rep1_filtered_noXY.cxb,mscbm_rep2_filtered_noXY.cxb \
        hob_rep1_filtered_noXY.cxb,hob_rep2_filtered_noXY.cxb
