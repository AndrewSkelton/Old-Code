#!/bin/bash

mkdir -p ./Cufflinks/transcripts ./Tophat/accepted_hits

##Standard Deviation and Mate inner distances for each sample.
###Possible Expansion : Auto calculation of this (Bowtie2 / Piccard Tools)
declare -a arr=("1" "2")
declare -a mate_inner_dist=(294 281)
declare -a std_dev=(52 51)
pointer=0

##bt2 Index Build - Option
##Would require a common folder with a common name for reference gtf, fa
##Output would need to be specified.

##Transcriptome Index Build - Option
##Would require a common folder with a common name for reference gtf, fa and bt2 index

##Tophat and Cufflinks calls
for i in "${arr[@]}"
do
      'hch_rep1_rd1.fastq.gz'
      i1='mscbm_rep'$i'_rd1.fastq.gz'
      i2='mscbm_rep'$i'_rd2.fastq.gz'
      dir='mscbm_rep'$i
      tophat_name='hch_rep'$i'_Tophat'

      qsub -N $tophat_name ./TOPHAT_RUNNER.sh $i1 $i2 ${mate_inner_dist[$pointer]} ${std_dev[$pointer]} $dir
      qsub -N 'Cufflinks'$dir -hold_jid $tophat_name ./CUFFLINKS_RUNNER.sh $dir

      pointer=$pointer+1

done

##Cuffmerge followed by cuffdiff
#qsub -N Cuffmerge_2 -hold_jid Cufflinks_2 ./cuffmerge_runner.sh
#qsub -N Cuffdiff_2  -hold_jid Cuffmerge_2 ./cuffdiff_runner.sh
