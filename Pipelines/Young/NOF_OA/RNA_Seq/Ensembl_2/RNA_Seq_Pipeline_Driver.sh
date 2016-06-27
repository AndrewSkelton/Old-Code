#!/bin/bash

mkdir -p ./Tophat/accepted_hits ./Cufflinks/transcripts ./Cuffmerge ./Cuffdiff

##Standard Deviation and Mate inner distances for each sample.
###Possible Expansion : Auto calculation of this (Bowtie2 / Piccard Tools)
declare -a arr=("A" "B")
declare -a mate_inner_dist=(29 35 28 28 40 39 51 43 36 42 37 41 49 27 53 54 )
declare -a std_dev=(33 35 38 42 39 36 34 33 36 42 39 61 54 41 50 48 )
pointer=0

##bt2 Index Build - Option
##Would require a common folder with a common name for reference gtf, fa
##Output would need to be specified.

##Transcriptome Index Build - Option
##Would require a common folder with a common name for reference gtf, fa and bt2 index

##Tophat and Cufflinks calls
for i in "${arr[@]}"
do
        for j in {1..8}
        do
                i1='Flowcell_'$i'_'$j'_Forward.fastq'
                i2='Flowcell_'$i'_'$j'_Reverse.fastq'
                dir='Flowcell_'$i'_'$j
                tophat_name='Flowcell_'$i'_'$j'_Tophat'

                qsub -N $tophat_name ./TOPHAT_RUNNER.sh $i1 $i2 ${mate_inner_dist[$pointer]} ${std_dev[$pointer]} $dir
                qsub -N 'Cufflinks' -hold_jid $tophat_name ./CUFFLINKS_RUNNER.sh $dir

                pointer=$pointer+1
        done
done

##Cuffmerge followed by cuffdiff
qsub -N Cuffmerge -hold_jid Cufflinks ./cuffmerge_runner.sh
qsub -N Cuffdiff  -hold_jid Cuffmerge ./cuffdiff_runner.sh
