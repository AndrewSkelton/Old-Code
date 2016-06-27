#!/bin/bash
tophat -o ${1} \
       --phred64-quals \
       -p 2 \
       --transcriptome-index ${3} \
       --mate-std-dev ${4} \
       --mate-inner-dist ${5} \
       --tmp-dir ${1}/tmp \
       --microexon-search \
       --report-secondary-alignments \
       --b2-very-sensitive \
       ${6} \
       ${7} \
       ${8}

mv ${1}/accepted_hits.bam ${9}/alignment/${10}.bam
mv ${1}/align_summary.txt ${9}/summary/${10}.txt
mv ${1}/deletions.bed ${9}/deletions/${10}.bed
mv ${1}/insertions.bed ${9}/insertions/${10}.bed
mv ${1}/junctions.bed ${9}/junctions/${10}.bed
mv ${1}/prep_reads.info ${9}/prep/${10}.info
mv ${1}/unmapped.bam ${9}/unmapped/${10}.bam

samtools index ${9}/alignment/${10}.bam ${9}/alignment/${10}.bam.bai

