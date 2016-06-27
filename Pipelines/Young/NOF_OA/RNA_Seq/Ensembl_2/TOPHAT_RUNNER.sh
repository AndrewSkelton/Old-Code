#!/bin/bash
#$ -cwd -V
#$ -pe smp 10
#$ -l h_vmem=50G
#$ -e ~/log/
#$ -o ~/log
source ~/.bash_profile

tophat -r $3 --mate-std-dev $4 --phred64-quals --tmp-dir $TMPDIR -p 10 --transcriptome-index /opt/databases/genomes/Ensembl_2/transcriptome_index/Ensembl \
		--b2-very-sensitive -o ./Tophat/$5 -G /opt/databases/genomes/Ensembl_2/genome.gtf \
		/opt/databases/genomes/Ensembl_2/genome	'./trimmed_raw_data/'$1 './trimmed_raw_data/'$2

mv ./Tophat/$5/accepted_hits.bam ./Tophat/accepted_hits/$5'.bam'
samtools index ./Tophat/accepted_hits/$5'.bam'
