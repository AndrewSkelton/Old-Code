#!/bin/bash

#-------------------------------------------------------------------------------------------#
#  Author      : Andrew J Skelton                                                           |
#  Language    : Bash                                                                       |
#  Study       : miRNA Sequencing of Pooled Bovine Cartilage Samples                        |
#  Data Owner  : Newcastle University - Prof. David Young, (Cardiff Collaboration           |
#  Description : Decompress and Trim Adaptors, Discard Reads Less than 17bp                 |
#-------------------------------------------------------------------------------------------#



##'Decompress All Fastq Files (Assuming GunZip Compression)
##'-----------------------------------------------------------------------------------------#
for i in *.gz
do
	pigz -d -p 2 ${i} &
done
##'-----------------------------------------------------------------------------------------#



##'Wait for all Child Processes
##'-----------------------------------------------------------------------------------------#
wait
##'-----------------------------------------------------------------------------------------#



##'Run Cutadapt
##'-----------------------------------------------------------------------------------------#
for i in *.fastq
do
        ~/.local/bin/cutadapt \
											--minimum-length 17 \
											-a AGATCGGAAGAGCACACG \
											-o ../trimmed_$i \
											$i &
done
##'-----------------------------------------------------------------------------------------#
