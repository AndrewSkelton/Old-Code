#!/bin/bash

#-------------------------------------------------------------------------------------------#
#  Author      : Andrew J Skelton                                                           |
#  Language    : Bash                                                                       |
#  Study       : Quantification of SRX544934 from ENA                                       |
#  Data Owner  : ENA - At the request of Colin Shepard                                      |
#  Description : Runs the Salmon Quantification process															        |
#-------------------------------------------------------------------------------------------#



##'Run Salmon Qunatification
##'Parameters:
##'						-g - Gene Map
##'						-p - Number of threads
##'						-i - Salmon Index
##'						-l - Library Type (IU - Inward Unstranded)
##'						-1 - Fastq File (Gunzipped)
##'					  -2 - Fastq File (Gunzipped)
##'						-o - Output Directory
##'-----------------------------------------------------------------------------------------#
salmon quant \
	    			-g ${5} \
            -p 5 \
            -i ${3} \
            -l 'IU' \
            -1 <(gunzip -c ${1}) \
            -2 <(gunzip -c ${2}) \
            -o ${4} &
##'-----------------------------------------------------------------------------------------#
