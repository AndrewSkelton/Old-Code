#!/bin/bash

#-------------------------------------------------------------------------------------------#
#  Author      : Andrew J Skelton                                                           |
#  Language    : Bash                                                                       |
#  Study       : NOF and OA RNA Seq Study                                                   |
#  Data Owner  : Newcastle University - Prof. David Young                                   |
#  Description : Extract Fasta Records for TCONS of a Tuxedo Pipeline Run                   |
#-------------------------------------------------------------------------------------------#



##'Extract for GRCh38
##'-----------------------------------------------------------------------------------------#
files=( "TCONS_00262175" "TCONS_00262176" "TCONS_00262177" "TCONS_00262178" "TCONS_00262179" \
        "TCONS_00262180" "TCONS_00262181" "TCONS_00262182" "TCONS_00262183" "TCONS_00262184" )
touch grch38_TCON_fa_Extract.fa
for i in "${files[@]}"
do
  grep ${i} grch38_merged.gtf > tmp.gtf
  gffread -w tmp.fa -g hs_grch38.fa tmp.gtf
  cat tmp.fa >> grch38_TCON_fa_Extract.fa
  rm tmp.*
done
##'-----------------------------------------------------------------------------------------#



##'Extract for GRCh37
##'-----------------------------------------------------------------------------------------#
files=( "TCONS_00260503" "TCONS_00260504" "TCONS_00260505" "TCONS_00260506" "TCONS_00260507" \
        "TCONS_00260508" "TCONS_00260509" "TCONS_00260510" )
touch grch37_TCON_fa_Extract.fa
for i in "${files[@]}"
do
  grep ${i} grch37_merged.gtf > tmp.gtf
  gffread -w tmp.fa -g hs_grch37.fa tmp.gtf
  cat tmp.fa >> grch37_TCON_fa_Extract.fa
  rm tmp.*
done
##'-----------------------------------------------------------------------------------------#
