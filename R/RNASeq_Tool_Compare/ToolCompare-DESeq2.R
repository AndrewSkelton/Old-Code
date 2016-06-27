#-------------------------------------------------------------------------------------------#
#  Author      : Andrew J Skelton                                                           |
#  Language    : R Statistical Programming Language                                         |
#  Study       : RNA Seq Quantification Tool Compare                                        |
#  Data Owner  : Worked Example Using NOF OA RNA Seq Data - Prof. David Young               |
#  Description : RNA Seq Quantification Tool Compare - DESeq2                               |
#-------------------------------------------------------------------------------------------#



##'Set the Working Directory, load essential packages, and import files
##'-----------------------------------------------------------------------------------------#
source("http://bioconductor.org/biocLite.R")
biocLite()

library(DESeq2)

setwd("~/Documents/Bioinformatics/Customers/Young/Test/htseq_count/grch38/")
##'-----------------------------------------------------------------------------------------#



##'DESeq2 - Get Normalised Counts
##'-----------------------------------------------------------------------------------------#
des                <- DESeqDataSetFromHTSeqCount(pheno_in,
                                                 directory="./",
                                                 design= ~ Treatment)
des                <- DESeq(des)
DESeq2_Gene_Counts <- counts(des, normalized=T)
##'-----------------------------------------------------------------------------------------#



##'Save Normalised Counts to save re-processing time
##'-----------------------------------------------------------------------------------------#
save(DESeq2_Gene_Counts, file="DESeq2_Data.Rdata")
##'-----------------------------------------------------------------------------------------#
