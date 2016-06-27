#-------------------------------------------------------------------------------------------#
#  Author      : Andrew J Skelton                                                           |
#  Language    : R Statistical Programming Language                                         |
#  Study       : RNA Seq Quantification Tool Compare                                        |
#  Data Owner  : Worked Example Using NOF OA RNA Seq Data - Prof. David Young               |
#  Description : RNA Seq Quantification Tool Compare - Tuxedo                               |
#-------------------------------------------------------------------------------------------#



##'Set the Working Directory, load essential packages, and import files
##'-----------------------------------------------------------------------------------------#
source("http://bioconductor.org/biocLite.R")
biocLite()

library(cummeRbund)

setwd("~/Documents/Bioinformatics/Customers/Young/Test/grch38_No_Novel_v2/")

xloc_map             <- read.table("genes.fpkm_tracking", sep="\t", header=T)
tcons_map            <- read.table("isoforms.fpkm_tracking", sep="\t", header=T)
cuff                 <- readCufflinks()
tuxedo_map           <- tcons_map[,c(1,3,4)]
##'-----------------------------------------------------------------------------------------#




##'Get Normalised Counts - Tuxedo - Gene Level
##'-----------------------------------------------------------------------------------------#
Tuxedo_Gene_Counts       <- repCountMatrix(genes(cuff))
# gene.counts.fpkm      <- repFpkmMatrix(genes(cuff))
##'-----------------------------------------------------------------------------------------#



##'Get Normalised Counts - Tuxedo - Transcript Level
##'-----------------------------------------------------------------------------------------#
Tuxedo_Transcript_Counts <- repCountMatrix(isoforms(cuff))
# isoform.counts.fpkm   <- repFpkmMatrix(isoforms(cuff))
##'-----------------------------------------------------------------------------------------#



##'Save Tables
##'-----------------------------------------------------------------------------------------#
save(Tuxedo_Gene_Counts, Tuxedo_Transcript_Counts, tuxedo_map, 
     file="Tuxedo_Data.Rdata")
##'-----------------------------------------------------------------------------------------#
