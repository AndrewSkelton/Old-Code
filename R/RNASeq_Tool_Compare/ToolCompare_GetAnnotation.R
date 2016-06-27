#-------------------------------------------------------------------------------------------#
#  Author      : Andrew J Skelton                                                           |
#  Language    : R Statistical Programming Language                                         |
#  Study       : RNA Seq Quantification Tool Compare                                        |
#  Data Owner  : Worked Example Using NOF OA RNA Seq Data - Prof. David Young               |
#  Description : RNA Seq Quantification Tool Compare - Get Annotation                       |
#-------------------------------------------------------------------------------------------#



##'Set the Working Directory, load essential packages, and import files
##'-----------------------------------------------------------------------------------------#
source("http://bioconductor.org/biocLite.R")
biocLite()

library(biomaRt)

setwd("~/Documents/Bioinformatics/Customers/Young/Test")

# load("Salmon_Data.Rdata")
##'-----------------------------------------------------------------------------------------#


##'Get Annotation for all Ensemble GRCh38 Known Transcripts
##'-----------------------------------------------------------------------------------------#
transcript_names <- rownames(Salmon_Trans_Counts)

ensembl          <- useMart("ensembl",
                            dataset="hsapiens_gene_ensembl")

annotation       <- getBM(attributes=c("ensembl_transcript_id",
                                       "ensembl_gene_id",
                                       "external_gene_name",
                                       "description"),
                          filters="ensembl_transcript_id",
                          values=transcript_names,
                          ensembl)

annotation       <- annotation[match(transcript_names,
                                     annotation$ensembl_transcript_id),]
##'-----------------------------------------------------------------------------------------#



##'Save Tables
##'-----------------------------------------------------------------------------------------#
save(annotation, file="Annotation_Data.Rdata")
##'-----------------------------------------------------------------------------------------#
