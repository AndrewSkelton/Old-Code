#-------------------------------------------------------------------------------------------#
#  Author      : Andrew J Skelton                                                           |
#  Language    : R Statistical Programming Language                                         |
#  Study       : RNA Seq Quantification Tool Compare                                        |
#  Data Owner  : Worked Example Using NOF OA RNA Seq Data - Prof. David Young               |
#  Description : RNA Seq Quantification Tool Compare - Kallisto                             |
#-------------------------------------------------------------------------------------------#



##'Set the Working Directory, load essential packages, and import files
##'-----------------------------------------------------------------------------------------#
source("http://bioconductor.org/biocLite.R")
biocLite()

setwd("~/Documents/Bioinformatics/Customers/Young/Test/kallisto/grch38/alignment/")

pheno_in <- read.table("../Pheno.txt",
                       sep="\t",
                       header=T,
                       row.names=1)
##'-----------------------------------------------------------------------------------------#



##'Get Normalised Counts - Kallisto - Gene Level
##'-----------------------------------------------------------------------------------------#
#Apparently Coming Soon
##'-----------------------------------------------------------------------------------------#



##'Get Normalised Counts - Kallisto - Transcript Level
##'-----------------------------------------------------------------------------------------#
df_kallisto <- matrix()

for(i in 1:length(list.files())) {
  sample_in  <- read.table(list.files()[i], sep="\t", header=T)
  sample_tmp <- sample_in[,c(1,4)]
  colnames(sample_tmp)[2] <- list.files()[i]
  df_kallisto <- cbind(df_kallisto, sample_tmp)
}

df_kallisto                          <- df_kallisto[,c(2,seq(3,ncol(df_kallisto),2))]
rownames(df_kallisto)                <- df_kallisto$target_id
Kallisto_Transcript_Counts           <- df_kallisto[,-1]
colnames(Kallisto_Transcript_Counts) <- pheno_in$Sample
##'-----------------------------------------------------------------------------------------#



##'Save Tables
##'-----------------------------------------------------------------------------------------#
save(Kallisto_Transcript_Counts, file="Kallisto_Data.Rdata")
##'-----------------------------------------------------------------------------------------#
