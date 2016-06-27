#-------------------------------------------------------------------------------------------#
#  Author      : Andrew J Skelton                                                           |
#  Language    : R Statistical Programming Language                                         |
#  Study       : RNA Seq Quantification Tool Compare                                        |
#  Data Owner  : Worked Example Using NOF OA RNA Seq Data - Prof. David Young               |
#  Description : RNA Seq Quantification Tool Compare - Salmon                               |
#-------------------------------------------------------------------------------------------#



##'Set the Working Directory, load essential packages, and import files
##'-----------------------------------------------------------------------------------------#
source("http://bioconductor.org/biocLite.R")
biocLite()

library(EBSeq)

setwd("~/Documents/Bioinformatics/Customers/Young/Test/salmon/grch38_gene")

pheno_in <- read.table("Pheno.txt",
                       sep="\t",
                       header=T,
                       row.names=1)
##'-----------------------------------------------------------------------------------------#


##'Utility Function - get_counts
##'                 - Description : Reads in a Salmon Output file and converts it to a
##'                                 useful dataframe structure
##'                 - Parameters  : file - path and filename of the Salmon .sf file
##'-----------------------------------------------------------------------------------------#
get_counts = function(file) {
  dat               <- read.table(file, skip=10, sep="\t", stringsAsFactors=F, row.names=1)
  raw_counts        <- dat[,4]
  names(raw_counts) <- rownames(dat)
  return(raw_counts)
}
##'-----------------------------------------------------------------------------------------#



##'Salmon - Gene Level - Get Normalised Counts
##'-----------------------------------------------------------------------------------------#
gene_files_in       <- list.files("./genes/",
                                  pattern="*.sf",
                                  full.names=T)
gene_counts         <- lapply(gene_files_in,
                              get_counts)
gene_count_table    <- matrix(nrow=length(gene_counts[[1]]),
                              ncol=length(gene_files_in))

for(i in 1:length(gene_counts)) {
  gene_count_table[,i] <- gene_counts[[i]]
}

rownames(gene_count_table) <- names(gene_counts[[1]])
colnames(gene_count_table) <- gsub(".sf", "", gene_files_in)
Sizes                      <- MedianNorm(gene_count_table)
Salmon_Gene_Counts         <- GetNormalizedMat(gene_count_table, Sizes)
##'-----------------------------------------------------------------------------------------#



##'Salmon - Transcript Level - Get Normalised Counts
##'-----------------------------------------------------------------------------------------#
transcript_files_in    <- list.files("../grch38/alignment/",
                                     pattern="*.sf",
                                     full.names=T)
transcript_counts      <- lapply(transcript_files_in,
                                 get_counts)
transcript_count_table <- matrix(nrow=length(transcript_counts[[1]]),
                                 ncol=length(transcript_files_in))

for(i in 1:length(transcript_counts)) {
  transcript_count_table[,i] <- transcript_counts[[i]]
}

rownames(transcript_count_table) <- names(transcript_counts[[1]])
colnames(transcript_count_table) <- gsub(".sf", "", transcript_files_in)
Sizes                            <- MedianNorm(transcript_count_table)
Salmon_Trans_Counts              <- GetNormalizedMat(transcript_count_table, Sizes)
##'-----------------------------------------------------------------------------------------#



##'Save Outputs to keep from processing again
##'-----------------------------------------------------------------------------------------#
save(pheno_in, Salmon_Gene_Counts, Salmon_Trans_Counts, file="Salmon_Data.Rdata")
##'-----------------------------------------------------------------------------------------#
