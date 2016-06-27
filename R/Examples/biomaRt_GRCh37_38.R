#-------------------------------------------------------------------------------------------#
#  Author      : Andrew J Skelton                                                           |
#  Language    : R Statistical Programming Language                                         |
#  Study       : Utility Code                                                               |
#  Data Owner  : NA                                                                         |
#  Description : Perform a query of GRCh37 or GRCh38                                        |
#-------------------------------------------------------------------------------------------#



##'Set the Working Directory, load essential packages, and import files
##'-----------------------------------------------------------------------------------------#
library(biomaRt)
##'-----------------------------------------------------------------------------------------#



##'Query GRCh38
##'-----------------------------------------------------------------------------------------#
ensembl_genes_in <- c("ENSG...", "ENSG...")

ensembl    <- useMart("ensembl",
                      dataset="hsapiens_gene_ensembl")

annotation <- getBM(attributes=c("ensembl_gene_id",
                                 "external_gene_name",
                                 "description"),
                    filters="ensembl_gene_id",
                    values=ensembl_genes_in,
                    ensembl)
##'-----------------------------------------------------------------------------------------#



##'Query GRCh37
##'-----------------------------------------------------------------------------------------#
ensembl_genes_in <- c("ENSG...", "ENSG...")

ensembl    <- useMart("ensembl",
                      dataset="hsapiens_gene_ensembl",
                      host="grch37.ensembl.org",
                      path="/biomart/martservice")

annotation <- getBM(attributes=c("ensembl_gene_id",
                                 "external_gene_name",
                                 "description"),
                    filters="ensembl_gene_id",
                    values=ensembl_genes_in,
                    ensembl)
##'-----------------------------------------------------------------------------------------#
