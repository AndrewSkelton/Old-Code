#-------------------------------------------------------------------------------------------#
#  Author      : Andrew J Skelton                                                           |
#  Language    : R Statistical Programming Language                                         |
#  Study       : Utility Code                                                               |
#  Data Owner  : NA                                                                         |
#  Description : Perform a query of dbSNP Identifiers                                       |
#-------------------------------------------------------------------------------------------#



##'Set the Working Directory, load essential packages, and import files
##'-----------------------------------------------------------------------------------------#
library(biomaRt)
snp.db     <- useMart("snp", dataset="hsapiens_snp")
##'-----------------------------------------------------------------------------------------#



##'Set the Working Directory, load essential packages, and import files
##'-----------------------------------------------------------------------------------------#
dbSNP_ID_In <- c("rs5931272")
nt.biomart  <- getBM(c("refsnp_id","allele","chr_name","chrom_start",
                       "chrom_strand","associated_gene",
                       "ensembl_gene_stable_id"),
                     filters="snp_filter",
                     values=dbSNP_ID_In,
                     mart=snp.db)
##'-----------------------------------------------------------------------------------------#
