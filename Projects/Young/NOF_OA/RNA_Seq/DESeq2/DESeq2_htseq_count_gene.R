#-------------------------------------------------------------------------------------------#
#  Author      : Andrew J Skelton                                                           |
#  Language    : R Statistical Programming Language                                         |
#  Study       : NOF and OA Microarray Study                                                |
#  Data Owner  : Newcastle University - Prof. David Young                                   |
#  Description : Illumina HT-12v3 Microarray of NOF and OA Samples. Grouped analysis and    |
#                Regression by Age.
#-------------------------------------------------------------------------------------------#



##'Set the Working Directory, load essential packages, and import files
##'-----------------------------------------------------------------------------------------#
setwd("~/Young/Pipeline_Analysis/htseq_count/grch38/")

source("http://bioconductor.org/biocLite.R")
biocLite()

library(ggplot2)
library(DESeq2)
library(BiocParallel)
library(biomaRt)
##'-----------------------------------------------------------------------------------------#



##'Build pheno table and set worker threads
##'-----------------------------------------------------------------------------------------#
sampleTable <- data.frame(
  row.names    = c("F1_S1","F1_S2","F1_S3","F1_S4","F1_S5","F1_S6","F1_S7","F1_S8",
                   "F2_S1","F2_S2","F2_S3","F2_S4","F2_S5","F2_S6","F2_S7","F2_S8"),
  sample.names = c("F1_S1","F1_S2","F1_S3","F1_S4","F1_S5","F1_S6","F1_S7","F1_S8",
                   "F2_S1","F2_S2","F2_S3","F2_S4","F2_S5","F2_S6","F2_S7","F2_S8"),
  countName    = list.files("./", full.names=T, pattern="*.count"),
  condition    = c("OA","NOF","OA","OA","NOF","OA","NOF","OA","OA","NOF","OA","OA",
                   "NOF","NOF","OA","OA"),
  libType      = c(rep("paired-end", 16))
)
BPPARAM = MulticoreParam(workers=20)



##'Build DESeq2 Dataset/Model, and run DESeq2
##'-----------------------------------------------------------------------------------------#
des <- DESeqDataSetFromHTSeqCount(sampleTable,
                                  directory="./",
                                  design= ~ condition)
des <- DESeq(des, BPPARAM=BPPARAM)
res <- results(des, BPPARAM=BPPARAM)
##'-----------------------------------------------------------------------------------------#



##'Extract count data
##'-----------------------------------------------------------------------------------------#
raw_counts        <- counts(des, normalized=F)
normalised_counts <- counts(des, normalized=T)
mean_counts       <- data.frame(NOF_mean_norm= apply(normalised_counts[,grep("NOF",sampleTable$condition)], 1, mean),
                                OA_mean_norm = apply(normalised_counts[,grep("OA", sampleTable$condition)], 1, mean))
res_counts        <- cbind(res, mean_counts)
##'-----------------------------------------------------------------------------------------#



##'Get biomaRt Annotation
##'-----------------------------------------------------------------------------------------#
ensembl_genes_in <- rownames(res)
ensembl          <- useMart("ensembl",
                            dataset  = "hsapiens_gene_ensembl")
annotation       <- getBM(attributes = c("ensembl_gene_id",
                                         "external_gene_name",
                                         "description"),
                          filters    = "ensembl_gene_id",
                          values     = ensembl_genes_in,
                          ensembl)
##'-----------------------------------------------------------------------------------------#



##'Mangle results with annotation
##'-----------------------------------------------------------------------------------------#
merged_output           <- merge(res_counts, annotation, by="ensembl_gene_id")
merged_order            <- merged_output[with(merged_output, order(padj)), ]
final_NOF_OA_DE         <- cbind(merged_order[,c(1,10,11)], merged_order[,c(2:9)])
final_NOF_OA_DE$abs.fc  <- abs(final_NOF_OA_DE$log2FoldChange)
final_NOF_OA_DE         <- final_NOF_OA_DE[,c(1:5, 12, 6:11)]
head(final_NOF_OA_DE)
# write.csv(final_NOF_OA_DE,
#           file="DESeq2_NOF_OA_DE_List_Whole.csv")
##'-----------------------------------------------------------------------------------------#
