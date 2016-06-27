#-------------------------------------------------------------------------------------------#
#  Author      : Andrew J Skelton                                                           |
#  Language    : R Statistical Programming Language                                         |
#  Study       : NOF and OA Microarray Study                                                |
#  Data Owner  : Newcastle University - Prof. David Young                                   |
#  Description : Illumina HT-12v3 Microarray of NOF and OA Samples. Grouped analysis and    |
#                Regression by Age.                                                         |
#-------------------------------------------------------------------------------------------#



##'Set the Working Directory, load essential packages, and import files
##'-----------------------------------------------------------------------------------------#
setwd("~/Young/Pipeline_Analysis/dexseq/grch38/")


source("http://bioconductor.org/biocLite.R")
biocLite()

library(DEXSeq)
library(ggplot2)
library(BiocParallel)
library(biomaRt)

BPPARAM     <- MulticoreParam(workers=20)
##'-----------------------------------------------------------------------------------------#



##'Build Pheno Table
##'-----------------------------------------------------------------------------------------#
sampleTable <- data.frame(
  row.names = c("F1_S1","F1_S2","F1_S3","F1_S4","F1_S5","F1_S6","F1_S7","F1_S8",
                "F2_S1","F2_S2","F2_S3","F2_S4","F2_S5","F2_S6","F2_S7","F2_S8"),
  condition = c("OA","NOF","OA","OA","NOF","OA","NOF","OA","OA","NOF","OA","OA",
                "NOF","NOF","OA","OA"),
  libType   = c(rep("paired-end", 16)),
  countName = list.files("./", full.names=T, pattern="*.count$")
)
##'-----------------------------------------------------------------------------------------#



##'Build a DEXSeq Dataset and run DEXSeq
##'-----------------------------------------------------------------------------------------#
dxd <- DEXSeqDataSetFromHTSeq(list.files("./", full.names=T, pattern="*.count"),
                              sampleData    = sampleTable,
                              design        = ~ sample + exon + condition:exon,
                              flattenedfile = "/data/genomes/GRCh38_Ensembl_Homo_Sapiens/hs_grch38_dexseq_ref/hs_grch38_dexseq.gff")
dxd  <- estimateSizeFactors(dxd)
dxd  <- estimateDispersions(dxd, BPPARAM=BPPARAM)
dxd  <- testForDEU(dxd, BPPARAM=BPPARAM)
dxd  <- estimateExonFoldChanges(dxd, fitExpToVar="condition", BPPARAM=BPPARAM)
dxr  <- DEXSeqResults(dxd)
##'-----------------------------------------------------------------------------------------#



##'Extract count data and inspect a gene
##'-----------------------------------------------------------------------------------------#
normalised_counts <- counts(dxr,normalized=T)
gene_of_interest  <- normalised_counts[grep("ENSG00000113580", rownames(normalised_counts)),]
gene_of_interest  <- cbind(rowMeans(gene_of_interest[,grep("OA", sampleTable$condition)]),
                           rowMeans(gene_of_interest[,grep("NOF", sampleTable$condition)]))
colnames(gene_of_interest) <- c("NOF", "OA")
# write.csv(gene_of_interest, file="NR3C1_Exon_Counts.csv")
##'-----------------------------------------------------------------------------------------#



##'Build a html table of all differentionally expressed exons
##'-----------------------------------------------------------------------------------------#
DEXSeqHTML(dxr, FDR=0.1, fitExpToVar="condition", BPPARAM=BPPARAM)

# png("NR3C1_Exon_Usage_and_Expression.png", width=10.5, height=12.53, units="in", res=600)
plotDEXSeq(dxr,
           "ENSG00000113580",
           expression=T,
           splicing=T,
           norCounts=T,
           legend=T,
           displayTranscripts=T,
           names=T,
           cex.axis=1,
           cex=1,
           lwd=1.2)
# dev.off()
##'-----------------------------------------------------------------------------------------#



##'Get biomaRt annotation
##'-----------------------------------------------------------------------------------------#
#GRCh37
ensembl <- useMart(biomart="ENSEMBL_MART_ENSEMBL",
                   dataset="hsapiens_gene_ensembl",
                   host="grch37.ensembl.org",
                   path="/biomart/martservice")

anno_df <- getBM(attributes=c("ensembl_gene_id",
                              "external_gene_name",
                              "description"),
                 filters="ensembl_gene_id",
                 values=substr(foo$groupID,1,15),
                 ensembl)

#GRCh38
ensembl    <- useMart("ensembl",
                      dataset="hsapiens_gene_ensembl")
annotation <- getBM(attributes=c("ensembl_gene_id",
                                 "external_gene_name",
                                 "description"),
                    filters="ensembl_gene_id",
                    values=ensembl_genes_in,
                    ensembl)
##'-----------------------------------------------------------------------------------------#
