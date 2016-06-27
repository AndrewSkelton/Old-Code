#-------------------------------------------------------------------------------------------#
#  Author      : Andrew J Skelton                                                           |
#  Language    : R Statistical Programming Language                                         |
#  Study       : miRNA Sequencing of Pooled Bovine Cartilage Samples                        |
#  Data Owner  : Newcastle University - Prof. David Young, (Cardiff Collaboration           |
#  Description : Experiment to detect the changes in miRNA expression under different       |
#                Pressures in Bovine Cartilage                                              |
#-------------------------------------------------------------------------------------------#



##'Set the Working Directory, load essential packages, and import files
##'-----------------------------------------------------------------------------------------#
setwd("~/Paulina/Resequencing/April/Alignment/mature/cow/")

source('http://bioconductor.org/biocLite.R')
biocLite()

library(stringr)
library(Rsamtools)
library(ShortRead)
library(DESeq2)
library(ggplot2)
##'-----------------------------------------------------------------------------------------#



##'Count miRNAs from alignment to Ensembl Annotation - Bos Taurus UMD 3.1 
##'-----------------------------------------------------------------------------------------#
hairpins        <- readDNAStringSet("../../../bowtie2Idx/DNA/mature/cow/cow_mature_dna.fa", 
                                    use.names=T)
names(hairpins) <- str_match(names(hairpins),
                             "^\\S+") 
hairpinGR       <- GRanges(names(hairpins), 
                           IRanges(1, 
                                   width(hairpins)), 
                           strand="*")
bamView         <- BamViews(list.files("./", 
                                       pattern="*.bam$"), 
                            bamRanges=hairpinGR)
bamcounts       <- countBam(bamView)
bamcounts_df    <- as.data.frame(bamcounts)
##'-----------------------------------------------------------------------------------------#



##'Build Count Table
##'-----------------------------------------------------------------------------------------#
out             <- data.frame(miRBase_ID=unique(bamcounts_df$space))
for(i in unique(bamcounts_df$group_name)) {
  tmp           <- bamcounts_df[bamcounts_df$group_name == i,]
  tmp           <- tmp[,c(3,8)]
  out           <- merge(out, 
                         tmp, 
                         by.x="miRBase_ID", 
                         by.y="space")
  colnames(out)[ncol(out)] <- i
}
rownames(out)   <- out$miRBase_ID
out             <- out[,-1]
##'-----------------------------------------------------------------------------------------#



##'DESeq2 - Pheno Prep
##'-----------------------------------------------------------------------------------------#
setwd("~/Paulina/Resequencing/April/DESeq2/")
sampleTable             <- read.table("pheno.txt", 
                                      header=T, 
                                      row.names=1, 
                                      sep="\t", 
                                      stringsAsFactors=T)
colnames(out)           <- rownames(sampleTable)
sampleTable_bck         <- sampleTable
sampleTable$UniqueClass <- paste0(sampleTable$Pressure, 
                                  sampleTable$TimePoint)
sampleTable$TimePoint   <- factor(sampleTable$TimePoint, 
                                  levels=c("2hr", "6hr", "24hr"))
sampleTable$Prefix      <- factor(sampleTable$Prefix, 
                                  labels=c("A","B","C","D","E","F","G","H","I"))
##'-----------------------------------------------------------------------------------------#



##'DESeq2 - Filter Pheno and Raw Count Inputs for Given parameter
##'-----------------------------------------------------------------------------------------#
out_in          <- out[, grep("24hr", sampleTable$TimePoint)]
sampleTable_in  <- sampleTable[sampleTable$TimePoint == "24hr",]
# out_in          <- out[, grep("7MPa", sampleTable$Pressure)]
# sampleTable_in  <- sampleTable[sampleTable$Pressure == "7MPa",]
##'-----------------------------------------------------------------------------------------#



##'DESeq2 - Run
##'-----------------------------------------------------------------------------------------#
dds2           <- DESeqDataSetFromMatrix(countData = out_in,
                                         colData   = sampleTable_in,
                                         design    = ~ Prefix + Pressure)
dds2           <- DESeq(dds2)

# dds2           <- DESeqDataSetFromMatrix(countData = out_in,
#                                          colData   = sampleTable,
#                                          design    = ~ TimePoint)
# dds2           <- DESeq(dds2)
##'-----------------------------------------------------------------------------------------#



##'DESeq2 - Useful
##'-----------------------------------------------------------------------------------------#
deseq_rld      <- rlog(dds2)
deseq_vst      <- varianceStabilizingTransformation(dds2)
deseq_rlogMat  <- assay(deseq_rld)
deseq_vstMat   <- assay(deseq_vst)

ggpca(assay(deseq_vst), sampleTable_in, 1, 8, F, T)
ggpca(assay(deseq_vst), sampleTable_in, 1, 7, T, T)
ggpca(assay(deseq_vst), sampleTable_in, 1, 6, F, T)
ggpca(assay(deseq_vst), sampleTable_in, 1, 5, T, T)

resultsNames(dds2)

deseq_normCoun <- counts(dds2, normalized=T)
deseq_rawCoun  <- counts(dds2, normalized=F)
##'-----------------------------------------------------------------------------------------#



##'DESeq2 - Pressure Contrasts
##'-----------------------------------------------------------------------------------------#
res_7_0  <- results(dds2)
res_7_0  <- res_7_0[order(res_7_0$padj),]
res_7_0  <- as.data.frame(res_7_0)

res_7_25 <- results(dds2, contrast=c("Pressure","7MPa","2.5MPa"))
res_7_25 <- res_7_25[order(res_7_25$padj),]
res_7_25 <- as.data.frame(res_7_25)

res_25_0 <- results(dds2, contrast=c("Pressure","2.5MPa","0MPa"))
res_25_0 <- res_25_0[order(res_25_0$padj),]
res_25_0 <- as.data.frame(res_25_0)

# P24hr_out <- list(res_7_0, res_7_25, res_25_0, deseq_normCoun, deseq_rawCoun)
# save(P24hr_out, file="P24hr.Rdata")
# res2[order(res2$padj),]
##'-----------------------------------------------------------------------------------------#



##'DESeq2 - Time Contrasts
##'-----------------------------------------------------------------------------------------#
# deseq_normCoun <- counts(dds2, normalized=T)
# deseq_rawCoun  <- counts(dds2, normalized=F)
# res_24_2  <- results(dds2)
# res_24_2  <- res_24_2[order(res_24_2$padj),]
# res_24_2  <- as.data.frame(res_24_2)
# 
# res_24_6  <- results(dds2, contrast=c("TimePoint","24hr","6hr"))
# res_24_6  <- res_24_6[order(res_24_6$padj),]
# res_24_6  <- as.data.frame(res_24_6)
# 
# res_6_2   <- results(dds2, contrast=c("TimePoint","6hr","2hr"))
# res_6_2   <- res_6_2[order(res_6_2$padj),]
# res_6_2   <- as.data.frame(res_6_2)
# 
# T7mpa_out <- list(res_24_2, res_24_6, res_6_2, deseq_normCoun, deseq_rawCoun)
# save(T7mpa_out, file="T7mpa.Rdata")
##'-----------------------------------------------------------------------------------------#


##'Visualise miRNA
##'-----------------------------------------------------------------------------------------#
plot_paired("bta-mir-148b", 
            deseq_normCoun,
            sampleTable_in,
            "24hr - 7Mpa Vs 2.5Mpa")
##'-----------------------------------------------------------------------------------------#

