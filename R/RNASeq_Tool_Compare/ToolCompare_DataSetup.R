#-------------------------------------------------------------------------------------------#
#  Author      : Andrew J Skelton                                                           |
#  Language    : R Statistical Programming Language                                         |
#  Study       : RNA Seq Quantification Tool Compare                                        |
#  Data Owner  : Worked Example Using NOF OA RNA Seq Data - Prof. David Young               |
#  Description : RNA Seq Quantification Tool Compare - Set Up Dataset                       |
#-------------------------------------------------------------------------------------------#



##'Set the Working Directory, load essential packages, and import files
##'-----------------------------------------------------------------------------------------#
source("http://bioconductor.org/biocLite.R")
biocLite()

library(reshape2)
library(scales)
library(ggplot2)

setwd("~/Documents/Repos/Bioinformatics_Sprt_Code/R/RNASeq_Tool_Compare")

load("RData/Annotation_Data.Rdata")
load("RData/DESeq2_Data.Rdata")
load("RData/Salmon_Data.Rdata")
load("RData/Kallisto_Data.Rdata")
load("RData/Tuxedo_Data.Rdata")

tuxedo_mega_map <- merge(tuxedo_map, annotation, by.x="nearest_ref_id", by.y="ensembl_transcript_id")
##'-----------------------------------------------------------------------------------------#


##'
##'-----------------------------------------------------------------------------------------#
head(Tuxedo_Transcript_Counts)
head(Tuxedo_Gene_Counts)
head(DESeq2_Gene_Counts)
head(Kallisto_Transcript_Counts)
head(Salmon_Gene_Counts)
head(Salmon_Trans_Counts)

colnames(Salmon_Gene_Counts)       <- pheno_in$Sample
colnames(Salmon_Trans_Counts)      <- pheno_in$Sample
colnames(Tuxedo_Transcript_Counts) <- c(as.vector(pheno_in[pheno_in$Treatment == "NOF",]$Sample),
                                        as.vector(pheno_in[pheno_in$Treatment == "OA",]$Sample))
colnames(Tuxedo_Gene_Counts)       <- c(as.vector(pheno_in[pheno_in$Treatment == "NOF",]$Sample),
                                        as.vector(pheno_in[pheno_in$Treatment == "OA",]$Sample))


Tux_Trans_melt           <- melt(as.matrix(Tuxedo_Transcript_Counts))
Tux_Trans_melt           <- merge(Tux_Trans_melt, tuxedo_map[,c(1,2)], by.x="Var1", by.y="tracking_id")
Tux_Trans_melt           <- Tux_Trans_melt[,c(4,2,3)]
colnames(Tux_Trans_melt) <- c("BioID", "SampleID", "Count")
Tux_Trans_melt$Type      <- pheno_in[match(Tux_Trans_melt$SampleID, pheno_in$Sample),]$Treatment
Tux_Trans_melt$Tool      <- "Tuxedo"
Tux_Trans_melt$unique    <- paste0(Tux_Trans_melt$BioID, "_", Tux_Trans_melt$SampleID)

gene_map                 <- unique(tuxedo_mega_map[,c(3,4)])
Tux_Gene_melt            <- melt(as.matrix(Tuxedo_Gene_Counts))
Tux_Gene_melt            <- merge(Tux_Gene_melt, gene_map, by.x="Var1", by.y="gene_id")
Tux_Gene_melt            <- Tux_Gene_melt[,c(4,2,3)]
colnames(Tux_Gene_melt)  <- c("BioID", "SampleID", "Count")
Tux_Gene_melt$Type       <- pheno_in[match(Tux_Gene_melt$SampleID, pheno_in$Sample),]$Treatment
Tux_Gene_melt$Tool       <- "Tuxedo"
Tux_Gene_melt$unique     <- paste0(Tux_Gene_melt$BioID, "_", Tux_Gene_melt$SampleID)

DES_Gene_melt            <- melt(as.matrix(DESeq2_Gene_Counts))
colnames(DES_Gene_melt)  <- c("BioID", "SampleID", "Count")
DES_Gene_melt$Type       <- pheno_in[match(DES_Gene_melt$SampleID, pheno_in$Sample),]$Treatment
DES_Gene_melt$Tool       <- "DESeq2"
DES_Gene_melt$unique     <- paste0(DES_Gene_melt$BioID, "_", DES_Gene_melt$SampleID)

Sal_Gene_melt            <- melt(as.matrix(Salmon_Gene_Counts))
colnames(Sal_Gene_melt)  <- c("BioID", "SampleID", "Count")
Sal_Gene_melt$Type       <- pheno_in[match(Sal_Gene_melt$SampleID, pheno_in$Sample),]$Treatment
Sal_Gene_melt$Tool       <- "Salmon"
Sal_Gene_melt$unique     <- paste0(Sal_Gene_melt$BioID, "_", Sal_Gene_melt$SampleID)

Sal_Trans_melt           <- melt(as.matrix(Salmon_Trans_Counts))
colnames(Sal_Trans_melt) <- c("BioID", "SampleID", "Count")
Sal_Trans_melt$Type      <- pheno_in[match(Sal_Trans_melt$SampleID, pheno_in$Sample),]$Treatment
Sal_Trans_melt$Tool      <- "Salmon"
Sal_Trans_melt$unique    <- paste0(Sal_Trans_melt$BioID, "_", Sal_Trans_melt$SampleID)

Kal_Trans_melt           <- melt(as.matrix(Kallisto_Transcript_Counts))
colnames(Kal_Trans_melt) <- c("BioID", "SampleID", "Count")
Kal_Trans_melt$Type      <- pheno_in[match(Kal_Trans_melt$SampleID, pheno_in$Sample),]$Treatment
Kal_Trans_melt$Tool      <- "Kallisto"
Kal_Trans_melt$unique    <- paste0(Kal_Trans_melt$BioID, "_", Kal_Trans_melt$SampleID)

foo <- merge(Kal_Trans_melt, Sal_Trans_melt, by="unique")
foo <- foo[,c(1:4,9,10)]
colnames(foo) <- c("unique", "BioID", "SampleID", "KallistoCount", "SalmonCount", "Type")

# df_in <- foo[foo$SampleID == "Flowcell_A_1",]
df_in <- foo
df_in$SalmonCount   <- log2(df_in$SalmonCount + 1)
df_in$KallistoCount <- log2(df_in$KallistoCount + 1)



png("test.png", width=4096, height=8096, units="px", res=300)
ggplot(df_in, aes(x=KallistoCount, y=SalmonCount)) +
  geom_point() +
  theme_bw() + 
  facet_grid(SampleID ~ Type) 
dev.off()



