#-------------------------------------------------------------------------------------------#
#  Author      : Andrew J Skelton                                                           |
#  Language    : R Statistical Programming Language                                         |
#  Study       : Globin Clear vs Full RNA                                                   |
#  Data Owner  : Newcastle University - Prof. Fai Ng, Shereen Al-Ali                        |
#  Description : Primary Analysis Pipeline for Illumina HT-12 v4 Microarray Data - Bespoke  |
#                for use in analysis the described study                                    |
#-------------------------------------------------------------------------------------------#




##'Set the Working Directory and load essential packages
##'-----------------------------------------------------------------------------------------#
setwd("/Users/andrew/Documents/Bioinformatics/Customers/Shereen/Globin_Clear")

source("http://bioconductor.org/biocLite.R")
biocLite(c("lumi", "gplots", "ggplot2", "limma", "annotate", "lumiHumanAll.db", "limma", "sva",
           "lumiHumanIDMapping"))
install.packages(c("scales", "reshape2"))
library(stringr)
library(sva)
library(lumi)
library(gplots)
library(ggplot2)
library(annotate)
library(lumiHumanAll.db)
library(limma)
library(scales)
library(reshape2)
library(lumiHumanIDMapping)
##'-----------------------------------------------------------------------------------------#




##'Load in Raw Data and attach pheno data
##'-----------------------------------------------------------------------------------------#
filename              <- "raw_data_SA.txt"
raw_data              <- lumiR(filename)
pheno_table           <- read.table("pheno.txt",
                                    header=T,
                                    sep="\t",
                                    row.names=1,
                                    stringsAsFactors=F)
pData(raw_data)       <- pheno_table
##'-----------------------------------------------------------------------------------------#




##'Explore the Detection P Values from the Scanner to look for technical failures
##'-----------------------------------------------------------------------------------------#
det <- melt(detection(raw_data))
dtp <- ggplot(data=det, aes(x=Var2, y=value)) +
              geom_boxplot(outlier.size=0.5,
                           size=0.2) +
              scale_x_discrete(name="") +
              scale_y_continuous(name="Amplitude") +
              theme_bw() +
              theme(axis.text.x=element_text(angle=90,
                                             vjust=0.5,
                                             size=6),
                    axis.text.y=element_text(size=6),
                    axis.title.y=element_text(size=6))

image_deploy(dtp, "DetectionPval_")
##'-----------------------------------------------------------------------------------------#




##'Choose what samples to include for normalisation. Seperated to see the technical
##'differences between probes passing detection threshold, and the fact that limma MAY use
##'all samples as an estimate for dispersion.
##'-----------------------------------------------------------------------------------------#
raw_data_det          <- raw_data
raw_data_det          <- raw_data[, grep("Globin_Clear", pData(raw_data)$class)]
# raw_data_det          <- raw_data[, grep("Full_RNA",     pData(raw_data)$class)]
##'-----------------------------------------------------------------------------------------#




##'Normalise using RSN and VST - reccomended for Illumina Microarray Data
##'-----------------------------------------------------------------------------------------#
vst_data         <- lumiT(raw_data_det, method='vst')
rsn_data         <- lumiN(vst_data, method = "rsn")
lumi.Q           <- lumiQ(rsn_data)

##'Detection Threshold Filtering
exprs_data       <- exprs(lumi.Q)
present_count    <- detectionCall(lumi.Q)
normalised_data  <- exprs_data[present_count > 0, ]
##'-----------------------------------------------------------------------------------------#




##'Probe Annotation Mapping
##'-----------------------------------------------------------------------------------------#
probe_list       <- rownames(normalised_data)
nuIDs            <- probeID2nuID(probe_list)[, "nuID"]
symbol           <- getSYMBOL(nuIDs, "lumiHumanAll.db")
name             <- unlist(lookUp(nuIDs, "lumiHumanAll.db", "GENENAME"))
anno_df          <- data.frame(ID = nuIDs, probe_list, symbol, name)
entrez_map       <- data.frame(nuID=as.vector(anno_df$ID),
                               EntrezID=nuID2EntrezID(as.vector(anno_df$ID),
                                                      "lumiHumanIDMapping"))
##'-----------------------------------------------------------------------------------------#




##'Limma Factorial Group Experimental Design
##'-----------------------------------------------------------------------------------------#
design              <- model.matrix(~0 + factor(pData(raw_data_det)$treatment,
                                                levels=c("SS", "control")))
colnames(design)    <- c("SS", "control")
num_parameters      <- ncol(design)
fit                 <- lmFit(normalised_data, design)
cont_mat            <- makeContrasts(SS-control, levels=c("SS", "control"))
fit2                <- contrasts.fit(fit, contrasts=cont_mat)
fit2                <- eBayes(fit2)
fit2$genes          <- anno_df
##'-----------------------------------------------------------------------------------------#




##'Differential Expression Tests
##'-----------------------------------------------------------------------------------------#
comparisons <- c("SS - control")
p_cut_off   <- 0.05
fold_change <- 1.2
i           <- 1

gene_list_unfiltered <- topTable(fit2,
                                 coef="SS - control",
                                 number=Inf,
                                 adjust.method="BH")

gene_list            <- topTable(fit2,
                                 coef="SS - control",
                                 p.value=p_cut_off,
                                 lfc=log2(fold_change),
                                 number=Inf,
                                 adjust.method="BH")
##'-----------------------------------------------------------------------------------------#
