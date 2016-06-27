#-------------------------------------------------------------------------------------------#
#  Author      : Andrew J Skelton                                                           |
#  Language    : R Statistical Programming Language                                         |
#  Study       : Secondary Array - Confirmation Study                                                   |
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

library(devtools)
install_github("drmjc/lumidat")
library(lumidat)
##'-----------------------------------------------------------------------------------------#




##'Load in Raw Data and attach pheno data
##'-----------------------------------------------------------------------------------------#
setwd("/Users/andrew/Documents/Bioinformatics/Customers/Shereen/Sj_Ly_Jan15")

raw.lumidat     <- lumiR.idat(files=list.files("raw_data/"),
                              path=file.path(getwd(), "raw_data/"),
                              probeID="ProbeID",
                              manifestfile="HumanHT-12_V4_0_R2_15002873_B.txt",
                              controls=F, detectionTh = 0.01,
                              backgroundCorrect=T, collapseMode="none",
                              QC=T, memory="-Xmx8080m", verbose=T)
pheno_table     <- read.table("pheno.txt", header=T, row.name=1, sep="\t", stringsAsFactors=F)
##'-----------------------------------------------------------------------------------------#




##'Explore the Detection P Values from the Scanner to look for technical failures
##'-----------------------------------------------------------------------------------------#
det <- melt(detection(raw.lumidat))
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




##'Choose what samples to include for normalisation. QC exclusions.
##'-----------------------------------------------------------------------------------------#
raw.lumidat.filtered <- raw.lumidat[,-(grep("9989522041_F$|9989522055_I$|9989522055_K$|9989522055_D$|9989522051_B$|9989522051_H$|9989522041_D$",
                                            colnames(raw.lumidat)))]
pheno_table.filtered <- pheno_table[-(grep("9989522041_F$|9989522055_I$|9989522055_K$|9989522055_D$|9989522051_B$|9989522051_H$|9989522041_D$",
                                           colnames(raw.lumidat))),]
pData(raw.lumidat.filtered) <- pheno_table.filtered
##'-----------------------------------------------------------------------------------------#




##'Normalise using RSN and VST - reccomended for Illumina Microarray Data
##'-----------------------------------------------------------------------------------------#
vst_data         <- lumiT(raw.lumidat.filtered, method='vst')
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
treatments          <- unique(pData(lumi.Q)$Treatment)
treatment_arrays    <- pData(lumi.Q)$Treatment
batchCorrected_data <- data.matrix(normalised_data)
design              <- model.matrix(~0 + factor(treatment_arrays,
                                                levels=treatments))
colnames(design)    <- treatments
num_parameters      <- ncol(design)
fit                 <- lmFit(normalised_data, design)

cont_mat            <- makeContrasts(SS-Control,
                                     lymphoma-Control,
                                     SS-lymphoma,
                                     levels=treatments)
fit2                <- contrasts.fit(fit, contrasts=cont_mat)
fit2                <- eBayes(fit2)
fit2$genes          <- anno_df
##'-----------------------------------------------------------------------------------------#




##'Differential Expression Tests
##'-----------------------------------------------------------------------------------------#
comparisons <- c("SS - Control", "lymphoma - Control", "SS - lymphoma")
p_cut_off   <- 0.05
fold_change <- 1.2
mtc         <- 'BH'
contrast    <- comparisons[3]

gene_list_unfiltered <- topTable(fit2,
                                 coef=contrast,
                                 number=Inf,
                                 adjust.method="BH")

gene_list            <- topTable(fit2,
                                 coef=contrast,
                                 p.value=p_cut_off,
                                 lfc=log2(fold_change),
                                 number=Inf,
                                 adjust.method="BH")
##'-----------------------------------------------------------------------------------------#
