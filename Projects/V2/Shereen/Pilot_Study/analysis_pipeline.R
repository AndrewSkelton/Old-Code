#-------------------------------------------------------------------------------------------#
#  Author      : Andrew J Skelton                                                           |
#  Language    : R Statistical Programming Language                                         |
#  Study       : Primary Array - Pilot Study - pSS, Lymphoma                                |
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
setwd("~/Documents/Bioinformatics/Customers/Shereen/Shereen_Sept_2014/")
filename              <- "Fai_Sample_Probe_Profile.txt"
raw_data              <- lumiR(filename)
##'Gender Discrepancies
raw_data              <- raw_data[, -c(2, 6)]
pheno_table           <- read.table("sample_info.txt",
                                    header=T,
                                    sep="\t",
                                    stringsAsFactors=F)
rownames(pheno_table) <- pheno_table$Sentrix_ID
pData(raw_data)       <- pheno_table
##'-----------------------------------------------------------------------------------------#




##'Filter Input Data
##'-----------------------------------------------------------------------------------------#


##'CLINICAL OUTLIERS
##'"IPS-004-1", - Detected Outlier
removal     <- c("IPS-004-1", "IPS-002-1", "NCL-130-0", "SWI-084-0", "LEE-060-0", "NCL-113-0", "TOR-006-1")
removal_pos <- match(rownames(pData(raw_data)[pData(raw_data)$SampleID %in% removal,]), colnames(raw_data))
raw_data_in <- raw_data[,-removal_pos]


##'RIN SCORE EXCLUSION
##'Arrays Less than 7
Seven <- c("LEE-060-0","BIR-033-1","NCL-097-0","NCL-007-1","TOR-006-1","BIR-039-1","LEE-034-1",
           "LEE-062-1","TOR-007-0","NOT-036-1","FIF-027-1","BIR-011-1","SWI-006-1","BIR-030-1",
           "IPS-002-1","DER-019-1","FIF-014-1","WIN-016-1","LEE-016-1","LEE-012-1","GLA-019-1",
           "BIR-041-1","SWI-031-1","NCL-053-1")
removal_pos <- match(rownames(pData(raw_data_in)[pData(raw_data_in)$SampleID %in% Seven,]), colnames(raw_data_in))
raw_data_in <- raw_data_in[,-removal_pos]


##'Arrays Less than 7 with lymphoma
SevenL <- c("LEE-060-0", "BIR-033-1", "NCL-097-0", "NCL-007-1", "TOR-006-1", "BIR-039-1", "LEE-034-1",
            "LEE-062-1", "TOR-007-0", "BIR-011-1", "SWI-006-1", "BIR-030-1",
            "IPS-002-1", "DER-019-1", "FIF-014-1", "WIN-016-1", "LEE-016-1", "LEE-012-1", "GLA-019-1",
            "BIR-041-1", "SWI-031-1", "NCL-053-1")
removal_pos <- match(rownames(pData(raw_data_in)[pData(raw_data_in)$SampleID %in% SevenL,]), colnames(raw_data_in))
raw_data_in <- raw_data_in[,-removal_pos]



##'Arrays Less than 5
Five <- c("NCL-097-0", "BIR-011-1", "SWI-006-1", "BIR-041-1", "SWI-031-1", "NCL-053-1")
removal_pos <- match(rownames(pData(raw_data_in)[pData(raw_data_in)$SampleID %in% Five,]), colnames(raw_data_in))
raw_data_in <- raw_data_in[,-removal_pos]


##'BATCH EXCLUSION
raw_data_in <- raw_data_in[,-grep(1, pData(raw_data_in)$Batch)]
pData(raw_data_in)$Batch <- pData(raw_data_in)$Batch - 1
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




##'Normalise using RSN and VST - reccomended for Illumina Microarray Data
##'-----------------------------------------------------------------------------------------#
vst_data         <- lumiT(raw_data_in, method='vst')
rsn_data         <- lumiN(vst_data, method = "rsn")
lumi.Q           <- lumiQ(rsn_data)

##'Detection Threshold Filtering
exprs_data       <- exprs(lumi.Q)
present_count    <- detectionCall(lumi.Q)
normalised_data  <- exprs_data[present_count > 0, ]
##'-----------------------------------------------------------------------------------------#




##'Batch Correct for a technical effect identified in RNA Amplification
##'-----------------------------------------------------------------------------------------#
batches             <- pData(raw_data_in)$Batch
pheno               <- data.frame(sample=c(1:ncol(normalised_data)),
                                  outcome=as.factor(pData(lumi.Q)$Group),
                                  batch=batches)
rownames(pheno)     <- colnames(normalised_data)
batch               <- pheno$batch
mod                 <- model.matrix(~as.factor(outcome), data = pheno)
batchCorrected_data <- ComBat(dat=normalised_data,
                              batch=batch,
                              mod=mod,
                              par.prior=T,
                              prior.plots=F)
##'-----------------------------------------------------------------------------------------#




##'Probe Annotation Mapping
##'-----------------------------------------------------------------------------------------#
probe_list       <- rownames(batchCorrected_data)
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
treatments          <- unique(pData(raw_data_in)$Group)
treatment_arrays    <- pData(raw_data_in)$Group
batchCorrected_data <- data.matrix(batchCorrected_data)
design              <- model.matrix(~0 + factor(treatment_arrays, levels = treatments))
colnames(design)    <- treatments
num_parameters      <- ncol(design)
fit                 <- lmFit(batchCorrected_data, design)

cont_mat            <- makeContrasts(SS-Control,
                                     Lymphoma-Control,
                                     Cancer-Control,
                                     PreMalignancy-Control,
                                     SS-Lymphoma,
                                     SS-Cancer,
                                     SS-PreMalignancy,
                                     Lymphoma-Cancer,
                                     Lymphoma-PreMalignancy,
                                     Cancer-PreMalignancy,
                                     levels=treatments)
fit2                <- contrasts.fit(fit, contrasts=cont_mat)
fit2                <- eBayes(fit2)
##'-----------------------------------------------------------------------------------------#




##'Differential Expression Tests
##'-----------------------------------------------------------------------------------------#
comparisons <- c("SS - Control", "Lymphoma - Control", "Cancer - Control",
                 "PreMalignancy - Control", "SS - Lymphoma", "SS - Cancer",
                 "SS - PreMalignancy", "Lymphoma - Cancer", "Lymphoma - PreMalignancy",
                 "Cancer - PreMalignancy")
p_cut_off   <- 0.05
fold_change <- 1.2
i           <- 5

gene_list_unfiltered <- topTable(fit2,
                                 coef=comparisons[i],
                                 number=Inf,
                                 adjust.method="BH")

gene_list            <- topTable(fit2,
                                 coef=comparisons[i],
                                 p.value=p_cut_off,
                                 lfc=log2(fold_change),
                                 number=Inf,
                                 adjust.method="BH")
##'-----------------------------------------------------------------------------------------#
