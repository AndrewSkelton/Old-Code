#-------------------------------------------------------------------------------------------#
#  Author      : Andrew J Skelton                                                           |
#  Language    : R Statistical Programming Language                                         |
#  Study       : Data Combination Project - Early Onset RA Study                            |
#  Data Owner  : Newcastle University - Dr. Arthur Pratt                                    |
#  Description : Bespoke Pipeline to Combine Pilot Study Data (Illumina WG-v3), and         |
#                Confirmation Study Data (Illumina HT-12v4). Accounting for Technical       |
#                Effects in all Arrays.                                                     |
#-------------------------------------------------------------------------------------------#



##'Set the Working Directory, load essential packages, and import files
##'-----------------------------------------------------------------------------------------#
source("http://bioconductor.org/biocLite.R")
biocLite()

library(inSilicoDb)
library(inSilicoMerging)
library(lumi)
library(lumidat)
library(sva)
library(ggplot2)
library(annotate)
library(lumiHumanAll.db)
library(limma)
library(reshape2)
library(scales)
library(pathview)

setwd("~/Documents/Bioinformatics/Customers/Pratt/Project_Master_Branch/")

pheno_selection <- read.table("Data_Combination/Code_Base/pheno_selection.txt",
                              header=T, sep="\t", stringsAsFactors=F)

signature       <- c("ro9U1LFU0.WLkLnbnI","ERRUHoLT147nuJIOpc",
                     "ce1oOUVE31P1O1XVIk","Tyao14HoeDAG63VHmU",
                     "cU6IqFSTKAlEx0l6Oo","0np1587tLe3XXfu1Co",
                     "lpcsVWbpc1HXRNREeQ","oevLuAy95Rd4ZFHlDk",
                     "BNbgFq61QicgHSgFUo","fr1fXPUJ9TJ55XL.f0",
                     "TUiTQIR9V11IpMeQRA","lqWI6.t85cye.O1RS4")
##'-----------------------------------------------------------------------------------------#


##'Load in Raw Data from Confirmation Study (Illumina HT-12v4)
##'Normalise using RSN VST Methods
##'Technical Correction for RNA Amplification - ComBat Method (Null Covariates, Batches Only)
##'-----------------------------------------------------------------------------------------#
raw_data_a            <- lumiR("Confirmation_Study/Raw_Data/probe_profile/bg_corrected Sample Probe Profile.txt")
pheno_table           <- read.table("Pheno_Data/Pheno_Confirmation_Study.txt",
                                    header=T, row.names=1, sep="\t", stringsAsFactors=F)
pheno_table$Source    <- "A"
pData(raw_data_a)     <- pheno_table
raw_data_a            <- raw_data_a[, -c(115,116,146,148,150,151,152,153,
                                         154,155,156,195,199,92,134)]

##'Optional: Remove UA Patients
raw_data_a <- raw_data_a[,grep("UA", pData(raw_data_a)$Baseline, invert=T)]

vst_data              <- lumiT(raw_data_a, method='vst')
rsn_data              <- lumiN(vst_data, method = "rsn")
data_alumi.Q          <- lumiQ(rsn_data)

normalised_data_a     <- exprs(data_alumi.Q)
pheno                 <- data.frame(sample=c(1:ncol(normalised_data_a)),
                                    outcome=pData(data_alumi.Q)$Baseline,
                                    batch=pData(data_alumi.Q)$Amplification_Batch)

rownames(pheno)       <- colnames(normalised_data_a)
batch                 <- pheno$batch
mod                   <- model.matrix(~as.factor(outcome), data=pheno)
batchCorrected_data_a <- ComBat(dat=normalised_data_a,
                                batch=batch,
                                mod=NULL,
                                #mod=mod,
                                par.prior=T,
                                prior.plots=F)
exprs(data_alumi.Q)   <- batchCorrected_data_a
##'-----------------------------------------------------------------------------------------#


##'Load in Raw Data from Pilot Study (Illumina HT-12v4) - Experiment 1
##'Normalise using RSN VST Methods
##'Technical Correction for RNA Amplification - ComBat Method (Null Covariates, Batches Only)
##'-----------------------------------------------------------------------------------------#
raw_data_b            <- lumiR("Pilot_Study/Raw_Data/Expt_1_raw_data.txt")
raw_data_b            <- raw_data_b[,-87]
raw_data_b            <- raw_data_b[-321,]
probe_list            <- rownames(raw_data_b)
raw_data_b            <- raw_data_b[-grep("iiipihinwAfdI0NuTU",
                                          as.vector(probeID2nuID(probe_list,
                                                                 chipVersion="HumanWG6_V3_0_R3_11282955_A",
                                                                 species="Human")[,"nuID"])),]
probe_list            <- rownames(raw_data_b)
rownames(raw_data_b)  <- as.vector(probeID2nuID(probe_list,
                                                chipVersion="HumanWG6_V3_0_R3_11282955_A",
                                                species="Human")[,"nuID"])
pheno_table           <- read.table("Pheno_Data/Pheno_PilotStudy_Ex1.txt",
                                    header=T, row.names=1, sep="\t", stringsAsFactors=F)
pheno_table$Source    <- "B"
pData(raw_data_b)     <- pheno_table

##'Optional: Remove UA Patients
raw_data_b <- raw_data_b[,grep("UA", pData(raw_data_b)$Baseline, invert=T)]

##'Optional: Remove Include Only Predefined Patients
# raw_data_b <- raw_data_b[,pData(raw_data_b)$EA_Number %in% pheno_selection$EA.No]

vst_data              <- lumiT(raw_data_b, method='vst')
rsn_data              <- lumiN(vst_data, method = "rsn")
data_blumi.Q          <- lumiQ(rsn_data)
normalised_data_b     <- exprs(data_blumi.Q)

pheno                 <- data.frame(sample=c(1:ncol(normalised_data_b)),
                                    outcome=pData(data_blumi.Q)$Baseline,
                                    batch=pData(data_blumi.Q)$Amplification_Batch)

rownames(pheno)       <- colnames(normalised_data_b)
batch                 <- pheno$batch
mod                   <- model.matrix(~as.factor(outcome), data=pheno)
batchCorrected_data_b <- ComBat(dat=normalised_data_b,
                                batch=batch,
                                mod=NULL,
                                #mod=mod,
                                par.prior=T,
                                prior.plots=F)
exprs(data_blumi.Q)   <- batchCorrected_data_b
##'-----------------------------------------------------------------------------------------#


##'Load in Raw Data from Pilot Study (Illumina HT-12v4) - Experiment 2
##'Normalise using RSN VST Methods
##'Technical Correction for RNA Amplification - ComBat Method (Null Covariates, Batches Only)
##'-----------------------------------------------------------------------------------------#
raw_data_c            <- lumiR("Pilot_Study/Raw_Data/Expt_2_raw_data.txt")
raw_data_c            <- raw_data_c[-321,]
probe_list            <- rownames(raw_data_c)
raw_data_c            <- raw_data_c[-grep("iiipihinwAfdI0NuTU",
                                          as.vector(probeID2nuID(probe_list,
                                                                 chipVersion="HumanWG6_V3_0_R3_11282955_A",
                                                                 species="Human")[,"nuID"])),]
probe_list            <- rownames(raw_data_c)
rownames(raw_data_c)  <- as.vector(probeID2nuID(probe_list,
                                                chipVersion="HumanWG6_V3_0_R3_11282955_A",
                                                species="Human")[, "nuID"])

pheno_table           <- read.table("Pheno_Data/Pheno_PilotStudy_Ex2.txt",
                                    header=T, row.names=1, sep="\t", stringsAsFactors=F)
pheno_table$Source    <- "C"
pData(raw_data_c)     <- pheno_table

##'Optional: Remove UA Patients
raw_data_c <- raw_data_c[,grep("UA", pData(raw_data_c)$Baseline, invert=T)]

##'Optional: Remove Include Only Predefined Patients
# raw_data_c <- raw_data_c[,pData(raw_data_c)$EA_Number %in% pheno_selection$EA.No]

vst_data              <- lumiT(raw_data_c, method='vst')
rsn_data              <- lumiN(vst_data, method = "rsn")
data_clumi.Q          <- lumiQ(rsn_data)
normalised_data_c     <- exprs(data_clumi.Q)

pheno                 <- data.frame(sample=c(1:ncol(normalised_data_c)),
                                    outcome=pData(data_clumi.Q)$Baseline,
                                    batch=pData(data_clumi.Q)$Amplification_Batch)

rownames(pheno)       <- colnames(normalised_data_c)
batch                 <- pheno$batch
mod                   <- model.matrix(~as.factor(outcome), data=pheno)
batchCorrected_data_c <- ComBat(dat=normalised_data_c,
                                batch=batch,
                                mod=NULL,
                                #mod=mod,
                                par.prior=T,
                                prior.plots=F)
exprs(data_clumi.Q)   <- batchCorrected_data_c
##'-----------------------------------------------------------------------------------------#


##'Combine using the inSilicoMerging Package
##'-----------------------------------------------------------------------------------------#
data_sets         <- list(data_alumi.Q, data_blumi.Q, data_clumi.Q)
combined_combat   <- inSilicoMerging::merge(data_sets, method="COMBAT")
combined_none     <- inSilicoMerging::merge(data_sets, method="NONE")
raw_data_det      <- combined_combat
##'-----------------------------------------------------------------------------------------#


##'Build Annotation
##'-----------------------------------------------------------------------------------------#
nuIDs             <- rownames(exprs(raw_data_det))
symbol            <- getSYMBOL(nuIDs, "lumiHumanAll.db")
name              <- unlist(lookUp(nuIDs, "lumiHumanAll.db", "GENENAME"))
anno_df           <- data.frame(ID = nuIDs, symbol, name)
##'-----------------------------------------------------------------------------------------#


##'Differential Expression Model - Limma
##'-----------------------------------------------------------------------------------------#
treatments        <- unique(pData(raw_data_det)$Update_Collapse_2)
treatment_arrays  <- pData(raw_data_det)$Update_Collapse_2
design            <- model.matrix(~0 + factor(treatment_arrays, levels = treatments))
colnames(design)  <- treatments
num_parameters    <- ncol(design)
fit               <- lmFit(exprs(raw_data_det), design)
cont_mat          <- makeContrasts(NRA-RA, levels=design)
fit2              <- contrasts.fit(fit, contrasts=cont_mat)
fit2              <- eBayes(fit2)
fit2$genes        <- anno_df
##'-----------------------------------------------------------------------------------------#


##'Test for Differential Expression
##'-----------------------------------------------------------------------------------------#
comparisons       <- c("NRA - RA")
p_cut_off         <- 0.05
fold_change       <- 1.2
mtc               <- 'BH'
gene_list         <- topTable(fit2, coef=1, p.value=p_cut_off, lfc=log2(fold_change),
                              number=Inf, adjust.method=mtc)
filtered_in       <- gene_list
##'-----------------------------------------------------------------------------------------#
