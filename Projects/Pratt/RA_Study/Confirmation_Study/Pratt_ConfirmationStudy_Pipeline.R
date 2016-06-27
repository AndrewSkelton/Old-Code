#-------------------------------------------------------------------------------------------#
#  Author      : Andrew J Skelton                                                           |
#  Language    : R Statistical Programming Language                                         |
#  Study       : Early Onset RA Study - Confirmation Study (CD4 T-Cells)                    |
#  Data Owner  : Newcastle University - Dr. Arthur Pratt                                    |
#  Description : Bespoke Pipeline to Analyse Illumina HT-12v4 Arrays                        |
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
library(gplots)

setwd("/Users/andrew/Documents/Bioinformatics/Customers/Pratt/Project_Master_Branch/Confirmation_Study")

pheno_table           <- read.table("../Pheno_Data/Ind_Pheno_Confirmation_Study.txt", 
                                    header=T, sep="\t", stringsAsFactors=F)
rownames(pheno_table) <- pheno_table$X
##'-----------------------------------------------------------------------------------------#



##'Read in Raw Data and Pheno Data
##'-----------------------------------------------------------------------------------------#
raw_data        <- lumiR("Raw_Data/probe_profile/bg_corrected Sample Probe Profile.txt", 
                         verbose=T)
# raw_data        <- addControlData2lumi("control_probe_profile.txt", raw_data)
# raw_data        <- lumiB(raw_data, method="bgAdjust", verbose=T)
pData(raw_data) <- pheno_table
##'-----------------------------------------------------------------------------------------#



##'Outlier Removal
##'-----------------------------------------------------------------------------------------#
##' 115/116 (879/948) - Technical Failure
##' 146/134 (835rpt/835) - Outliers (possible sample contamination), vst/rsn Distribution off
##'
##' Repeat Samples
##' 148/150/151/152/153/154/155/156
##' (836rpt/837rpt/832rpt/838rpt/833rpt/839rpt/834rpt/840rpt)
##' 
##' PCA Outliers - 923/926 (195/199)
##' Hierarchical Clustering Outliers - 9477874168_H - 92 - 2233
raw_data_det <- raw_data[, -c(115,116, 146,148,150,151,152,153,154,155,156, 195,199, 92, 134)]
##'-----------------------------------------------------------------------------------------#


##'Optional - UA Exclusion
##'-----------------------------------------------------------------------------------------#
raw_data_det <- raw_data_det[, pData(raw_data_det)$Baseline != "UA"]
##'-----------------------------------------------------------------------------------------#


##'Outlier Removal
##'-----------------------------------------------------------------------------------------#
vst_data         <- lumiT(raw_data_det, method='vst')
rsn_data         <- lumiN(vst_data,     method='rsn')
lumi.Q           <- lumiQ(rsn_data)
##'-----------------------------------------------------------------------------------------#



##'Detection P Value Filtering
##'-----------------------------------------------------------------------------------------#
exprs_data      <- exprs(lumi.Q)
present_count   <- detectionCall(lumi.Q)
normalised_data <- exprs_data[present_count > 0, ]
##'-----------------------------------------------------------------------------------------#



##'Technical Correction - ComBat
##'-----------------------------------------------------------------------------------------#
pheno               <- data.frame(sample=c(1:ncol(normalised_data)),
                                  outcome=pData(lumi.Q)$Baseline, 
                                  batch=pData(lumi.Q)$RNA.Batch)
rownames(pheno)     <- colnames(normalised_data)
batch               <- pheno$batch
mod                 <- model.matrix(~as.factor(outcome), 
                                    data = pheno)
batchCorrected_data <- ComBat(dat=normalised_data, 
                              batch=batch, 
                              mod=NULL,
                              # mod=mod, 
                              par.prior=T, 
                              prior.plots=F)
##'-----------------------------------------------------------------------------------------#



##'Probe Annotation
##'-----------------------------------------------------------------------------------------#
nuIDs            <- rownames(batchCorrected_data)
symbol           <- getSYMBOL(nuIDs, "lumiHumanAll.db")
name             <- unlist(lookUp(nuIDs, "lumiHumanAll.db", "GENENAME"))
anno_df          <- data.frame(ID=nuIDs, probe_list=nuIDs, symbol, name)
##'-----------------------------------------------------------------------------------------#



##'Experimental Design - Limma Model Fit
##'-----------------------------------------------------------------------------------------#
treatments          <- unique(pData(lumi.Q)$Baseline)
treatment_arrays    <- pData(lumi.Q)$Baseline
batchCorrected_data <- data.matrix(batchCorrected_data)
design              <- model.matrix(~0 + factor(treatment_arrays, 
                                                levels=treatments))
colnames(design)    <- treatments
num_parameters      <- ncol(design)
fit                 <- lmFit(batchCorrected_data, design)
cont_mat            <- makeContrasts(NRA-RA, levels=treatments)
fit2                <- contrasts.fit(fit, contrasts=cont_mat)
fit2                <- eBayes(fit2)
fit2$genes          <- anno_df
##'-----------------------------------------------------------------------------------------#



##'Differential Expression Test
##'-----------------------------------------------------------------------------------------#
comparisons <- c("NRA - RA")
p_cut_off   <- 0.05
fold_change <- 1.2
mtc         <- 'none'

gene_list   <- topTable(fit2, 
                        coef=comparisons[1], 
                        p.value=p_cut_off, 
                        lfc=log2(fold_change),
                        number=Inf, 
                        adjust.method=mtc,
                        sort.by="P")
filtered_in <- gene_list
##'-----------------------------------------------------------------------------------------#



##'12 Genes
##'-----------------------------------------------------------------------------------------#
signature       <- c("ro9U1LFU0.WLkLnbnI","ERRUHoLT147nuJIOpc",
                     "ce1oOUVE31P1O1XVIk","Tyao14HoeDAG63VHmU",
                     "cU6IqFSTKAlEx0l6Oo","0np1587tLe3XXfu1Co",
                     "lpcsVWbpc1HXRNREeQ","oevLuAy95Rd4ZFHlDk",
                     "BNbgFq61QicgHSgFUo","fr1fXPUJ9TJ55XL.f0",
                     "TUiTQIR9V11IpMeQRA","lqWI6.t85cye.O1RS4")



##'-----------------------------------------------------------------------------------------#












