#-------------------------------------------------------------------------------------------#
#  Author      : Andrew J Skelton                                                           |
#  Language    : R Statistical Programming Language                                         |
#  Optimised   : CD4 T-Cell RA Microarray Data                                              |
#  Data Owner  : Newcastle University - Dr. A. Pratt                                        |
#  Description : A Microarray analysis pipeline for IlluminaHT-12 array data                |
#                Consisting of two distinct batch corrections.                              |
#-------------------------------------------------------------------------------------------#

##'MONOLITH_DIR
setwd("/data/customers/Musculoskeletal/Pratt/Rheumatoid_Arthritis_CD4_T-Cell_Signature//Raw_Data")
list.files()
##'RAPTOR_DIR
setwd("/Users/andrew/Documents/Bioinformatics/Customers/Pratt/Rheumatoid_Arthritis_CD4_T-Cell_Signature/Raw_Data")
list.files()
raw_files <- c("Raw Data Phase I.txt", "Raw Data Phase II.txt")

##'PACKAGES---------------------------------------------------------------------------------------------------#
install.packages("gplots")
source("http://bioconductor.org/biocLite.R")
biocLite()
biocLite("lumi")
biocLite("arrayQualityMetrics")
biocLite("lumiHumanAll.db")
biocLite("annotate")
biocLite("limma")
biocLite("GOstats")
biocLite("RamiGO")
biocLite("pathview")
biocLite("sva")

library(lumi)
library(gplots)
library(arrayQualityMetrics)
library(lumiHumanAll.db)
library(annotate)
library(limma)
library(ggplot2)
library(GOstats)
library(RamiGO)
library(GOstats)
library(RamiGO)
library(pathview)
library(sva)
##'END_PACKAGES-----------------------------------------------------------------------------------------------#

##'READ_IN_AND_NORMALISATION----------------------------------------------------------------------------------#
##'REMOVE ARRAY 87 - DECISION MADE BY DR. A. PRATT [DATA OWNER]
raw_data            <- lumiR.batch(raw_files)
removal_elements    <- c(87)
removal             <- sampleNames(raw_data)[removal_elements]
raw_data            <- raw_data[, !(sampleNames(raw_data) %in% removal)]
vst_data            <- lumiT(raw_data, method="vst")
#rsn_data            <- lumiN(vst_data, method="rsn")
##'QUANTILE_NORMALISATION
quantile_data      <- lumiN(vst_data, method="quantile")
analysis_ready_data <- lumiQ(quantile_data)
#analysis_ready_data <- lumiQ(rsn_data)
##'END_READ_IN_AND_NORMALISATION------------------------------------------------------------------------------#

##'APPLY_ARRAY_NAMES------------------------------------------------------------------------------------------#
new_Array_Names <- c("177_D_I_i_T_Mn",   "304_C_I_i_T_Mn",   "182_B_I_i_UA_Mn",  "174_A_I_i_UA_Mn",
                     "239_A_I_i_UA_Mp",  "264_B_I_i_T_Mp",   "151_C_I_i_T_Mn",   "175_D_I_i_T_Mn",
                     "259_D_I_i_UA_Mn",  "314_C_I_i_T_Mn",   "214_B_I_i_T_Mn",   "322_A_I_i_T_Mn",
                     "273_A_I_i_T_Mp",   "315_B_I_i_T_Mn",   "307_C_I_i_T_Mn",   "277_D_I_i_UA_Mn",
                     "192_D_I_i_T_Mn",   "193_C_I_i_UA_Mn",  "236_B_I_i_UA_Mp",  "184_B_I_i_UA_Mp",
                     "241_A_I_i_T_Mp",   "221_B_I_i_T_Mp",   "263_C_I_i_T_Mn",   "270_D_I_ii_T_Mn",
                     "210_C_I_ii_UA_Mn", "153_C_I_ii_T_Mn",  "195_B_I_ii_UA_Mn", "253_A_I_ii_T_Mn",
                     "254_A_I_ii_UA_Mp", "294_B_I_ii_UA_Mp", "231_C_I_ii_UA_Mn", "230_D_I_ii_T_Mn",
                     "248_D_I_ii_T_Mn",  "191_C_I_ii_UA_Mn", "218_B_I_ii_T_Mp",  "196_A_I_ii_UA_Mn",
                     "246_B_I_ii_T_Mp",  "149_C_I_ii_T_Mn",  "237_C_I_ii_T_Mn",  "156_D_I_ii_T_Mn",
                     "260_D_I_ii_T_Mn",  "181_C_I_ii_UA_Mn", "256_B_I_ii_T_Mp",  "220_A_I_ii_T_Mp",
                     "305_A_I_ii_UA_Mn", "133_B_I_ii_UA_Mp", "272_C_I_ii_T_Mn",  "157_D_I_iii_T_Mn",
                     "291_D_I_iii_UA_Mn","240_C_I_iii_T_Mn", "268_B_I_iii_T_Mn", "302_A_I_iii_T_Mn",
                     "243_C_I_iii_T_Mn", "125_C_I_iii_T_Mn", "188_C_I_iii_T_Mn", "247_B_I_iii_T_Mn",
                     "269_B_I_iii_T_Mn", "209_A_I_iii_UA_Mn","142_A_I_iii_UA_Mp","261_B_I_iii_UA_Mp",
                     "303_C_I_iii_T_Mn", "211_B_I_iii_UA_Mp","176_A_I_iii_UA_Mn","251_B_I_iii_T_Mn",
                     "250_C_I_iii_T_Mn", "306_B_I_iii_T_Mp", "292_D_I_iii_UA_Mn","172_B_I_iii_T_Mp",
                     "134_C_I_iii_T_Mn", "137_C_I_iii_T_Mn", "204_B_I_iii_T_Mp", "266_C_I_iii_T_Mn",
                     "183_B_I_iv_UA_Mp", "120_A_I_iv_T_Mn",  "144_C_I_iv_T_Mn",  "145_C_I_iv_UA_Mn",
                     "207_D_I_iv_T_Mn",  "141_D_I_iv_T_Mn",  "118_C_I_iv_UA_Mn", "249_B_I_iv_T_Mn",
                     "203_B_I_iv_T_Mp",  "131_A_I_iv_T_Mp",  "299_A_I_iv_UA_Mn", "163_B_I_iv_UA_Mn",
                     "170_C_I_iv_T_Mn",  "159_D_I_iv_T_Mn",  "161_D_I_iv_T_Mn",  "136_C_I_iv_T_Mn",
                     "300_D_I_iv_T_Mn",  "232_A_I_iv_UA_Mn", "288_B_I_iv_T_Mn",  "296_A_I_iv_T_Mn",
                     "212_D_I_iv_UA_Mn", "132_C_I_iv_T_Mn",  "135_C_I_iv_T_Mn",  "388_A_II_v_T_Mn",
                     "297_C_II_v_T_Mn",  "154_C_II_v_UA_Mn", "333_B_II_v_UA_Mp", "364_A_II_v_T_Mn",
                     "316_C_II_v_T_Mn",  "336_D_II_v_T_Mn",  "227_E_II_v_UA_Mn", "442_C_II_v_UA_Mn",
                     "358_B_II_v_UA_Mn", "155_E_II_v_T_Mn",  "382_E_II_v_T_Mn",  "284_E_II_v_T_Mn",
                     "355_E_II_v_T_Mn",  "242_E_II_v_T_Mn",  "337_E_II_v_T_Mn",  "205_E_II_v_T_Mn",
                     "286_E_II_v_UA_Mn", "165_D_II_v_T_Mn",  "301_B_II_v_T_Mp",  "258_D_II_v_T_Mn",
                     "430_E_II_v_T_Mn",  "323_D_II_v_T_Mn",  "401_A_II_v_T_Mn",  "349_B_II_v_UA_Mn",
                     "441_D_II_v_T_Mn",  "391_C_II_v_T_Mn",  "431_C_II_v_UA_Mn", "354_C_II_v_UA_Mn",
                     "126_D_II_v_UA_Mn", "213_C_II_v_UA_Mn", "262_E_II_v_T_Mn",  "238_C_II_v_T_Mn",
                     "197_B_II_v_T_Mp",  "367_E_II_v_T_Mn",  "310_E_II_v_T_Mn",  "255_E_II_vi_T_Mn",
                     "311_C_II_vi_UA_Mn","397_A_II_vi_T_Mn", "283_D_II_vi_T_Mn", "392_D_II_vi_UA_Mn",
                     "366_A_II_vi_UA_Mn","353_C_II_vi_UA_Mn","346_D_II_vi_T_Mn", "309_E_II_vi_T_Mn",
                     "287_B_II_vi_UA_Mn","282_C_II_vi_UA_Mn","334_E_II_vi_UA_Mn","344_A_II_v_UA_Mp",
                     "313_C_II_v_T_Mn",  "340_D_II_v_T_Mn",  "200_E_II_v_UA_Mn", "387_C_II_v_T_Mn",
                     "406_C_II_v_T_Mn",  "433_D_II_v_UA_Mn", "276_A_II_v_T_Mn",  "389_A_II_v_UA_Mn",
                     "369_A_II_v_T_Mp",  "359_B_II_v_UA_Mp", "285_E_II_v_T_Mn",  "293_C_II_v_UA_Mn",
                     "244_A_II_v_T_Mn",  "393_B_II_v_T_Mn",  "439_A_II_v_T_Mp",  "326_E_II_v_T_Mn",
                     "412_C_II_v_UA_Mn", "350_A_II_v_T_Mn",  "405_C_II_v_UA_Mn", "198_D_II_v_UA_Mn",
                     "351_C_II_v_UA_Mn", "320_C_II_v_T_Mn",  "378_B_II_v_T_Mn",  "223_C_II_v_UA_Mn",
                     "374_B_II_v_T_Mp",  "173_E_II_v_T_Mn",  "328_D_II_v_UA_Mn", "216_D_II_v_T_Mn",
                     "222_B_II_v_T_Mn")
sampleNames(raw_data) <- new_Array_Names
##'----------------------------------------------------------------------------------------------------------#

##'QC_FUNCTION - GRAPHS--------------------------------------------------------------------------------------#
##'----------------------------------------------------------------------------------------------------------#

##'ARRAY_REMOVAL_&_RE-NORMALISATION--------------------------------------------------------------------------#
removal   <- c("241_A_I_i_T_Mp",  "133_B_I_ii_UA_Mp", "157_D_I_iii_T_Mn", "142_A_I_iii_UA_Mp",
               "258_D_II_v_T_Mn", "126_D_II_v_UA_Mn", "213_C_II_v_UA_Mn", "238_C_II_v_T_Mn",
               "328_D_II_v_UA_Mn","222_B_II_v_T_Mn" )
raw_data  <- raw_data[, !(sampleNames(raw_data) %in% removal)]
length(sampleNames(raw_data))
vst_data            <- lumiT(raw_data, method="vst")
#rsn_data            <- lumiN(vst_data, method = "rsn")
##'QUANTILE_NORMALISATION
quantile_data      <- lumiN(vst_data, method="quantile")
analysis_ready_data <- lumiQ(quantile_data)
#analysis_ready_data <- lumiQ(rsn_data)
##'----------------------------------------------------------------------------------------------------------#

##'DETECTION_CALL_THRESHOLD_REMOVAL--------------------------------------------------------------------------#
exprs_data                   <- exprs(analysis_ready_data)
present_count                <- detectionCall(analysis_ready_data)
filtered_analysis_ready_data <- exprs_data[present_count > 0, ]
##'----------------------------------------------------------------------------------------------------------#

##'PHASE_AND_AMPLIFICATION_BATCH_CLASSIFICATION--------------------------------------------------------------#
phase_batch <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,
                 2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,
                 2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
amp_batch   <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,
                 2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,
                 4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
                 5,5,5,5,5,5,5,5,5,5,5,5,6,6,6,6,6,6,6,6,6,6,6,6,5,5,5,5,5,5,5,5,5,5,5,5,5,
                 5,5,5,5,5,5,5,5,5,5,5,5,5,5,5)

amp_batch_1 <- amp_batch[1:sum(phase_batch == 1)]
amp_batch_2 <- amp_batch[(sum(phase_batch == 1)+1):(length(amp_batch))]
amp_batch_2 <- gsub(5, 1, amp_batch_2); amp_batch_2 <- gsub(6, 2, amp_batch_2)
amp_batch_2 <- as.numeric(amp_batch_2)

treatments  <- c("A", "B", "C", "D", "E")
array_names <- c("D","C","B","A","A","B","C","D","D","C","B","A","A","B","C","D","D","C",
                 "B","B","B","C","D","C","C","B","A","A","B","C","D","D","C","B","A","B",
                 "C","C","D","D","C","B","A","A","C","D","C","B","A","C","C","C","B","B",
                 "A","B","C","B","A","B","C","B","D","B","C","C","B","C","B","A","C","C",
                 "D","D","C","B","B","A","A","B","C","D","D","C","D","A","B","A","D","C",
                 "C","A","C","C","B","A","C","D","E","C","B","E","E","E","E","E","E","E",
                 "E","D","B","E","D","A","B","D","C","C","C","E","B","E","E","E","C","A",
                 "D","D","A","C","D","E","B","C","E","A","C","D","E","C","C","D","A","A",
                 "A","B","E","C","A","B","A","E","C","A","C","D","C","C","B","C","B","E",
                 "D")
##'----------------------------------------------------------------------------------------------------------#

##'PHASE_BATCH_SEPARATION------------------------------------------------------------------------------------#
phase_1 <- filtered_analysis_ready_data[,1:length(amp_batch_1)]
phase_2 <- filtered_analysis_ready_data[,(length(amp_batch_1)+1):length(phase_batch)]
##'----------------------------------------------------------------------------------------------------------#

##'AMPLIFICATION_BATCH_CORRECTION----------------------------------------------------------------------------#
pheno            <- data.frame(sample = c(1:length(amp_batch_1)),
                               outcome = array_names[1:length(amp_batch_1)],
                               batch = amp_batch_1)
rownames(pheno)  <- colnames(phase_1)
batch                = pheno$batch
mod                  = model.matrix(~as.factor(outcome), data=pheno)
phase_1_combat_edata = ComBat(dat=phase_1, batch=batch, mod=mod, numCovs=NULL,
                              par.prior=TRUE, prior.plots=FALSE)

pheno            <- data.frame(sample = c(1:length(amp_batch_2)),
                               outcome = array_names[(length(amp_batch_1)+1):length(phase_batch)],
                               batch = amp_batch_2)
rownames(pheno)  <- colnames(phase_2)
batch                = pheno$batch
mod                  = model.matrix(~as.factor(outcome), data=pheno)
phase_2_combat_edata = ComBat(dat=phase_2, batch=batch, mod=mod, numCovs=NULL,
                              par.prior=TRUE, prior.plots=FALSE)
##'----------------------------------------------------------------------------------------------------------#

##'POST_AMPLIFICATION_CORRECTION_JOIN DATASETS---------------------------------------------------------------#
merged_data_post_amp <- cbind(phase_1_combat_edata, phase_2_combat_edata)
##'----------------------------------------------------------------------------------------------------------#

##'PHASE_BATCH_CORRECTION------------------------------------------------------------------------------------#
pheno            <- data.frame(sample = c(1:163), outcome = array_names,
                               batch = phase_batch)
rownames(pheno)  <- colnames(merged_data_post_amp)
batch                  = pheno$batch
mod                    = model.matrix(~as.factor(outcome), data=pheno)
merged_data_post_phase = ComBat(dat=merged_data_post_amp, batch=batch,
                                mod=mod, numCovs=NULL, par.prior=TRUE,
                                prior.plots=FALSE)
##'----------------------------------------------------------------------------------------------------------#

##PROBE_ANNOTATION-------------------------------------------------------------------------------------------#
probe_list       <- rownames(merged_data_post_phase)
nuIDs            <- probeID2nuID(probe_list)[, "nuID"]
nuIDs[182]       <- probeID2nuID(3450064)[, "nuID"]

getSYMBOL(probeID2nuID(3840358)[, "nuID"], "lumiHumanAll.db")
lookUp(probeID2nuID(3840358)[, "nuID"], "lumiHumanAll.db", "GENENAME")

symbol           <- getSYMBOL(nuIDs, "lumiHumanAll.db")
name             <- unlist(lookUp(nuIDs, "lumiHumanAll.db", "GENENAME"))
anno_df          <- data.frame(ID = nuIDs, probe_list, symbol, name)

design           <- model.matrix(~0 + factor(array_names, levels = treatments))
colnames(design) <- treatments
num_parameters   <- ncol(design)
fit              <- lmFit(merged_data_post_phase, design)

cont_mat         <- makeContrasts(E-A, E-B, E-C, E-D,
                                  levels=treatments)
fit2             <- contrasts.fit(fit, contrasts=cont_mat)
fit2             <- eBayes(fit2)
fit2$genes = anno_df
##'----------------------------------------------------------------------------------------------------------#

##PROBE_ANNOTATION_TRAINING_COHORT---------------------------------------------------------------------------#

merged_data_post_phase_T <- merged_data_post_phase[,grep("_T", colnames(merged_data_post_phase))]

probe_list       <- rownames(merged_data_post_phase_T)
nuIDs            <- probeID2nuID(probe_list)[, "nuID"]
nuIDs[grep("NA", nuIDs)]  <- probeID2nuID(grep("NA", nuIDs))[, "nuID"]

getSYMBOL(probeID2nuID(3840358)[, "nuID"], "lumiHumanAll.db")
lookUp(probeID2nuID(3840358)[, "nuID"], "lumiHumanAll.db", "GENENAME")

getSYMBOL(probeID2nuID(2060615)[, "nuID"], "lumiHumanAll.db")
lookUp(probeID2nuID(2060615)[, "nuID"], "lumiHumanAll.db", "GENENAME")
symbol_2 <- c()
for( i in 1:length(nuIDs) )
{
  symbol_2<- c(symbol_2, getSYMBOL(nuIDs[i], "lumiHumanAll.db"))
}

nuIDs_2 <- c()
for( i in 1:length(rownames(merged_data_post_phase_T)) )
{
  nuIDs_2 <- c(nuIDs_2, probeID2nuID(rownames(merged_data_post_phase_T)[i])[, "nuID"])
}

symbol           <- getSYMBOL(nuIDs, "lumiHumanAll.db")
name             <- unlist(lookUp(nuIDs, "lumiHumanAll.db", "GENENAME"))
anno_df          <- data.frame(ID = nuIDs, probe_list, symbol, name)

treatments_T         <- c("AB", "CD", "E")
array_names_T      <- array_names[grep("_T", colnames(merged_data_post_phase_T))]
array_names_T[array_names_T == "B" | array_names_T == "A"] <- "AB"
array_names_T[array_names_T == "C" | array_names_T == "D"] <- "CD"
design_T           <- model.matrix(~0 + factor(array_names_T, levels = treatments_T))
colnames(design_T) <- treatments_T
num_parameters     <- ncol(design_T)
fit_T              <- lmFit(merged_data_post_phase_T, design_T)

cont_mat_T         <- makeContrasts(AB-CD,levels=treatments_T)
fit2_T             <- contrasts.fit(fit_T, contrasts=cont_mat_T)
fit2_T             <- eBayes(fit2_T)
fit2_T$genes = anno_df ##Anno_Df
##'----------------------------------------------------------------------------------------------------------#

##TOP_TABLES-------------------------------------------------------------------------------------------------#
topTable(fit2_T, coef="AB - CD", p.value=0.05, lfc=log2(1.2), adjust.method="BH")
gene_list_Coef_A_unfiltered <- topTable(fit2_T, coef = "AB - CD", number = nrow(anno_df))
gene_list_Coef_A <- gene_list_Coef_A_unfiltered[abs(gene_list_Coef_A_unfiltered$logFC) > log2(1.2) &
                                                  gene_list_Coef_A_unfiltered$P.Value < 0.05,]

controls      <- merged_data_post_phase_T[as.character(gene_list_Coef_A$probe_list),
                                          grep("AB", array_names_T)]
conditions    <- merged_data_post_phase_T[as.character(gene_list_Coef_A$probe_list),
                                          grep("CD", array_names_T)]

gene_list_Coef_A$ControlMean   <- apply(controls, 1, mean)
gene_list_Coef_A$ConditionMean <- apply(conditions, 1, mean)
##'----------------------------------------------------------------------------------------------------------#

##VOLCANO_PLOT-----------------------------------------------------------------------------------------------#
gene_list_Coef_A_unfiltered$threshold = as.factor(abs(gene_list_Coef_A_unfiltered$logFC) > log2(1.2) &
                                                    gene_list_Coef_A_unfiltered$P.Value < 0.05)
g = ggplot(data=gene_list_Coef_A_unfiltered, aes(x=gene_list_Coef_A_unfiltered$logFC,
                                                 y=-log10(gene_list_Coef_A_unfiltered$P.Value),
                                                 colour=gene_list_Coef_A_unfiltered$threshold)) +
  geom_point(alpha=0.4, size=1.75) +
  theme(legend.position = "none") +
  #xlim(c(-2, 2)) + ylim(c(0, 1)) +
  xlab("log2 fold change") + ylab("-log10 p-value") +
  ggtitle("Volcano Plot of Fold Change vs P-Value (No Multiple Test Correction)")
g
##'----------------------------------------------------------------------------------------------------------#
