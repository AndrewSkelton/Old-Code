#-------------------------------------------------------------------------------------------#
#  Author      : Andrew J Skelton                                                           |
#  Language    : R Statistical Programming Language                                         |
#  Optimised   : CD4 T-Cell RA Microarray Data                                              |
#  Data Owner  : Newcastle University - Dr. A. Pratt                                        |
#  Description : A Microarray analysis pipeline for downstream analysis of a                |
#                CSV file                                                                   |
#-------------------------------------------------------------------------------------------#

library(stringr)
library(sva)
library(lumi)
library(gplots)
library(annotate)
library(lumiHumanAll.db)
library(limma)

##LOAD_DATA_&_PREPARE----------------------------------------------------------------------------------------#
setwd("/Users/andrew/Documents/Bioinformatics/Customers/Pratt/Rheumatoid_Arthritis_CD4_T-Cell_Signature/Raptor")
list.files()
raw_data <- read.csv("phase_Correction_R_readin.csv", strip.white=TRUE, header=T, stringsAsFactors=FALSE)
rownames(raw_data) <- raw_data[,1]
raw_data <- raw_data[,-1]
colnames(raw_data) <- substring(colnames(raw_data), 2)
##'----------------------------------------------------------------------------------------------------------#

phase <- colnames(raw_data)
phase[grep("_1_", colnames(raw_data))] <- 1
phase[grep("_2_", colnames(raw_data))] <- 2
phase

##'BATCH CORRECTION (PHASE)----------------------------------------------------------------------------------#
#'hard coded based on sample naming convention
array_names      <- str_sub(colnames(raw_data), -5, -5)
pheno            <- data.frame(sample = c(1:length(colnames(raw_data))),
                               outcome = array_names, batch = phase)
rownames(pheno)  <- colnames(raw_data)
batch            <- pheno$batch
mod              <- model.matrix(~as.factor(outcome), data = pheno)
merged_data_post_phase = ComBat(dat = raw_data, batch = batch, mod = mod, numCovs = NULL,
                                par.prior = TRUE, prior.plots = FALSE)
##'----------------------------------------------------------------------------------------------------------#

##'PCA_PLOTS-------------------------------------------------------------------------------------------------#
#PCA_data <- data.matrix(raw_data)
PCA_data <- data.matrix(merged_data_post_phase)
colnames(PCA_data) <- phase
plotSampleRelation(PCA_data, method= 'mds' ,     cv.Th=0, main="")
#plotSampleRelation(PCA_data, method= 'cluster' , cv.Th=0, main="", cex=.4)
##'----------------------------------------------------------------------------------------------------------#

##COHORT_SUBSET----------------------------------------------------------------------------------------------#
normalised_dataset_T <- merged_data_post_phase[,grep("_T_", colnames(merged_data_post_phase))]
##'----------------------------------------------------------------------------------------------------------#

##PROBE_ANNOTATION-------------------------------------------------------------------------------------------#
probe_list       <- rownames(normalised_dataset_T)
nuIDs            <- probeID2nuID(probe_list)[, "nuID"]
which(is.na(nuIDs)) #321
nuIDs[321]       <- probeID2nuID(3450064)[, "nuID"]
symbol           <- getSYMBOL(nuIDs, "lumiHumanAll.db")
name             <- unlist(lookUp(nuIDs, "lumiHumanAll.db", "GENENAME"))
anno_df          <- data.frame(ID = nuIDs, probe_list, symbol, name)

treatments  <- c("AB", "CD")
array_names <- colnames(normalised_dataset_T)
array_names[grep("_A_|_B_", colnames(normalised_dataset_T))] <- "AB"
array_names[grep("_C_|_D_", colnames(normalised_dataset_T))] <- "CD"

normalised_dataset_T <- data.matrix(normalised_dataset_T)
design           <- model.matrix(~0 + factor(array_names, levels = treatments))
colnames(design) <- treatments
num_parameters   <- ncol(design)
fit              <- lmFit(normalised_dataset_T, design)

cont_mat         <- makeContrasts(AB-CD, levels=treatments)
fit2             <- contrasts.fit(fit, contrasts=cont_mat)
fit2             <- eBayes(fit2)
fit2$genes = anno_df
##'----------------------------------------------------------------------------------------------------------#


##TOP_TABLES-------------------------------------------------------------------------------------------------#
topTable(fit2, coef="AB - CD", p.value=0.05, lfc=log2(1.2), adjust.method="BH")
gene_list_Coef_A_unfiltered <- topTable(fit2, coef = "AB - CD", number = nrow(anno_df), adjust.method="none")
gene_list_Coef_A <- gene_list_Coef_A_unfiltered[abs(gene_list_Coef_A_unfiltered$logFC) > log2(1.2) &
                                                gene_list_Coef_A_unfiltered$P.Value < 0.05,]

controls      <- normalised_dataset_T[as.character(gene_list_Coef_A$probe_list),
                                        grep("AB", array_names)]
conditions    <- normalised_dataset_T[as.character(gene_list_Coef_A$probe_list),
                                        grep("CD", array_names)]

gene_list_Coef_A$ControlMean   <- apply(controls,   1, mean)
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
    xlim(c(-4, 4)) + ylim(c(0, 4.5)) +
    xlab("log2 fold change") + ylab("-log10 p-value") +
    ggtitle("Volcano Plot of Fold Change vs P-Value (No Multiple Test Correction)")
  g
##'----------------------------------------------------------------------------------------------------------#

write.csv(gene_list_Coef_A, "gene_list_post_phase_correction.csv")
