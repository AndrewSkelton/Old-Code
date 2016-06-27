#-------------------------------------------------------------------------------------------#
#  Author      : Andrew J Skelton                                                           |
#  Language    : R Statistical Programming Language                                         |
#  Optimised   : Rheumatoid Arthritis CD4 T-Cell Bio Signature                              |
#  Data Owner  : Newcastle University - Arthur Pratt                                        |
#-------------------------------------------------------------------------------------------#

source("http://bioconductor.org/biocLite.R")
biocLite()
biocLite(c("lumiHumanAll.db", "RamiGO", "pathview"))
library(stringr)
library(sva)
library(lumi)
library(gplots)
library(ggplot2)
library(annotate)
library(lumiHumanIDMapping)
library(lumiHumanAll.db)
library(limma)
library(pathview)
library(GOstats)
library(RamiGO)
library(affycoretools)

##LOAD_DATA_&_PREPARE----------------------------------------------------------------------------------------#
setwd("/Users/andrew/Documents/Bioinformatics/Customers/Pratt/Illumina_Arrays_July_2014/Raw_data/")
list.files()
filename <- "pratt_all_sample_probe_profile.txt"

raw_data <- lumiR(filename, verbose=T)
raw_data <- addControlData2lumi("control_probe_profile.txt", raw_data)

det <- detection(raw_data)
png("Raw_Detection_PValues.png",width=12.53,height=6.98,units="in",res=600) 
  boxplot((det), main="Distribution of Detection P Values", xaxt="n", cex=.2)
dev.off()

png("Raw_Intensities.png",width=12.53,height=6.98,units="in",res=600) 
  density(raw_data, legend=NA, main="Density of log2 intensities")
dev.off()

png("879_948_Intensities.png",width=12.53,height=6.98,units="in",res=600) 
  density(raw_data[,115:116], legend=NA, main="Density of log2 intensities - 879|948")
dev.off()

png("879_948_Removed_Intensities.png",width=12.53,height=6.98,units="in",res=600) 
  density(raw_data_det, legend=NA, main="Density of log2 intensities - 879|948 Removed")
dev.off()

array_names <- c("1072_RA_9477874134_A_SB2_NR", "1085_UA_9477874134_B_SB2_ND", "1076_NRA_9477874134_C_SB2_ND", 
                 "1086_RA_9477874134_D_SB2_GR", "1079_NRA_9477874134_E_SB2_ND", "1087_NRA_9477874134_F_SB2_ND", 
                 "1080_NRA_9477874134_G_SB2_ND", "1088_UA_9477874134_H_SB2_ND", "1083_RA_9477874134_I_SB2_ND", 
                 "1094_NRA_9477874134_J_SB2_NR", "1084_NRA_9477874134_K_SB2_ND", "1096_NRA_9477874134_L_SB2_ND", 
                 "1045_NRA_9477874138_A_SB2_ND", "1059_NRA_9477874138_B_SB2_ND", "1050_NRA_9477874138_C_SB2_ND", 
                 "1067_RA_9477874138_D_SB2_NR", "1051_RA_9477874138_E_SB2_GR", "1068_NRA_9477874138_F_SB2_ND", 
                 "1054_UA_9477874138_G_SB2_ND", "1069_NRA_9477874138_H_SB2_ND", "1056_NRA_9477874138_I_SB2_ND", 
                 "1070_NRA_9477874138_J_SB2_ND", "1058_UA_9477874138_K_SB2_ND", "1071_NRA_9477874138_L_SB2_ND", 
                 "1097_RA_9477874140_A_SB2_ND", "2019_NRA_9477874140_B_SB2_ND", "2002_NRA_9477874140_C_SB2_ND", 
                 "2021_NRA_9477874140_D_SB2_ND", "2010_RA_9477874140_E_SB2_ND", "2022_NRA_9477874140_F_SB2_ND", 
                 "2012_NRA_9477874140_G_SB2_ND", "2025_NRA_9477874140_H_SB2_ND", "2013_RA_9477874140_I_SB2_ND", 
                 "2028_UA_9477874140_J_SB2_ND", "2014_NRA_9477874140_K_SB2_ND", "2029_RA_9477874140_L_SB2_GR", 
                 "2055_NRA_9477874145_A_SB3_ND", "2075_NRA_9477874145_B_SB3_ND", "2062_NRA_9477874145_C_SB3_ND", 
                 "2078_NRA_9477874145_D_SB3_ND", "2063_NRA_9477874145_E_SB3_ND", "2086_UA_9477874145_F_SB3_ND", 
                 "2065_NRA_9477874145_G_SB3_ND", "2087_NRA_9477874145_H_SB3_ND", "2067_UA_9477874145_I_SB3_ND", 
                 "2088_NRA_9477874145_J_SB3_ND", "2072_RA_9477874145_K_SB3_NR", "2090_RA_9477874145_L_SB3_ND", 
                 "2094_NRA_9477874149_A_SB3_ND", "2106_RA_9477874149_B_SB3_NR", "2096_NRA_9477874149_C_SB3_ND", 
                 "2110_RA_9477874149_D_SB3_MR", "2099_NRA_9477874149_E_SB3_ND", "2124_NRA_9477874149_F_SB3_ND", 
                 "2101_NRA_9477874149_G_SB3_ND", "2128_NRA_9477874149_H_SB3_ND", "2103_NRA_9477874149_I_SB3_ND", 
                 "2131_NRA_9477874149_J_SB3_ND", "2104_NRA_9477874149_K_SB3_ND", "2135_NRA_9477874149_L_SB3_ND", 
                 "950_NRA_9477874152_A_SB5_ND", "1019_UA_9477874152_B_SB5_ND", "997_RA_9477874152_C_SB5_NR", 
                 "1020_UA_9477874152_D_SB5_ND", "1003_RA_9477874152_E_SB5_ND", "1028_NRA_9477874152_F_SB5_ND", 
                 "1006_NRA_9477874152_G_SB5_ND", "1029_RA_9477874152_H_SB5_GR", "1008_NRA_9477874152_I_SB5_ND", 
                 "1033_NRA_9477874152_J_SB5_ND", "1018_NRA_9477874152_K_SB5_ND", "1036_NRA_9477874152_L_SB5_ND", 
                 "2030_RA_9477874158_A_SB2_GR", "2044_NRA_9477874158_B_SB2_ND", "2033_NRA_9477874158_C_SB2_ND", 
                 "2045_RA_9477874158_D_SB2_GR", "2034_UA_9477874158_E_SB2_ND", "2047_NRA_9477874158_F_SB2_ND", 
                 "2036_NRA_9477874158_G_SB2_ND", "2052_UA_9477874158_H_SB2_ND", "2040_UA_9477874158_I_SB2_ND", 
                 "2053_NRA_9477874158_J_SB2_ND", "2042_RA_9477874158_K_SB2_NR", "2054_NRA_9477874158_L_SB2_ND", 
                 "2140_UA_9477874168_A_SB3_ND", "2154_RA_9477874168_B_SB3_ND", "2144_NRA_9477874168_C_SB3_ND", 
                 "2222_NRA_9477874168_D_SB3_ND", "2145_NRA_9477874168_E_SB3_ND", "2231_RA_9477874168_F_SB3_GR", 
                 "2146_UA_9477874168_G_SB3_ND", "2233_NRA_9477874168_H_SB3_ND", "2148_NRA_9477874168_I_SB3_ND", 
                 "2234_UA_9477874168_J_SB3_ND", "2149_NRA_9477874168_K_SB3_ND", "2255_RA_9477874168_L_SB3_NR", 
                 "811_NRA_9477874172_A_SB4_ND", "827_RA_9477874172_B_SB4_GR", "812_NRA_9477874172_C_SB4_ND", 
                 "828_NRA_9477874172_D_SB4_ND", "813_UA_9477874172_E_SB4_ND", "832_NRA_9477874172_F_SB4_ND", 
                 "814_UA_9477874172_G_SB4_ND", "833_NRA_9477874172_H_SB4_ND", "815_UA_9477874172_I_SB4_ND", 
                 "834_RA_9477874172_J_SB4_MR", "826_UA_9477874172_K_SB4_ND", "836_NRA_9477874172_L_SB4_ND", 
                 "1037_NRA_9477874174_A_SB5_ND", "985_NRA_9477874174_B_SB5_ND", "857_NRA_9477874174_C_SB5_ND", 
                 "1005_RA_9477874174_D_SB5_ND", "873_NRA_9477874174_E_SB5_ND", "1014_NRA_9477874174_F_SB5_ND", 
                 "879_UA_9477874174_G_SB5_ND", "948_UA_9477874174_H_SB5_GR", "932_RA_9477874174_I_SB5_ND", 
                 "984_NRA_9477874174_J_SB5_ND", "944_UA_9477874174_K_SB5_ND", "1010_RA_9477874174_L_SB5_NR", 
                 "837_UA_9477874177_A_SB4_ND", "844_RA_9477874177_B_SB4_NR", "838_RA_9477874177_C_SB4_GR", 
                 "849_NRA_9477874177_D_SB4_ND", "839_UA_9477874177_E_SB4_ND", "855_NRA_9477874177_F_SB4_ND", 
                 "840_RA_9477874177_G_SB4_GR", "914_NRA_9477874177_H_SB4_ND", "841_UA_9477874177_I_SB4_ND", 
                 "915_UA_9477874177_J_SB4_MR", "842_NRA_9477874177_K_SB4_ND", "925_NRA_9477874177_L_SB4_ND", 
                 "2257_UA_9477874180_A_SB3_ND", "835_RA_9477874180_B_SB3_ND", "2261_NRA_9477874180_C_SB3_ND", 
                 "874_UA_9477874180_D_SB3_ND", "2264_NRA_9477874180_E_SB3_ND", "878_NRA_9477874180_F_SB3_ND", 
                 "808_NRA_9477874180_G_SB3_ND", "969_UA_9477874180_H_SB3_ND", "829_NRA_9477874180_I_SB3_ND", 
                 "809_NRA_9477874180_J_SB3_ND", "830_NRA_9477874180_K_SB3_ND", "810_RA_9477874180_L_SB3_GR", 
                 "1032_RA_9477874192_A_SB5_MR", "835rpt_RA_9477874192_B_SB5_ND", "1040_NRA_9477874192_C_SB5_ND", 
                 "836rpt_NRA_9477874192_D_SB5_ND", "1042_RA_9477874192_E_SB5_GR", "837rpt_UA_9477874192_F_SB5_ND", 
                 "832rpt_NRA_9477874192_G_SB5_ND", "838rpt_RA_9477874192_H_SB5_GR", "833rpt_NRA_9477874192_I_SB5_ND", 
                 "839rpt_UA_9477874192_J_SB5_ND", "834rpt_RA_9477874192_K_SB5_MR", "840rpt_RA_9477874192_L_SB5_GR", 
                 "938_UA_9479187009_A_SB4_ND", "1012_NRA_9479187009_B_SB4_ND", "954_UA_9479187009_C_SB4_ND", 
                 "1015_NRA_9479187009_D_SB4_ND", "962_NRA_9479187009_E_SB4_ND", "1022_NRA_9479187009_F_SB4_ND", 
                 "965_UA_9479187009_G_SB4_ND", "941_RA_9479187009_H_SB4_MR", "968_RA_9479187009_I_SB4_ND", 
                 "995_UA_9479187009_J_SB4_NR", "1000_UA_9479187009_K_SB4_ND", "996_UA_9479187009_L_SB4_NR", 
                 "843_NRA_9479628113_A_SB1_ND", "861_NRA_9479628113_B_SB1_ND", "845_NRA_9479628113_C_SB1_ND", 
                 "862_NRA_9479628113_D_SB1_ND", "850_NRA_9479628113_E_SB1_ND", "863_NRA_9479628113_F_SB1_ND", 
                 "854_RA_9479628113_G_SB1_ND", "881_NRA_9479628113_H_SB1_ND", "856_NRA_9479628113_I_SB1_ND", 
                 "882_RA_9479628113_J_SB1_GR", "860_NRA_9479628113_K_SB1_ND", "883_NRA_9479628113_L_SB1_ND", 
                 "884_UA_9479628130_A_SB1_ND", "897_NRA_9479628130_B_SB1_ND", "890_NRA_9479628130_C_SB1_ND", 
                 "898_UA_9479628130_D_SB1_ND", "891_UA_9479628130_E_SB1_ND", "899_NRA_9479628130_F_SB1_ND", 
                 "892_RA_9479628130_G_SB1_NR", "905_NRA_9479628130_H_SB1_ND", "893_NRA_9479628130_I_SB1_ND", 
                 "906_NRA_9479628130_J_SB1_ND", "896_RA_9479628130_K_SB1_NR", "912_NRA_9479628130_L_SB1_ND", 
                 "913_NRA_9479628144_A_SB1_ND", "930_RA_9479628144_B_SB1_ND", "923_UA_9479628144_C_SB1_ND", 
                 "931_RA_9479628144_D_SB1_NR", "924_NRA_9479628144_E_SB1_ND", "934_UA_9479628144_F_SB1_ND", 
                 "926_NRA_9479628144_G_SB1_NR", "935_NRA_9479628144_H_SB1_ND", "927_NRA_9479628144_I_SB1_ND", 
                 "937_NRA_9479628144_J_SB1_ND", "929_RA_9479628144_K_SB1_MR", "945_UA_9479628144_L_SB1_ND", 
                 "946_UA_9481417007_A_SB1_ND", "975_RA_9481417007_B_SB1_ND", "957_RA_9481417007_C_SB1_ND", 
                 "978_RA_9481417007_D_SB1_ND", "959_NRA_9481417007_E_SB1_ND", "980_UA_9481417007_F_SB1_ND", 
                 "967_NRA_9481417007_G_SB1_ND", "983_NRA_9481417007_H_SB1_ND", "973_UA_9481417007_I_SB1_ND", 
                 "992_RA_9481417007_J_SB1_GR", "974_RA_9481417007_K_SB1_MR", "1011_NRA_9481417007_L_SB1_ND")

sampleNames(raw_data) <- array_names
grep("rpt", sampleNames(raw_data), value=T)#, 
##'115/116 (879/948) - Technical Failure
##'146/134 (835rpt/835) - Outliers (possible sample contamination)
##'
##'148/150/151/152/153/154/155/156
##'(836rpt/837rpt/832rpt/838rpt/833rpt/839rpt/834rpt/840rpt)
##'Repeat Samples
raw_data_det <- raw_data[, -c(115,116,146,134,148,150,151,152,153,154,155,156)]
##'----------------------------------------------------------------------------------------------------------#

##'NORMALISATION---------------------------------------------------------------------------------------------#
raw_bck          <- lumiB(raw_data_det, method="bgAdjust", verbose=T)
vst_data         <- lumiT(raw_bck, method='vst')
rsn_data         <- lumiN(vst_data, method = "rsn")
lumi.Q           <- lumiQ(rsn_data)
##'QUANTILE_NORMALISATION - COMMENT OUT RSN NORMALISATION LINE IF YOU'RE USING QUANTILE
#quantile_data       <- lumiN(vst_data, method="quantile")
#normalised_data     <- lumiQ(quantile_data)


range_in <- grep('9477874174', colnames(raw_data))
png("array_9477874174_boxplot.png",width=12.53,height=6.98,units="in",res=600) 
  boxplot(raw_data[,range_in], main="Boxplot of Array 9477874174")
dev.off()
density(raw_data[,range_in])

png("density_Plot_Without879|948|835rpt|835.png",width=12.53,height=6.98,units="in",res=600) 
  density(lumi.Q, legend=F, main='Normalised Sample Intensities - Without 879|948|835rpt|835')
dev.off()

density(raw_data_det, legend=F)
density(raw_data, legend=F)
plot(vst_data, what='density')
plot(rsn_data, what='density')
plot(lumi.Q, what='density')

summary(lumi.Q, 'QC')
plot(lumi.Q, what='pair')
plotCDF(lumi.Q, reverse=TRUE)
pairs(lumi.Q)
pairs(lumi., smoothScatter=T)


##'----------------------------------------------------------------------------------------------------------#

##'DETECTION_CALL_THRESHOLD_REMOVAL--------------------------------------------------------------------------#
exprs_data                   <- exprs(lumi.Q)
present_count                <- detectionCall(lumi.Q)
normalised_data              <- exprs_data[present_count > 0, ]
##'----------------------------------------------------------------------------------------------------------#

##'BATCH CORRECTION (PHASE)----------------------------------------------------------------------------------#
# batches <- colnames(normalised_data)
# if(subset == TRUE)
# {
#     batches[grep('_2', batches)] <- 1
#     batches[grep('_3', batches)] <- 2
#     batches[grep('_4', batches)] <- 3
# }else
# {
#     batches[grep('_1', batches)] <- 1
#     batches[grep('_2', batches)] <- 2
#     batches[grep('_3', batches)] <- 3
#     batches[grep('_4', batches)] <- 4
# }
# 
# classification   <- sapply(strsplit(colnames(normalised_data), split='_'), function(x) x[[2]])
# pheno            <- data.frame(sample = c(1:length(colnames(normalised_data))),
#                                outcome = classification, batch = batches)
# rownames(pheno)  <- colnames(normalised_data)
# batch            <- pheno$batch
# mod              <- model.matrix(~as.factor(outcome), data = pheno)
# batchCorrected_data = ComBat(dat = normalised_data, batch = batch, mod = mod, numCovs = NULL,
#                                       par.prior = TRUE, prior.plots = TRUE)
##'----------------------------------------------------------------------------------------------------------#

##PROBE_ANNOTATION-------------------------------------------------------------------------------------------#
probe_list       <- rownames(normalised_data)
nuIDs            <- probeID2nuID(probe_list)[, "nuID"]
symbol           <- getSYMBOL(nuIDs, "lumiHumanAll.db")
name             <- unlist(lookUp(nuIDs, "lumiHumanAll.db", "GENENAME"))
anno_df          <- data.frame(ID = nuIDs, probe_list, symbol, name)
##'----------------------------------------------------------------------------------------------------------#

##EXPERIMENTAL_DESIGN----------------------------------------------------------------------------------------#
treatments  <- c("RA", "NRA", "UA")
array_names <- Treatment

design              <- model.matrix(~0 + factor(array_names, levels = treatments))
colnames(design)    <- treatments
num_parameters      <- ncol(design)
fit                 <- lmFit(normalised_data, design)

cont_mat         <- makeContrasts(RA-NRA, RA-UA, NRA-UA, levels=treatments)
fit2             <- contrasts.fit(fit, contrasts=cont_mat)
fit2             <- eBayes(fit2)
fit2$genes = anno_df
##'----------------------------------------------------------------------------------------------------------#

##GENE_LISTS-------------------------------------------------------------------------------------------------#
comparisons <- c("RA - NRA", "RA - UA", "NRA - UA")
p_cut_off   <- 0.05
fold_change <- 1.2
i           <- 1

topTable(fit2, coef = comparisons[i], p.value = p_cut_off, lfc = log2(fold_change), adjust.method="none")
gene_list_unfiltered <- topTable(fit2, coef = comparisons[i], number = nrow(anno_df), adjust.method="none")
gg_volcano(gene_list_unfiltered, afc=1.2, pval=0.05)


gene_list            <- topTable(fit2, coef = comparisons[i], p.value = p_cut_off, lfc = log2(fold_change),
                                 number = nrow(anno_df), adjust.method = "BH")

comparison_split <- strsplit(comparisons[i], split=" ")
controls      <- batchCorrected_data[as.character(gene_list$probe_list),
                                     grep(comparison_split[[1]][1], array_names)]
conditions    <- batchCorrected_data[as.character(gene_list $probe_list),
                                     grep(comparison_split[[1]][3], array_names)]

gene_list$ControlMean   <- apply(controls,   1, mean)
gene_list$ConditionMean <- apply(conditions, 1, mean)

write.csv(gene_list, paste("gene_list_", comparison_split[[1]][1], "-", comparison_split[[1]][3], ".csv",
                           sep="",collapse=""))
##'----------------------------------------------------------------------------------------------------------#

##'VOLCANO_PLOT----------------------------------------------------------------------------------------------#
gene_list_unfiltered$threshold = as.factor(abs(gene_list_unfiltered$logFC) > log2(fold_change) &
                                             gene_list_unfiltered$adj.P.Val < p_cut_off)
g = ggplot(data=gene_list_unfiltered, aes(x=gene_list_unfiltered$logFC,
                                          y=-log10(gene_list_unfiltered$P.Value),
                                          colour=gene_list_unfiltered$threshold)) +
  geom_point(alpha=0.4, size=1.75) +
  theme(legend.position = "none") +
  #xlim(c(-4, 4)) + ylim(c(0, 4.5)) +
  xlab("log2 fold change") + ylab("-log10 p-value") +
  ggtitle(paste("Volcano Plot of Fold Change vs P-Value", comparisons[i], sep=" "))
print(g)
#}#
#dev.off()#
##'----------------------------------------------------------------------------------------------------------#

##'PATHWAY_ANALYSIS----------------------------------------------------------------------------------------------#
##'TO BE REPLACED BY REACTOME.
sig_probes      <- as.character(gene_list$probe_list)
entrez          <- unique(unlist(lookUp(nuIDs[sig_probes],
                                        "lumiHumanAll.db", "ENTREZID")))
entrez          <- as.character(entrez[!is.na(entrez)])
entrez_universe <- unique(unlist(lookUp(nuIDs,"lumiHumanAll.db", "ENTREZID")))
entrez_universe <- as.character(entrez_universe[!is.na(entrez_universe)])

kegg_params <- new("KEGGHyperGParams",
                   geneIds=entrez,
                   universeGeneIds=entrez_universe,
                   annotation="lumiHumanAll.db",
                   pvalueCutoff= 0.01,
                   testDirection="over")
kegg_hyperg_result <- hyperGTest(kegg_params)
print(kegg_hyperg_result)

sig_logFCs <- gene_list$logFC
fc_entrez <- unlist(lookUp(nuIDs[sig_probes],
                           "lumiHumanAll.db", "ENTREZID"))
names(sig_logFCs) <- fc_entrez
sig_logFCs <- sig_logFCs[!is.na(entrez)]

summary(kegg_hyperg_result)

for(j in 1:length(sigCategories(kegg_hyperg_result)))
{
  pv <- pathview(gene.data=sig_logFCs,
                 pathway.id=sigCategories(kegg_hyperg_result)[j],
                 species="hsa",
                 limit=list(gene=c(round(min(sig_logFCs)),
                                   round(max(sig_logFCs))), cpd=1))
}
##'----------------------------------------------------------------------------------------------------------#

##GO ANNOTATION----------------------------------------------------------------------------------------------#
sig_probes      <- as.character(gene_list$probe_list)
entrez          <- unique(unlist(lookUp(nuIDs[sig_probes],
                                        "lumiHumanAll.db", "ENTREZID")))
entrez          <- as.character(entrez[!is.na(entrez)])
entrez_universe <- unique(unlist(lookUp(nuIDs,"lumiHumanAll.db", "ENTREZID")))
entrez_universe <- as.character(entrez_universe[!is.na(entrez_universe)])

params <- new("GOHyperGParams", geneIds=entrez,universeGeneIds=entrez_universe,annotation="lumiHumanAll.db",
              ontology="BP", pvalueCutoff= 0.01, conditional=FALSE, testDirection="over")
hyperg_result <- hyperGTest(params)
print(hyperg_result)

pval_go     <- pvalues(hyperg_result)
go_fdr      <- p.adjust(pvalues(hyperg_result), method="fdr")
sig_go_id   <- names(go_fdr[go_fdr < 0.01])
if(length(sig_go_id) > 0)
{
  sig_go_term <- getGOTerm(sig_go_id)#[["BP"]]
  GO_Table <- data.frame(GO.ID = sig_go_id, Adj.P.Val = go_fdr[go_fdr < 0.01], GO.Term = unlist(sig_go_term))
  length(GO_Table$GO.ID)
  write.table(GO_Table, paste(comparisons[i], "_P0.01_GO.txt", sep=""), sep="\t")
  amigo_tree = getAmigoTree(sig_go_id, pvalues=go_fdr[go_fdr < 0.01], webserver="http://amigo1.geneontology.org/cgi-bin/amigo/visualize",
                          filename=paste(comparisons[i], "_P0.001_upregulated_RamiGo.png", sep=""))
}
##'----------------------------------------------------------------------------------------------------------#
