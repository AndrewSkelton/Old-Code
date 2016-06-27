#-------------------------------------------------------------------------------------------#

library(stringr)
library(sva)
library(lumi)
library(gplots)
library(ggplot2)
library(annotate)
library(lumiHumanAll.db)
library(limma)
library(pathview)
library(GOstats)
library(RamiGO)

library(microRNA)
library(targetscan.Hs.eg.db)
library(MmPalateMiRNA)

##LOAD_DATA_&_PREPARE----------------------------------------------------------------------------------------#
setwd("/Users/andrew/Documents/Bioinformatics/Customers/Young/SR Microarray/Steven & Ruddy/Raw_data")
list.files()
filename <- "Raw_Data.txt"

raw_data <- lumiR(filename)
array_names <- c("S_mic_a",             "S_mic_c",               "S_3p_b",  "S_iso1_a", "S_iso1_c",
                 "S_iso2_b",            "S_HPC_a",               "S_HPC_c", "S_HP140_b","R_DO_L1_Line_1_D0_13",
                 "R_D0_L6_line_6_D0_21","R_D21_L3_Line_3_D21_13","S_mic_b", "S_3p_a",   "S_3p_c",
                 "S_iso1_b",            "S_iso2_a",              "S_iso2_c","S_HPC_b",  "S_HP140_a",
                 "S_HP140_c",           "R_D0_L3_Line_3_D0_13",  "R_D21_L1_Line_1_D21_13",
                 "R_D21_L6_line_6_D21_21")

sampleNames(raw_data) <- array_names

raw_data_in <- raw_data[,grep("_mic_|_3p_|_HPC_|_HP140_", colnames(raw_data))]
raw_data_in <- raw_data_in[, -c(3)]
raw_data_in <- raw_data_in[, -c(3,10)]
raw_data_in <- raw_data_in[, -c(1,3,10)]
sampleNames(raw_data_in)

png("boxplot_detection.png",width=12.53,height=6.98,units="in",res=600)
  boxplot(detection(raw_data_in),xaxt="n", cex=.2)
dev.off()

png("density_raw.png",width=12.53,height=6.98,units="in",res=600)
  plot(raw_data_in, main="")
dev.off()

#Normalisation - VST RSN
vst_data         <- lumiT(raw_data_in)
rsn_data         <- lumiN(vst_data, method = "rsn")
normalised_data  <- lumiQ(rsn_data)

png("density_raw.png",width=12.53,height=6.98,units="in",res=600)
  plot(normalised_data)
dev.off()

##'DETECTION_CALL_THRESHOLD_REMOVAL--------------------------------------------------------------------------#
exprs_data                   <- exprs(normalised_data)
present_count                <- detectionCall(normalised_data)
normalised_data              <- exprs_data[present_count > 0, ]
##'----------------------------------------------------------------------------------------------------------#

Treatment <- c('mic', 'mic', 'threep', 'HPC', 'HPC', 'HP140', 'mic', 'threep', 'threep', 'HPC', 'HP140', 'HP140')
Treatment <- c('mic', 'mic', 'HPC', 'HPC', 'HP140', 'mic', 'threep', 'threep', 'HPC', 'HP140', 'HP140')
Treatment <- c('mic', 'mic', 'HPC', 'HPC', 'HP140', 'mic', 'threep', 'threep', 'HP140', 'HP140')
#Treatment <- c('mic', 'HPC', 'HPC', 'HP140', 'mic', 'threep', 'threep', 'HP140', 'HP140')

png("PCA_wo3|10.png",width=12.53,height=6.98,units="in",res=600) 
    plotPCA(normalised_data, groups=as.numeric(as.factor(Treatment)), 
            groupnames=levels(as.factor(Treatment)), addtext=sampleNames(raw_data_in), pcs=c(1,2), 
            plot3d=F, main="")
dev.off()

pp <- plot_profile(normalised_data, Treatment, T)
png("pp_wo3|10.png",width=12.53,height=6.98,units="in",res=600)
  print(pp)
dev.off()

library(ggdendro)
hc = hclust(dist(t(normalised_data)))
g <- ggdendrogram(hc, theme_dendro = FALSE)
png("dendro_wo3|10.png",width=12.53,height=6.98,units="in",res=600)
  print(g)
dev.off()

png("boxplot_raw.png",width=12.53,height=6.98,units="in",res=600)
  boxplot(normalised_data)
dev.off()
##'----------------------------------------------------------------------------------------------------------#

##PROBE_ANNOTATION-------------------------------------------------------------------------------------------#
probe_list       <- rownames(normalised_data)
nuIDs            <- probeID2nuID(probe_list)[, "nuID"]
symbol           <- getSYMBOL(nuIDs, "lumiHumanAll.db")
name             <- unlist(lookUp(nuIDs, "lumiHumanAll.db", "GENENAME"))
anno_df          <- data.frame(ID = nuIDs, probe_list, symbol, name)
##'----------------------------------------------------------------------------------------------------------#

##EXPERIMENTAL_DESIGN----------------------------------------------------------------------------------------#
treatments  <- unique(Treatment)
array_names <- Treatment

design              <- model.matrix(~0 + factor(array_names, levels = treatments))
colnames(design)    <- treatments
num_parameters      <- ncol(design)

fit                 <- lmFit(normalised_data, design)
# fit                 <- lmFit(normalised_data[grep(paste(ID_Subset, collapse="|"), rownames(normalised_data)),], design)

cont_mat         <- makeContrasts(mic-threep, HPC-HP140, levels=treatments)
fit2             <- contrasts.fit(fit, contrasts=cont_mat)
fit2             <- eBayes(fit2)
fit2$genes = anno_df
##'----------------------------------------------------------------------------------------------------------#

##'DIFFERENTIAL EXPRESSION-----------------------------------------------------------------------------------#

comparisons <- c("mic - threep", "HPC - HP140")
p_cut_off   <- 0.01
fold_change <- 1.5
i           <- 1
mtc         <- 'none'

gene_list_1            <- topTable(fit2, coef = comparisons[1], p.value=p_cut_off, lfc=log2(fold_change),
                                 nrow(fit2$p.value), adjust.method=mtc)
length(gene_list_1$ID)
gene_list_unfiltered_1 <- topTable(fit2, coef=comparisons[1], number=nrow(anno_df), adjust.method=mtc)
gg_volcano(gene_list_unfiltered_1, fold_change, p_cut_off)
# ggv_volcano(gene_list_unfiltered_1, p_cut_off, fold_change)

gene_list_2            <- topTable(fit2, coef = comparisons[2], p.value = p_cut_off, lfc = log2(fold_change),
                                   nrow(fit2$p.value), adjust.method=mtc)
length(gene_list_2$ID)
gene_list_unfiltered_2 <- topTable(fit2, coef=comparisons[2], number=nrow(anno_df), adjust.method=mtc)
gg_volcano(gene_list_unfiltered_2, fold_change, p_cut_off)
# ggv_volcano(gene_list_unfiltered_2, p_cut_off, fold_change)

##'----------------------------------------------------------------------------------------------------------#

##'MIRNA TARGETING-------------------------------------------------------------------------------------------#
entrez_map_1 <- na.omit(nuID2EntrezID(gene_list_1$ID, lib.mapping='lumiHumanIDMapping', 
                                    filterTh = c(Strength1 = 95, Uniqueness = 95)))
entrez_map_1 <- entrez_map_1[entrez_map_1 != ""]
entrez_map_1 <- entrez_map_1[as.vector(unlist(lapply(entrez_map_1, 
                                                     function(x) exists(x, targetscan.Hs.egTARGETS)))) == T]
gene_targets_1 <- as.vector(unlist(mget(entrez_map_1, targetscan.Hs.egTARGETS)))

##'
##'
entrez_map_2 <- na.omit(nuID2EntrezID(gene_list_2$ID, lib.mapping='lumiHumanIDMapping', 
                                      filterTh = c(Strength1 = 95, Uniqueness = 95)))
entrez_map_2 <- entrez_map_2[entrez_map_2 != ""]
entrez_map_2 <- entrez_map_2[as.vector(unlist(lapply(entrez_map_2, 
                                                     function(x) exists(x, targetscan.Hs.egTARGETS)))) == T]
gene_targets_2 <- as.vector(unlist(mget(entrez_map_2, targetscan.Hs.egTARGETS)))


length(intersect(gene_targets_1, gene_targets_2))
length(as.vector(unlist(gene_targets_1)))
length(as.vector(unlist(gene_targets_2)))



