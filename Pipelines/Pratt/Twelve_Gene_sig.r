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
raw_data <- lumiB(raw_data, method="bgAdjust", verbose=T)

pheno_table           <- read.table("pheno_table.txt", header=T, sep="\t", stringsAsFactors=F)
rownames(pheno_table) <- pheno_table$X
# sampleNames(raw_data) <- pheno_table$X
pData(raw_data)       <- pheno_table

det <- detection(raw_data)
png("Raw_Detection_PValues.png",width=12.53,height=6.98,units="in",res=600) 
boxplot((det), main="Distribution of Detection P Values", xaxt="n", cex=.2)
dev.off()

grep("830_NRA_9477874180_K_SB3_RA3_ND", sampleNames(raw_data), value=F)#, 
##'115/116 (879/948) - Technical Failure
##'146/134 (835rpt/835) - Outliers (possible sample contamination)
##'
##'148/150/151/152/153/154/155/156
##'(836rpt/837rpt/832rpt/838rpt/833rpt/839rpt/834rpt/840rpt)
##'Repeat Samples
##'
##'923/926 (195/199)
##'
##' 2110_RA_9477874149_D_SB3_RA3_MR
##' 830_NRA_9477874180_K_SB3_RA3_ND
raw_data_det <- raw_data[, -c(52,143, 115,116,146,134,148,150,151,152,153,154,155,156,195,199)]
# raw_data_det <- raw_data_det[, (pData(raw_data_det)$Pool_ID != "UA" & pData(raw_data_det)$Sex == "F")]
# raw_data_det <- raw_data_det[, (pData(raw_data_det)$Pool_ID != "UA" & pData(raw_data_det)$RNA.Batch == "2")]
raw_data_det <- raw_data_det[, (pData(raw_data_det)$Pool_ID == "UA" & pData(raw_data_det)$RNA.Batch == "2")]
# raw_data_det <- raw_data_det[, (pData(raw_data_det)$Sex == "M")]
##'----------------------------------------------------------------------------------------------------------#

##'NORMALISATION---------------------------------------------------------------------------------------------#
# raw_bck          <- lumiB(raw_data_det, method="bgAdjust", verbose=T)
vst_data         <- lumiT(raw_data_det, method='vst')
rsn_data         <- lumiN(vst_data, method = "rsn")
lumi.Q           <- lumiQ(rsn_data)
##'QUANTILE_NORMALISATION - COMMENT OUT RSN NORMALISATION LINE IF YOU'RE USING QUANTILE
#quantile_data       <- lumiN(vst_data, method="quantile")
#lumi.Q              <- lumiQ(quantile_data)
foo <- estimateLumiCV(lumi.Q, 'measurement', ifPlot=T)
plot(NULL, xlim=c(-4,10), ylim=c(0,1))
for(i in 1:nrow(foo)) {
  lines(density(log2(foo[i,])), col=i)
}
plot(density(foo))
plot(lumi.Q)
##'----------------------------------------------------------------------------------------------------------#

##'DETECTION_CALL_THRESHOLD_REMOVAL--------------------------------------------------------------------------#
exprs_data                   <- exprs(lumi.Q)
present_count                <- detectionCall(lumi.Q)
normalised_data              <- exprs_data[present_count > 0, ]
##'----------------------------------------------------------------------------------------------------------#

##'BATCH CORRECTION (Amplification)--------------------------------------------------------------------------#
batches <- pData(lumi.Q)$RNA.Batch
# batches[grep('RA1', batches)] <- 1
# batches[grep('RA2', batches)] <- 2
# batches[grep('RA3', batches)] <- 3
# batches[grep('RA4', batches)] <- 4

pheno            <- data.frame(sample=c(1:ncol(normalised_data)),
                               outcome=pData(lumi.Q)$Pool_ID, batch=batches)
rownames(pheno)  <- colnames(normalised_data)
batch            <- pheno$batch
mod              <- model.matrix(~as.factor(outcome), data = pheno)
batchCorrected_data = ComBat(dat=normalised_data, batch=batch, mod=mod, numCovs=NULL,
                             par.prior=T, prior.plots=F)
##'----------------------------------------------------------------------------------------------------------#

##PROBE_ANNOTATION-------------------------------------------------------------------------------------------#
probe_list       <- rownames(batchCorrected_data)
nuIDs            <- probeID2nuID(probe_list)[, "nuID"]
symbol           <- getSYMBOL(nuIDs, "lumiHumanAll.db")
name             <- unlist(lookUp(nuIDs, "lumiHumanAll.db", "GENENAME"))
anno_df          <- data.frame(ID = nuIDs, probe_list, symbol, name)
##'----------------------------------------------------------------------------------------------------------#

##EXPERIMENTAL_DESIGN----------------------------------------------------------------------------------------#
treatments  <- c("RA", "NRA", "UA")
array_names <- pData(lumi.Q)$Pool_ID

design              <- model.matrix(~0 + factor(array_names, levels = treatments))
colnames(design)    <- treatments
num_parameters      <- ncol(design)
fit                 <- lmFit(batchCorrected_data, design)

cont_mat         <- makeContrasts(NRA-RA, levels=treatments)
fit2             <- contrasts.fit(fit, contrasts=cont_mat)
fit2             <- eBayes(fit2)
fit2$genes = anno_df
##'----------------------------------------------------------------------------------------------------------#

##'GENE SIGNATURE EXTRACT------------------------------------------------------------------------------------#
##'-------------------------------------NUID MAPPING
wg6v3_ID            <- c(6330725,4230102,3130301,3190609,6280170,6370082,
                         4730523,0670576,7650026,3190092,5550066,6770603)
wg6v3_nuid          <- probeID2nuID(wg6v3_ID)[, "nuID"]
ht12_wg6v3_nuid_map <- anno_df[anno_df$ID %in% wg6v3_nuid,]
##'-------------------------------------

##'-------------------------------------SYMBOL MAPPING
wg6v3_sig  <- c("BCL3","SOCS3","PIM1","SBNO2","PDCD1","GPRIN3",
                "IGFL2","LOC731186","MUC1","LDHA","CMAHP","NOG")
symbol_map <- anno_df[anno_df$symbol %in% wg6v3_sig,]
wg6v3_symb <- getSYMBOL(wg6v3_nuid, "lumiHumanAll.db")
##'-------------------------------------

##'-------------------------------------CONSENSUS SIGNATURE
selection        <- rbind(ht12_wg6v3_nuid_map, symbol_map)
selection_unique <- selection[!duplicated(selection[c("ID","probe_list")]),]
selection_unique <- selection_unique[order(rownames(selection_unique)), ]
selection_IDs    <- rownames(selection_unique)
##'-------------------------------------

##'-------------------------------------SUBSETTING
normalised_subset_RA <- (batchCorrected_data[grep(paste(selection_IDs, collapse="|"), rownames(batchCorrected_data)), grep("_RA_", colnames(batchCorrected_data))]) #
normalised_subset_RA <- normalised_subset_RA[order(rownames(normalised_subset_RA)), ]
# write.csv(normalised_subset, file="log2_all.csv")

normalised_subset_NRA <- (batchCorrected_data[grep(paste(selection_IDs, collapse="|"), rownames(batchCorrected_data)), grep("_NRA_", colnames(batchCorrected_data))]) #
normalised_subset_NRA <- normalised_subset_NRA[order(rownames(normalised_subset_NRA)), ]

normalised_subset <- 2^normalised_subset
##'-------------------------------------

##'-------------------------------------PROBE SEQ EXTRACT
foo_ID <- c()
foo_Seq <- c()
for(i in 1:length(selection_unique$probe_list)) {
  foo_ID <- c(foo_ID , (as.vector(selection_unique$probe_list)[i]))
  foo_Seq <- c(foo_Seq, id2seq(as.vector(selection_unique$ID)[i]))
}
foo <- cbind(foo_ID, foo_Seq)
##'-------------------------------------
##'----------------------------------------------------------------------------------------------------------#
fooB <- apply(normalised_subset, 1, median_zero)
median_zero <- function(x){
  z <- median(normalised_subset[1,])
  x <- x/z
  return(x)
}
##'----------------------------------------------------------------------------------------------------------#
##'
library(gplots)
tmp <- normalised_subset[, grep('_RA', colnames(normalised_subset))]
rownames(tmp) <- selection_unique$symbol
png("heatmap_12Gene_RA_and_NRA_median_centred.png",width=12.53,height=6.98,units="in",res=600)

#breaks for the core of the distribution
breaks=seq(-6, 6, by=0.2) #41 values
#now add outliers
breaks=append(breaks, 10)
breaks=append(breaks, -10, 0)
#create colour panel with length(breaks)-1 colours
mycol <- colorpanel(n=length(breaks)-1,low="green",mid="black",high="red")
  
  z <- t(scale(t(normalised_subset)))
  heatmap.2((z), scale="row", cexRow=0.5,  key=TRUE, col=redgreen(75),
            symkey=FALSE, density.info="none", trace="none", cexCol=.4, dendrogram=c("column"))
dev.off()

png("pp_12Gene.png",width=12.53,height=6.98,units="in",res=600) 
print(plot_profile(log2(normalised_subset), Treatment, T))
dev.off()
##'----------------------------------------------------------------------------------------------------------#

##PROBE_ANNOTATION-------------------------------------------------------------------------------------------#
probe_list       <- rownames(normalised_subset)
nuIDs            <- probeID2nuID(probe_list)[, "nuID"]
symbol           <- getSYMBOL(nuIDs, "lumiHumanAll.db")
name             <- unlist(lookUp(nuIDs, "lumiHumanAll.db", "GENENAME"))
anno_df_ss       <- data.frame(ID = nuIDs, probe_list, symbol, name)
##'----------------------------------------------------------------------------------------------------------#

##EXPERIMENTAL_DESIGN----------------------------------------------------------------------------------------#
treatments  <- c("RA", "NRA", "UA")
array_names <- Treatment

design              <- model.matrix(~0 + factor(array_names, levels = treatments))
colnames(design)    <- treatments
num_parameters      <- ncol(design)
fit_ss              <- lmFit(normalised_subset, design)

cont_mat            <- makeContrasts(NRA-RA, levels=treatments)
fit2_ss             <- contrasts.fit(fit_ss, contrasts=cont_mat)
fit2_ss             <- eBayes(fit2_ss)
fit2_ss$genes = anno_df_ss
##'----------------------------------------------------------------------------------------------------------#

##GENE_LISTS-------------------------------------------------------------------------------------------------#
comparisons <- c("RA - NRA", "RA - UA", "NRA - UA")
comparisons <- c("NRA - RA")
p_cut_off   <- 0.05
fold_change <- 1.2
i           <- 1
mtc         <- 'none'

gene_list_unfiltered <- topTable(fit2, coef=comparisons[i], number=nrow(anno_df), adjust.method=mtc)
gene_list_filtered   <- topTable(fit2, coef=comparisons[i], p.value=p_cut_off, 
                                 number=nrow(anno_df), lfc=log2(fold_change), adjust.method=mtc)
gg_volcano(gene_list_unfiltered, afc=1.2, pval=0.05)

ggv_volcano(gene_list_unfiltered, pval=0.05)






gene_list_unfiltered <- topTable(fit2, coef = comparisons[i], number=nrow(anno_df), adjust.method=mtc)
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