#-------------------------------------------------------------------------------------------#
#  Author      : Andrew J Skelton                                                           |
#  Language    : R Statistical Programming Language                                         |
#  Optimised   : CD4 T-Cell RA Microarray Data                                              |
#  Data Owner  : Newcastle University - Dr. A. Pratt                                        |
#  Description : A Microarray analysis pipeline for downstream analysis of a normalised     |
#                CSV file                                                                   |
#-------------------------------------------------------------------------------------------#

setwd("/Users/andrew/Documents/Bioinformatics/Customers/Pratt/Rheumatoid_Arthritis_CD4_T-Cell_Signature/Raptor")
list.files()
normalised_dataset <- read.csv("normalised_R_readin.csv", strip.white=TRUE, header=T, stringsAsFactors=FALSE)
rownames(normalised_dataset) <- normalised_dataset[,1]
normalised_dataset <- normalised_dataset[,-1]
colnames(normalised_dataset) <- substring(colnames(normalised_dataset), 2)
normalised_dataset[normalised_dataset == "0.01*"] <- 0.0111
normalised_dataset <- data.matrix(normalised_dataset)
#present_count                <- detectionCall(normalised_dataset)
#adj_dataset           <- normalised_dataset[present_count > 0, ]
#(nrow(normalised_dataset)/nrow(normalised_dataset_raw))*100

##COHORT_SUBSET----------------------------------------------------------------------------------------------#
normalised_dataset_T <- normalised_dataset[,grep("_T_", colnames(normalised_dataset))]
##'----------------------------------------------------------------------------------------------------------#

##PROBE_ANNOTATION-------------------------------------------------------------------------------------------#
probe_list       <- rownames(normalised_dataset_T)
nuIDs            <- probeID2nuID(probe_list)[, "nuID"]
which(is.na(nuIDs)) #15740
nuIDs[15740]     <- probeID2nuID(3450064)[, "nuID"]
symbol           <- getSYMBOL(nuIDs, "lumiHumanAll.db")
name             <- unlist(lookUp(nuIDs, "lumiHumanAll.db", "GENENAME"))
anno_df          <- data.frame(ID = nuIDs, probe_list, symbol, name)

treatments  <- c("AB", "CD")
array_names <- colnames(normalised_dataset_T)
array_names[grep("_A_|_B_", colnames(normalised_dataset_T))] <- "AB"
array_names[grep("_C_|_D_", colnames(normalised_dataset_T))] <- "CD"

#normalised_dataset_UA2 <- data.matrix(normalised_dataset_UA)
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
topTable(fit2, coef="AB - CD", p.value=0.05, lfc=log2(1.2), adjust.method="none")
gene_list_Coef_A_unfiltered <- topTable(fit2, coef = "AB - CD", number = nrow(anno_df), adjust.method="none")
#gene_list_Coef_A <- gene_list_Coef_A_unfiltered[abs(gene_list_Coef_A_unfiltered$logFC) > log2(1.2) &
#                                                  gene_list_Coef_A_unfiltered$P.Value < 0.05,]
gene_list_Coef_A <- topTable(fit2, coef = "AB - CD", adjust.method="BH", p.value=0.05,
                             lfc=log2(1.2), number=nrow(anno_df))

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
  #xlim(c(-2.75, 2.75)) + ylim(c(0, 4.5)) +
  xlab("log2 fold change") + ylab("-log10 p-value") +
  ggtitle("Volcano Plot of Fold Change vs P-Value (No Multiple Test Correction)")
g
##'----------------------------------------------------------------------------------------------------------#


write.csv(gene_list_Coef_A, "gene_list_fully_Normalised.csv")
