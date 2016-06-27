#-------------------------------------------------------------------------------------------#
#  Author      : Andrew J Skelton                                                           |
#  Language    : R Statistical Programming Language                                         |
#  Optimised   :                                                                            |
#  Data Owner  : Newcastle University - Ruddy                                               |
#  Description : A Microarray analysis pipeline - Mesenchymal stem cells vs Osteoblasts     |
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

##LOAD_DATA_&_PREPARE----------------------------------------------------------------------------------------#
setwd("/Users/andrew/Documents/Bioinformatics/Customers/Steven & Ruddy/Raw_data")
list.files()
filename <- "Raw_Data.txt"

raw_data <- lumiR(filename)
paste(sampleNames(raw_data), collapse=" ")

array_names <- c("S_mic_a",             "S_mic_c",               "S_3p_b",  "S_iso1_a", "S_iso1_c",
                 "S_iso2_b",            "S_HPC_a",               "S_HPC_c", "S_HP140_b","R_DO_L1_Line_1_D0_13",
                 "R_D0_L6_line_6_D0_21","R_D21_L3_Line_3_D21_13","S_mic_b", "S_3p_a",   "S_3p_c",
                 "S_iso1_b",            "S_iso2_a",              "S_iso2_c","S_HPC_b",  "S_HP140_a",
                 "S_HP140_c",           "R_D0_L3_Line_3_D0_13",  "R_D21_L1_Line_1_D21_13",
                 "R_D21_L6_line_6_D21_21")

sampleNames(raw_data) <- array_names
raw_data <- raw_data[, -grep('S_', sampleNames(raw_data))]

#Normalisation - VST RSN
vst_data         <- lumiT(raw_data)
rsn_data         <- lumiN(vst_data, method = "rsn")
normalised_data  <- lumiQ(rsn_data)

#QC Plots
boxplot(normalised_data)
plot(normalised_data)

pca_sample_names <- c("1_D0","6_D0","3_D21","3_D0","1_D21","6_D21")
pca_data <- (normalised_data)
colnames(pca_data) <- pca_sample_names

library(ggdendro)
hc = hclust(dist(t(exprs(normalised_data))))
g <- ggdendrogram(hc, theme_dendro = FALSE)

pdf("Raw_Data_QC.pdf")
boxplot(normalised_data)
plot(normalised_data)
plot(pca_data, cex=.4,what="sampleRelation", method="mds")
print(g)
dev.off()

###ARRAY QUALITY METRICS
biocLite("arrayQualityMetrics")
library("arrayQualityMetrics")
arrayQualityMetrics(expressionset=normalised_data, outdir="QC_1")
###END OF ARRAY QUALITY METRICS

##'DETECTION_CALL_THRESHOLD_REMOVAL--------------------------------------------------------------------------#
exprs_data                   <- exprs(normalised_data)
present_count                <- detectionCall(normalised_data)
normalised_data              <- exprs_data[present_count > 0, ]
##'----------------------------------------------------------------------------------------------------------#

##PROBE_ANNOTATION-------------------------------------------------------------------------------------------#
probe_list       <- rownames(normalised_data)
nuIDs            <- probeID2nuID(probe_list)[, "nuID"]
symbol           <- getSYMBOL(nuIDs, "lumiHumanAll.db")
name             <- unlist(lookUp(nuIDs, "lumiHumanAll.db", "GENENAME"))
anno_df          <- data.frame(ID = nuIDs, probe_list, symbol, name)
##'----------------------------------------------------------------------------------------------------------#

##EXPERIMENTAL_DESIGN----------------------------------------------------------------------------------------#
pairs  <- c(1,6,3,3,1,6)
treatments <- c("M", "O")
treatment  <- c("D0","D0","D21","D0","D21","D21")
treatment[grep('D0', treatment)]  <- "M"
treatment[grep('D21', treatment)] <- "O"

exp_design <- data.frame(Names=colnames(normalised_data), Pairs=pairs, Treat=treatment)
Treat  <- factor(exp_design$Treat, levels=c("M","O"))
Pairs  <- factor(exp_design$Pairs)
design <- model.matrix(~Pairs+Treat)
design

fit                 <- lmFit(normalised_data, design)
fit2                <- eBayes(fit)
fit2$genes = anno_df

gene_list_unfiltered <- topTable(fit2, coef="TreatO", number = nrow(anno_df), adjust.method="BH")
gene_list            <- topTable(fit2, coef="TreatO", p.value=0.05, lfc = log2(1.2),
                                 number = nrow(anno_df), adjust.method = "BH")

controls      <- normalised_data[as.character(gene_list$probe_list),
                                 grep("M", treatment)]
conditions    <- normalised_data[as.character(gene_list$probe_list),
                                 grep("O", treatment)]

gene_list$ControlMean   <- apply(controls,   1, mean)
gene_list$ConditionMean <- apply(conditions, 1, mean)

write.csv(gene_list, paste("gene_list_Ruddy","-", "0.05P_1.2FC.csv",
                           sep="",collapse=""))
##'-------------------

comparisons <- c("TreatO")
p_cut_off   <- 0.05
fold_change <- 1.2
i           <- 1

pdf("volcano_plots.pdf")
#for(i in 1:length(comparisons))#
#{#

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
dev.off()#
##'----------------------------------------------------------------------------------------------------------#

##'PATHWAY_ANALYSIS----------------------------------------------------------------------------------------------#

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
                            filename=paste(comparisons[i], "_P0.001_RamiGo.png", sep=""))
}
##'----------------------------------------------------------------------------------------------------------#


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

write.table(summary(kegg_hyperg_result), "Pathway_Summary.txt")

for(j in 1:length(sigCategories(kegg_hyperg_result)))
{
  pv <- pathview(gene.data=sig_logFCs,
                 pathway.id=sigCategories(kegg_hyperg_result)[j],
                 species="hsa",
                 limit=list(gene=c(round(min(sig_logFCs)),
                                   round(max(sig_logFCs))), cpd=1))
}
##'----------------------------------------------------------------------------------------------------------#
