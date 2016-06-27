source("http://bioconductor.org/biocLite.R")
biocLite()
biocLite("lumi")
library(lumi)
install.packages("gplots")

###RAW DATASET
raw_data            <- lumiR("Raw_data.txt")
vst_data            <- lumiT(raw_data, method="vst")
rsn_data            <- lumiN(vst_data, method = "rsn")
analysis_ready_data <- lumiQ(rsn_data)
quan_data           <- lumiN(vst_data, method = "quantile")

pdf("Raw_Dataset.pdf")
  plot(raw_data,  main="Density Plot of Raw Dataset")
  boxplot(raw_data,  main="Boxplot of Raw Dataset")
  plot(vst_data,  main="Density Plot of Raw Dataset - Post VST")
  boxplot(vst_data,  main="Boxplot of Raw Dataset - Post VST")
  plot(rsn_data,  main="Density Plot of Raw Dataset - Post RSN")
  boxplot(rsn_data,  main="Boxplot of Raw Dataset - Post RSN")
  plot(quan_data, main="Density Plot of Raw Dataset - Post Quantile Analysis - No RSN")
  boxplot(quan_data,  main="Boxplot of Raw Dataset - Post Quantile Analysis - No RSN")
  plot(raw_data[, c(1, 2, 3, 4)], what="MAplot", main=paste(sampleNames(analysis_ready_data)[1:4], collapse=", "), cex=0.5)
  plot(analysis_ready_data[, c(31:36)], cex=0.9, what="MAplot", main=paste(sampleNames(analysis_ready_data)[31:36], collapse=", "))
  plot(analysis_ready_data, what="sampleRelation", method="mds")
dev.off()
###END OF RAW DATASET

#ARRAYS TO BE ELIMINATED
removal_elements <- c(1,2,4,19,24,31:36)
removal          <- sampleNames(analysis_ready_data)[removal_elements]
raw_data_post    <- raw_data[, !(sampleNames(raw_data) %in% removal)]
sampleNames(raw_data_post)
#END OF ARRAYS TO BE ELIMINATED

###OPTIMAL DATASET
vst_data_post            <- lumiT(raw_data_post, method="vst")
rsn_data_post            <- lumiN(vst_data_post, method="rsn")
analysis_ready_data_post <- lumiQ(rsn_data_post)

pdf("Optimal_Dataset_Final.pdf")
  plot(raw_data_post,  main="Density Plot of Optimal Dataset")
  boxplot(raw_data_post,  main="Boxplot of Optimal Dataset")
  plot(vst_data_post,  main="Density Plot of Optimal Dataset - Post VST")
  boxplot(vst_data_post,  main="Boxplot of Optimal Dataset - Post VST")
  plot(rsn_data_post,  main="Density Plot of Optimal Dataset - Post RSN")
  boxplot(rsn_data_post,  main="Boxplot of Optimal Dataset - Post RSN")
  plot(analysis_ready_data_post, main="Density Plot of Optimal Dataset - Post QC")
  boxplot(analysis_ready_data_post,  main="Boxplot of Optimal Dataset - Post QC")
  plot(raw_data_post[, c(1, 2, 3, 4)], what="MAplot", main=paste(sampleNames(analysis_ready_data_post)[1:4], collapse=", "))
  plot(analysis_ready_data_post[, c(26:30)], cex=0.9, what="MAplot", main=paste(sampleNames(analysis_ready_data_post)[26:30], collapse=", "))
  plot(analysis_ready_data_post, what="sampleRelation", method="mds")
dev.off()
###END OF OPTIMAL DATASET

###ARRAY QUALITY METRICS
biocLite("arrayQualityMetrics")
library("arrayQualityMetrics")

arrayQualityMetrics(expressionset=analysis_ready_data_post, outdir="QC_Optimal_2")
#arrayQualityMetrics(expressionset=analysis_ready_data_post, outdir="QC_Optimal")
arrayQualityMetrics(expressionset=analysis_ready_data, outdir="QC_With_Samples")
###END OF ARRAY QUALITY METRICS

###DETECTION CALL THRESHOLD REMOVAL
exprs_data                   <- exprs(analysis_ready_data_post)
present_count                <- detectionCall(analysis_ready_data_post) 
filtered_analysis_ready_data <- exprs_data[present_count > 0, ]
###PERCENTAGE OF PROBES REMOVED [43.25394%]
(nrow(filtered_analysis_ready_data)/nrow(analysis_ready_data_post))*100
###END OF DETECTION CALL THRESHOLD REMOVAL

###ARRAY ANNOTATION
biocLite("lumiHumanAll.db")
biocLite("annotate")
library(lumiHumanAll.db)
library(annotate)
probe_list <- rownames(filtered_analysis_ready_data)
nuIDs      <- probeID2nuID(probe_list)[, "nuID"]
symbol     <- getSYMBOL(nuIDs, "lumiHumanAll.db")
name       <- unlist(lookUp(nuIDs, "lumiHumanAll.db", "GENENAME"))
anno_df <- data.frame(ID = nuIDs, probe_list, symbol, name)
###END OF ARRAY ANNOTATION




###DIFFERENTIAL EXPRESSION ANALYSIS 
biocLite("limma")
library(limma) 
treatments  <- c("NOF", "OA")
array_names <- c("NOF","NOF","NOF","OA","OA","OA","NOF","NOF","OA","OA","NOF",
                 "OA","OA","NOF","NOF","OA","NOF","NOF","NOF","OA","OA","OA",
                 "OA","NOF","NOF")

#array_names <- c("NOF","NOF","NOF","OA","NOF","NOF","OA","OA","OA","NOF","NOF",
#                 "OA","OA","NOF","OA","OA","NOF","NOF","NOF","OA","NOF","NOF",
#                 "NOF","NOF","OA","OA","OA","OA","NOF","NOF")

design <- model.matrix(~0 + factor(array_names, levels = treatments))
colnames(design) <- treatments
num_parameters <- ncol(design)
fit <- lmFit(analysis_ready_data_post, design)

cont_mat <- makeContrasts(NOF-OA, levels=treatments)
fit2 <- contrasts.fit(fit, contrasts=cont_mat)
fit2 <- eBayes(fit2)
#fit2$genes = anno_df

## Filter by fold change (1.5x) and p-value (0.05) cutoffs
## Adjusted using Benjamini-Hochberg False Discovery Rate
topTable(fit2, coef="NOF - OA", p.value=0.01, lfc=log2(2))
#topTable(fit2, coef="OA - NOF", p.value=0.05, lfc=log2(1.5))
#topTable(fit2, coef="OA + NOF", p.value=0.05, lfc=log2(1.5))

gene_list   = topTable(fit2, coef = "NOF - OA", number = nrow(anno_df))
#gene_list_1 = topTable(fit2, coef = "OA + NOF", number = nrow(anno_df))
#gene_list_2 = topTable(fit2, coef = "OA - NOF", number = nrow(anno_df))

pdf("Volcano_Plots.pdf")

plot(gene_list$logFC, -log10(gene_list$adj.P.Val),
     #col=1+(abs(gene_list$logFC) > 1 & gene_list$adj.P.Val < 0.01), cex=0.4, 
     col="red", cex=0.4, 
     main="Volcano Plot: NOF - OA")
     abline(v=c(-1,1), col="red");abline(h=1.25, col="red")

#plot(gene_list_2$logFC, -log10(gene_list_2$adj.P.Val),
#     col=1+(abs(gene_list_2$logFC) > 1 & gene_list_2$adj.P.Val < 0.05), cex=0.4, 
#     pch=2, main="Volcano Plot: OA - NOF")
#     abline(v=c(-1,1), col="red");abline(h=1.25, col="red")

###GGPLOT2
library(ggplot2)
gene_list$threshold = as.factor(abs(gene_list$logFC) > 1 & gene_list$P.Value < 0.01/nrow(gene_list))#/nrow(anno_df))
g = ggplot(data=gene_list, aes(x=gene_list$logFC, y=-log10(gene_list$adj.P.Val), colour=gene_list$threshold)) +
  geom_point(alpha=0.4, size=1.75) +
  theme(legend.position = "none") +
  xlim(c(-5, 5)) + ylim(c(0, 15)) +
  xlab("log2 fold change") + ylab("-log10 p-value")
g
###END GGPLOT2

dev.off()

results = classifyTestsP(fit2, p.value = 0.01, method = "fdr")
vennDiagram(results[, c("NOF - OA", "OA - NOF")])

gene_list_complete = topTable(fit2, coef = "NOF - OA", p.value=0.01, lfc=log2(2), number = nrow(anno_df))
###END OF DIFFERENTIAL EXPRESSION ANALYSIS
gene_list_complete[1:10,]
###FUNCTIONAL ANALYSIS

###DIFFERENTIAL EXPRESSION HEATMAP
fa_data <- exprs(analysis_ready_data_post)
dist_measure <- dist(t(exprs(analysis_ready_data_post)), method="euclidian")
clust_by_samples <- hclust(dist_measure)
plot(clust_by_samples)

kmeans_by_samples <- kmeans(t(exprs(analysis_ready_data_post)), centers=8)
print(kmeans_by_samples$cluster)

#EXTRACT NOF-OA DIFFERENTIALLY EXPRESSED GENE DATA
nof_oa <- topTable(fit2, coef="NOF - OA", p.value=0.01,
                    lfc=log2(2), number=nrow(anno_df))

#REPLACE <NA> FACTORS WITH NO_DATA STRINGS
nof_oa_2 <- nof_oa
nof_oa_2 <- data.frame(lapply(nof_oa_2, as.character), stringsAsFactors=FALSE)
nof_oa_2[is.na(nof_oa_2)] <- "No_Data"
rownames(nof_oa_2) <- rownames(nof_oa)
colnames(nof_oa_2) <- colnames(nof_oa)
nof_oa_2 <- transform(nof_oa_2, logFC = as.numeric(logFC), AveExpr = as.numeric(AveExpr),t = as.numeric(t), 
                      P.Value = as.numeric(P.Value), adj.P.Val = as.numeric(adj.P.Val), B = as.numeric(B))
nof_oa_2$P.Value <- format(nof_oa_2$P.Value, scientific = FALSE)
nof_oa_2$adj.P.Val <- format(nof_oa_2$adj.P.Val, scientific = FALSE)

sapply(nof_oa_2, mode)
nof_oa_2[1:10,]

diff_exp <- fa_data[rownames(nof_oa_2),1:length(colnames(fa_data))]
length(diff_exp)
library("gplots")
heatmap.2(diff_exp)
heatmap.2(diff_exp, col="redblue", density.info="none",
          trace="none", cexRow=0.2)
###END OF DIFFERENTIAL EXPRESSION HEATMAP

###GO ANNOTATION

#tf1_sig         <- topTable(fit2, coef="TF1 - Ctrl", p.value=0.05, lfc=log2(1.5), number=nrow(select_data))
sig_probes      <- as.character(nof_oa_2$probe_list)
foo <- unique(unlist(sig_probes))
foo <- foo[-1]
foo
foo_1 <-nuIDs[sig_probes]
foo_2 <- foo_1[!is.na(foo_1)]
entrez <- lookUp(foo_2,"lumiHumanAll.db", "ENTREZID")
#entrez          <- unique(unlist(lookUp(nuIDs[sig_probes],"lumiHumanAll.db", "ENTREZID")))
entrez          <- as.character(entrez[!is.na(entrez)])
entrez_universe <- unique(unlist(lookUp(nuIDs,"lumiHumanAll.db", "ENTREZID")))
entrez_universe <- as.character(entrez_universe[!is.na(entrez_universe)])

biocLite("GOstats")
biocLite("RamiGO")
library("GOstats")
library("RamiGO")

params <- new("GOHyperGParams", geneIds=entrez,universeGeneIds=entrez_universe,annotation="lumiHumanAll.db",ontology="BP",
        pvalueCutoff= 0.01,conditional=FALSE,testDirection="over")
hyperg_result <- hyperGTest(params)
print(hyperg_result)

pval_go <- pvalues(hyperg_result)
go_fdr <- p.adjust(pvalues(hyperg_result), method="fdr")
sig_go_id <- names(go_fdr[go_fdr < 0.01])
sig_go_term <- getGOTerm(sig_go_id)[["BP"]]

top_5_go = sig_go_id[1:7]
top_5_p = go_fdr[go_fdr < 0.01][1:7]
amigo_tree = getAmigoTree(top_5_go, pvalues=top_5_p, filename='RamiGo_Final_Whole_7.png')

###END OF GO ANNOTATION

###KEGG
biocLite("pathview")
library('pathview')

kegg_params <- new("KEGGHyperGParams",
                     geneIds=entrez,
                     universeGeneIds=entrez_universe,
                     annotation="lumiHumanAll.db",
                     pvalueCutoff= 0.05,
                     testDirection="over")
kegg_hyperg_result <- hyperGTest(kegg_params)
print(kegg_hyperg_result)

sig_logFCs <- nof_oa_2$logFC
fc_entrez <- unlist(lookUp(nuIDs[sig_probes],
                             "lumiHumanAll.db", "ENTREZID"))
names(sig_logFCs) <- fc_entrez
sig_logFCs <- sig_logFCs[!is.na(entrez)]

pv <- pathview(gene.data=sig_logFCs,
               pathway.id=sigCategories(kegg_hyperg_result)[1],
               species="hsa",
               limit=list(gene=c(round(min(sig_logFCs)),
                                 round(max(sig_logFCs))), cpd=1))

############################
sig_logFCs   <- nof_oa_2$logFC
sig_probes_2 <- sig_probes[sig_probes!="No_Data"]
fc_entrez    <- unlist(lookUp(nuIDs[sig_probes_2],
                             "lumiHumanAll.db", "ENTREZID"))
names(sig_logFCs) <- names(fc_entrez)
sig_logFCs <- sig_logFCs[!is.na(entrez)]

sig_logFCs <- sig_logFCs[!is.na(names(sig_logFCs))]

pv <- pathview(gene.data=sig_logFCs,
               pathway.id=sigCategories(kegg_hyperg_result)[3],
               species="hsa",
               limit=list(gene=c(round(min(sig_logFCs)),
                                 round(max(sig_logFCs))), cpd=1))

###END OF KEGG

###END OF FUNCTIONAL ANALYSIS

write.csv(filtered_analysis_ready_data, "NOF_OA_Microarray_data.csv")
mean_exp_data <- c()
mean_exp_data <- cbind(mean_exp_data, rowMeans(filtered_analysis_ready_data[,grep('NOF', array_names)]))
mean_exp_data <- cbind(mean_exp_data, rowMeans(filtered_analysis_ready_data[,grep('OA', array_names)]))
colnames(mean_exp_data) <- c('NOF', 'OA')
write.csv(mean_exp_data, "NOF_OA_Mean_Microarray_data.csv")

