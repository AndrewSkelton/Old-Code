list.files()

##'RAW DATASET NORMALISATION
raw_data            <- lumiR("Raw_data.txt")

#ARRAYS TO BE ELIMINATED
removal_elements <- c(1,2,4,19,24,31:36)
removal          <- sampleNames(raw_data)[removal_elements]
raw_data         <- raw_data[, !(sampleNames(raw_data) %in% removal)]
sampleNames(raw_data)
#END OF ARRAYS TO BE ELIMINATED

vst_data            <- lumiT(raw_data, method = "vst")
rsn_data            <- lumiN(vst_data, method = "rsn")
analysis_ready_data <- lumiQ(rsn_data)
arrayQualityMetrics(expressionset=analysis_ready_data, outdir="QC_INITIAL_NORMALISED_DATASET")


pdf("Raw_Dataset.pdf")
  plot(raw_data,  main="Density Plot of Raw Dataset")
  boxplot(raw_data,  main="Boxplot of Raw Dataset")
  plot(vst_data,  main="Density Plot of Raw Dataset - Post VST")
  boxplot(vst_data,  main="Boxplot of Raw Dataset - Post VST")
  plot(rsn_data,  main="Density Plot of Raw Dataset - Post RSN")
  boxplot(rsn_data,  main="Boxplot of Raw Dataset - Post RSN")
  #plot(quan_data, main="Density Plot of Raw Dataset - Post Quantile Analysis - No RSN")
  #boxplot(quan_data,  main="Boxplot of Raw Dataset - Post Quantile Analysis - No RSN")
  plot(raw_data[, c(1, 2, 3, 4)], what="MAplot", main=paste(sampleNames(analysis_ready_data)[1:4], collapse=", "), cex=0.5)
  plot(analysis_ready_data[, c(31:36)], cex=0.9, what="MAplot", main=paste(sampleNames(analysis_ready_data)[31:36], collapse=", "))
  plot(analysis_ready_data, what="sampleRelation", method="mds")
dev.off()

pdf("Outlier_Detection.pdf")
plotSampleRelation(analysis_ready_data, method= 'cluster' , cv.Th=0, main="Hierarchical Clustering", cex=.4)
plotSampleRelation(analysis_ready_data, method= 'mds' , cv.Th=0, main="PCA Plot")
dev.off()

###DETECTION CALL THRESHOLD REMOVAL
exprs_data                   <- exprs(analysis_ready_data)
present_count                <- detectionCall(analysis_ready_data) 
filtered_analysis_ready_data <- exprs_data[present_count > 0, ]
###PERCENTAGE OF PROBES REMOVED [40.519%]
(nrow(filtered_analysis_ready_data)/nrow(analysis_ready_data))*100
###END OF DETECTION CALL THRESHOLD REMOVAL

###EXPERIMENTAL SETUP
treatments  <- c("NOF", "OA")
array_names <- c("NOF", "NOF", "NOF", "OA", "OA", "OA", "NOF", "NOF", "OA", "OA", "NOF", "OA",
                 "OA", "NOF", "NOF", "OA", "NOF", "NOF", "NOF", "OA", "OA", "OA", "OA", "NOF", "NOF")

###ARRAY ANNOTATION
probe_list <- rownames(filtered_analysis_ready_data)#analysis_ready_data)
nuIDs      <- probeID2nuID(probe_list)[, "nuID"]
symbol     <- getSYMBOL(nuIDs, "lumiHumanAll.db")
name       <- unlist(lookUp(nuIDs, "lumiHumanAll.db", "GENENAME"))
anno_df    <- data.frame(ID = nuIDs, probe_list, symbol, name)
###END OF ARRAY ANNOTATION

design <- model.matrix(~0 + factor(array_names, levels = treatments))
colnames(design) <- treatments
num_parameters <- ncol(design)
fit <- lmFit(filtered_analysis_ready_data, design)

cont_mat <- makeContrasts(NOF-OA, levels=treatments)
fit2 <- contrasts.fit(fit, contrasts=cont_mat)
fit2 <- eBayes(fit2)

#NOF Vs OA
topTable(fit2, coef="NOF - OA", p.value=0.01, lfc=log2(2), adjust.method="BH")
gene_list_NOF_OA_unfiltered <- topTable(fit2, coef = "NOF - OA", number = nrow(anno_df), adjust.method="BH")
gene_list_NOF_OA            <- topTable(fit2, coef = "NOF - OA", p.value=0.01, lfc=log2(2), number = nrow(anno_df), 
                                                        adjust.method="BH")
controls      <- filtered_analysis_ready_data[as.character(gene_list_NOF_OA$probe_list), grep("NOF", array_names)]
conditions    <- filtered_analysis_ready_data[as.character(gene_list_NOF_OA$probe_list), grep("OA", array_names)]
gene_list_NOF_OA$ControlMean   <- apply(controls, 1, mean)
gene_list_NOF_OA$ConditionMean <- apply(conditions, 1, mean) 

corr <- cor(probe_level_data, method="pearson")

write.table(gene_list_NOF_OA,    "NOF_OA.txt",    sep="\t")
write.table(probe_level_data,    "NOF_OA_Probes.txt",    sep="\t")
write.table(corr,    "NOF_OA_corr.txt",    sep="\t")

##HIERARCHICAL CLUSTERING
pdf("bw_Microarray_pca.pdf")
fa_data          <- filtered_analysis_ready_data
dist_measure     <- dist(t(filtered_analysis_ready_data), method="euclidian")
plotSampleRelation(fa_data, method= 'mds' , cv.Th=0, main="PCA Plot of NOF and OA Microarray Samples")
colnames(fa_data) <- array_names
plotSampleRelation(fa_data, method= 'mds' , cv.Th=0, main="PCA Plot of NOF and OA Microarray Samples")
dev.off()


install.packages("ggdendro")
library(ggdendro)
pdf("bw_Microarray_Dendogram.pdf")
hc = hclust(dist(t(fa_data)))
hc1 = hclust(dist(t(filtered_analysis_ready_data)))
gg <- ggdendrogram(hc, rotate = TRUE, size = 4, theme_dendro = FALSE) +
labs(title = "Dendogram of NOF and OA Microarray Samples")
print(gg)
gg <- ggdendrogram(hc1, rotate = TRUE, size = 4, theme_dendro = FALSE) +
labs(title = "Dendogram of NOF and OA Microarray Samples")
print(gg)
dev.off()
plot(princomp(fa_data))

##VOLCANO PLOT
pdf("bw_Microarray_Volcano_Plot.pdf")
gene_list <- eval(parse(text="gene_list_NOF_OA_unfiltered"))  

gene_list$threshold = as.factor(abs(gene_list$logFC) > log2(2) & gene_list$adj.P.Val < 0.01)#/nrow(gene_list))#/nrow(anno_df))

g = ggplot(data=gene_list, aes(x=gene_list$logFC, y=-log10(gene_list$adj.P.Val), )) +#colour=gene_list$threshold)) +
  geom_point(alpha=0.4, size=1.75) +
  labs(title = "NOF - OA") + 
  theme(legend.position = "none") +
  scale_fill_manual(my_palette) +
  geom_vline(xintercept=log2(2)) + geom_vline(xintercept=-log2(2)) + geom_hline(yintercept=-log10(0.01)) +
  xlab("log2 fold change") + ylab("-log10 p-value")
print(g)
dev.off()

##HEATMAP BW
pdf("bw_Microarray_Heatmap.pdf")
par(cex.main=1)
diff_exp <- filtered_analysis_ready_data[rownames(gene_list_NOF_OA),]
colnames(diff_exp) <- array_names
my_palette <- colorRampPalette(c("white", "grey", "black"))(n = 1000)
heatmap.2(diff_exp, col=my_palette, density.info="none",
          trace="none", cexRow=0.2, cexCol=0.8,
          main="NOF Vs OA Differentially Expressed Genes")
dev.off()

##GO ANNOTATION
gene_list1      <- eval(parse(text="gene_list_NOF_OA"))  
sig_probes      <- as.character(gene_list1$probe_list)
entrez          <- unique(unlist(lookUp(nuIDs[sig_probes],
                                        "lumiHumanAll.db", "ENTREZID")))
entrez          <- as.character(entrez[!is.na(entrez)])
entrez_universe <- unique(unlist(lookUp(nuIDs,"lumiHumanAll.db", "ENTREZID")))
entrez_universe <- as.character(entrez_universe[!is.na(entrez_universe)])

params <- new("GOHyperGParams", geneIds=entrez,universeGeneIds=entrez_universe,annotation="lumiHumanAll.db",ontology="BP",
              pvalueCutoff= 0.001, conditional=FALSE, testDirection="over")
hyperg_result <- hyperGTest(params)
print(hyperg_result)
#str(summary(hyperg_result))

pval_go     <- pvalues(hyperg_result)
go_fdr      <- p.adjust(pvalues(hyperg_result), method="fdr")
sig_go_id   <- names(go_fdr[go_fdr < 0.001])
if(length(sig_go_id) > 0)
{
  sig_go_term <- getGOTerm(sig_go_id)#[["BP"]]
  GO_Table <- data.frame(GO.ID = sig_go_id, Adj.P.Val = go_fdr[go_fdr < 0.001], GO.Term = unlist(sig_go_term))
  length(GO_Table$GO.ID)
  write.table(GO_Table, paste("NOF_OA", "_P0.01_GO.txt", sep=""), sep="\t")
  amigo_tree = getAmigoTree(sig_go_id, pvalues=go_fdr[go_fdr < 0.001], filename=paste("NOF_OA", "_P0.001_upregulated_RamiGo.png", sep=""))
  #amigo_tree = getAmigoTree(top_5_go, pvalues=top_5_p, filename=paste(gene_lists_1[i], "__RamiGo.png", sep=" "))
}

#KEGG PATHWAYS
gene_list1      <- eval(parse(text="gene_list_NOF_OA"))

sig_probes      <- as.character(gene_list1$probe_list)
entrez          <- unique(unlist(lookUp(nuIDs[sig_probes],
                                        "lumiHumanAll.db", "ENTREZID")))
entrez          <- as.character(entrez[!is.na(entrez)])
entrez_universe <- unique(unlist(lookUp(nuIDs,"lumiHumanAll.db", "ENTREZID")))
entrez_universe <- as.character(entrez_universe[!is.na(entrez_universe)])

kegg_params <- new("KEGGHyperGParams",
                   geneIds=entrez,
                   universeGeneIds=entrez_universe,
                   annotation="lumiHumanAll.db",
                   pvalueCutoff= 0.001,
                   testDirection="over")
kegg_hyperg_result <- hyperGTest(kegg_params)
print(kegg_hyperg_result)

sig_logFCs <- gene_list1$logFC
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

corr <- cor(probe_level_data, method="pearson")

install.packages("qgraph")
install.packages("fdrtool")
library("qgraph")
library("fdrtool")
grouping <- list()
grouping$NOF <- c(1:13)
grouping$OA <- c(14:25)
q <- qgraph(corr,minimum=0.45,cut=1,vsize=4, groups = grouping,legend=TRUE,borders=FALSE)
#title("Big 5 correlations",line=-2,cex.main=2)
qgraph(q, gray = TRUE, layout = "spring")
qgraph(q, graph = "sig")
qgraph(corr,minimum=0.45,cut=1,vsize=8, legend=TRUE,borders=FALSE, graph = "sig")

qgraph(corr, minimum = 0.25, cut = 0.4, vsize = 1.5, 
       legend = TRUE, borders = FALSE, scores = as.integer(corr[sample(1:20, 1),]), scores.range = c(1, 5))

corr,minimum=0.45,cut=1,vsize=8, legend=TRUE,borders=FALSE


