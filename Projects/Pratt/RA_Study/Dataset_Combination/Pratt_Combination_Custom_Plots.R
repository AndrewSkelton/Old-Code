#-------------------------------------------------------------------------------------------#
#  Author      : Andrew J Skelton                                                           |
#  Language    : R Statistical Programming Language                                         |
#  Study       : Data Combination Project - Early Onset RA Study                            |
#  Data Owner  : Newcastle University - Dr. Arthur Pratt                                    |
#  Description : Custom Plots                                                               |
#-------------------------------------------------------------------------------------------#



##'Set the Working Directory, load essential packages, and import files
##'-----------------------------------------------------------------------------------------#
library(gplots)
library(ggplot2)
library(reshape2)
library(scales)
library(ComplexHeatmap)
library(circlize)
##'-----------------------------------------------------------------------------------------#


##'PCA - Confirmation Study - Normalised, No Technical Correction
##'-----------------------------------------------------------------------------------------#
#png("PCA_dsA_Norm_Uncorr.png", width=4096, height=3096, units="px", res=300)
pca          <- prcomp(t(normalised_data_a))
d            <- as.data.frame(pca$x)
d            <- cbind(d, pData(data_alumi.Q))
d$Diagnosis  <- d$Update_Collapse_2
d$Amplification_Batch <- as.factor(d$Amplification_Batch)

ggplot(d, aes(x=PC1, y=PC2, shape=Amplification_Batch)) +
  geom_point(aes(colour=Amplification_Batch), size=3) +
  # geom_text(label=d$Sample_ID, size=4, vjust=1.2, hjust=-0.2) +
  theme_bw() +
  ggtitle("PCA of Dataset 'A' - Normalised, Uncorrected" )
  theme(axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15))
#dev.off()
##'-----------------------------------------------------------------------------------------#


##'PCA - Confirmation Study - Normalised, Technical Correction Applied
##'-----------------------------------------------------------------------------------------#
#png("PCA_dsA_Norm_Corr.png", width=4096, height=3096, units="px", res=300)
pca          <- prcomp(t(exprs(data_alumi.Q)))
d            <- as.data.frame(pca$x)
d            <- cbind(d, pData(data_alumi.Q))
d$Diagnosis  <- d$Update_Collapse_2
d$Amplification_Batch <- as.factor(d$Amplification_Batch)

ggplot(d, aes(x=PC1, y=PC2, shape=Amplification_Batch)) +
  geom_point(aes(colour=Amplification_Batch), size=3) +
  # geom_text(label=d$Sample_ID, size=4, vjust=1.2, hjust=-0.2) +
  theme_bw() +
  ggtitle("PCA of Dataset 'A' - Normalised, Corrected" )
  theme(axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15))
#dev.off()
##'-----------------------------------------------------------------------------------------#


##'PCA - Pilot Study 1 - Normalised, No Technical Correction
##'-----------------------------------------------------------------------------------------#
#png("PCA_dsB_Norm_Uncorr.png", width=4096, height=3096, units="px", res=300)
pca          <- prcomp(t(normalised_data_b))
d            <- as.data.frame(pca$x)
d            <- cbind(d, pData(data_blumi.Q))
d$Diagnosis  <- d$Update_Collapse_2
d$Amplification_Batch <- as.factor(d$Amplification_Batch)

ggplot(d, aes(x=PC1, y=PC2, shape=Amplification_Batch)) +
  geom_point(aes(colour=Amplification_Batch), size=3) +
  # geom_text(label=d$Sample_ID, size=4, vjust=1.2, hjust=-0.2) +
  theme_bw() +
  ggtitle("PCA of Dataset 'B' - Normalised, Uncorrected" )
theme(axis.title.x = element_text(size=15),
      axis.title.y = element_text(size=15))
#dev.off()
##'-----------------------------------------------------------------------------------------#


##'PCA - Pilot Study 1 - Normalised, Technical Correction Applied
##'-----------------------------------------------------------------------------------------#
#png("PCA_dsB_Norm_Corr.png", width=4096, height=3096, units="px", res=300)
pca          <- prcomp(t(exprs(data_blumi.Q)))
d            <- as.data.frame(pca$x)
d            <- cbind(d, pData(data_blumi.Q))
d$Diagnosis  <- d$Update_Collapse_2
d$Amplification_Batch <- as.factor(d$Amplification_Batch)

ggplot(d, aes(x=PC1, y=PC2, shape=Amplification_Batch)) +
  geom_point(aes(colour=Amplification_Batch), size=3) +
  # geom_text(label=d$Sample_ID, size=4, vjust=1.2, hjust=-0.2) +
  theme_bw() +
  ggtitle("PCA of Dataset 'B' - Normalised, Corrected" )
  theme(axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15))
#dev.off()
##'-----------------------------------------------------------------------------------------#


##'PCA - Pilot Study 2 - Normalised, No Technical Correction
##'-----------------------------------------------------------------------------------------#
#png("PCA_dsC_Norm_Uncorr.png", width=4096, height=3096, units="px", res=300)
pca          <- prcomp(t(normalised_data_c))
d            <- as.data.frame(pca$x)
d            <- cbind(d, pData(data_clumi.Q))
d$Diagnosis  <- d$Update_Collapse_2
d$Amplification_Batch <- as.factor(d$Amplification_Batch)

ggplot(d, aes(x=PC1, y=PC2, shape=Amplification_Batch)) +
  geom_point(aes(colour=Amplification_Batch), size=3) +
  # geom_text(label=d$Sample_ID, size=4, vjust=1.2, hjust=-0.2) +
  theme_bw() +
  ggtitle("PCA of Dataset 'C' - Normalised, Uncorrected" )
  theme(axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15))
#dev.off()
##'-----------------------------------------------------------------------------------------#


##'PCA - Pilot Study 2 - Normalised, Technical Correction Applied
##'-----------------------------------------------------------------------------------------#
#png("PCA_dsV_Norm_Corr.png", width=4096, height=3096, units="px", res=300)
pca          <- prcomp(t(exprs(data_clumi.Q)))
d            <- as.data.frame(pca$x)
d            <- cbind(d, pData(data_clumi.Q))
d$Diagnosis  <- d$Update_Collapse_2
d$Amplification_Batch <- as.factor(d$Amplification_Batch)

ggplot(d, aes(x=PC1, y=PC2, shape=Amplification_Batch)) +
  geom_point(aes(colour=Amplification_Batch), size=3) +
  # geom_text(label=d$Sample_ID, size=4, vjust=1.2, hjust=-0.2) +
  theme_bw() +
  ggtitle("PCA of Dataset 'C' - Normalised, Corrected" )
  theme(axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15))
#dev.off()
##'-----------------------------------------------------------------------------------------#


##'PCA - Combination - No Correction
##'-----------------------------------------------------------------------------------------#
#png("PCA_Combined_No_Correction.png", width=4096, height=3096, units="px", res=300)
pca          <- prcomp(t(exprs(combined_none)))
d            <- as.data.frame(pca$x)
d            <- cbind(d, pData(combined_none))
d$Diagnosis  <- d$Update_Collapse_2


ggplot(d, aes(x=PC1, y=PC2, shape=Source)) +
  geom_point(aes(colour=Source), size=3) +
  # geom_text(label=d$Sample_ID, size=4, vjust=1.2, hjust=-0.2) +
  theme_bw() +
  ggtitle("PCA of 3 Datasets Combined Using No Correction")
  theme(axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15))
#dev.off()
##'-----------------------------------------------------------------------------------------#


##'PCA - Combination - Correction Applied
##'-----------------------------------------------------------------------------------------#
#png("PCA_Combined_Correction.png", width=4096, height=3096, units="px", res=300)
pca          <- prcomp(t(exprs(combined_combat)))
d            <- as.data.frame(pca$x)
d            <- cbind(d, pData(combined_combat))
d$Diagnosis  <- d$Update_Collapse_2

ggplot(d, aes(x=PC1, y=PC2, shape=Source)) +
  geom_point(aes(colour=Source), size=3) +
  # geom_text(label=d$Sample_ID, size=4, vjust=1.2, hjust=-0.2) +
  theme_bw() +
  ggtitle("PCA of 3 Datasets Combined Using ComBat")
  theme(axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15))
#dev.off()
##'-----------------------------------------------------------------------------------------#


##'Volcano - Combination
##'-----------------------------------------------------------------------------------------#
png("Volcano_Combination_DE_Included.png", width=4096, height=3096, units="px", res=300)
afc                <- 1.2
pval               <- 0.05
filtered_in        <- topTable(fit2, coef=1, p.value=pval, lfc=log2(afc),
                               number=Inf, adjust.method="BH")
unfiltered_in      <- topTable(fit2, coef=1, number=Inf,
                               adjust.method="BH")
unfiltered_in$pass <- unfiltered_in$adj.P.Val < pval & abs(unfiltered_in$logFC) > log2(afc)
unfiltered_in$lab  <- ""
unfiltered_in[match(filtered_in$ID, unfiltered_in$ID),]$lab <- as.vector(filtered_in$symbol)

signature_mapping  <- anno_df[match(signature, anno_df$ID),]
unfiltered_in[match(signature_mapping$ID, unfiltered_in$ID),]$lab <- as.vector(signature_mapping$symbol)
unfiltered_in$Signature <- "False"
unfiltered_in[match(signature_mapping$ID, unfiltered_in$ID),]$Signature <- "True"

ggplot(unfiltered_in, aes(x=logFC, y=-log(adj.P.Val, 10))) +
  geom_point(aes(colour=Signature, shape=pass), show_guide=F) +
  scale_colour_manual(values=c(alpha('grey', 0.5), 'red')) +
  geom_hline(yintercept=-log(pval, 10), colour="black", linetype=2) +
  geom_vline(xintercept=-log2(afc), colour="black", linetype=2) +
  geom_vline(xintercept=log2(afc), colour="black", linetype=2) +
  geom_text(label=unfiltered_in$lab, size=4, vjust=-1, hjust=0.5) +
  ggtitle("Combined Datasets Differential Expression: NRA-RA\nAdj.P < 0.05 | FC > log2(1.2) | BH Corrected") +
  theme_bw()
dev.off()
##'-----------------------------------------------------------------------------------------#


##'Heatmap - Signature
##'-----------------------------------------------------------------------------------------#
png("Heatmap_Signature.png", width=4096, height=3096, units="px", res=300)
anno_df$symbol <- as.vector(anno_df$symbol)
anno_df[match(signature[c(11,12)], anno_df$ID),]$symbol <- c("LOC731186", "GPRIN3")

#Useful
#nuID2IlluminaID(signature)

signature_exp                 <- exprs(raw_data_det)[rownames(raw_data_det) %in% signature,
                                                     grep("RA", pData(raw_data_det)$Update_Collapse_2)]
gene_names                    <- as.vector(anno_df[match(rownames(signature_exp),
                                                   anno_df$ID),]$symbol)
gene_names[is.na(gene_names)] <- "No Annotation"
rownames(signature_exp)       <- gene_names
foo                           <- as.vector(pData(raw_data_det)$Update_Collapse_2)
foo                           <- foo[grep("RA", foo)]
foo[foo == "RA"]              <- "black"
foo[foo == "NRA"]             <- "grey"
foo[foo == "UA"]              <- "white"
foo[foo == "OA"]              <- "grey"
par(cex.main=1.5)
heatmap.2(signature_exp,
          scale='row',
          dendrogram="column",
          trace='none',
          Rowv=NA,
          col=rev(redblue(100)),
          ColSideColors=foo,
          cexRow=1,
          labCol="",
          breaks=seq(-2,2,length.out=101),
          keysize = .9,
          density.info=c("none"),
          main="Data Combination - Signature Genes")
dev.off()
##'-----------------------------------------------------------------------------------------#


##'Heatmap - Signature - Key Players Only
##'-----------------------------------------------------------------------------------------#
png("Heatmap_Key_Players.png", width=4096, height=3096, units="px", res=300)
signature_exp                 <- signature_exp[match(c("BCL3", "SOCS3", "PIM1"),
                                               rownames(signature_exp)),]
par(cex.main=1.5)
heatmap.2(signature_exp,
          scale='row',
          dendrogram="column",
          trace='none',
          Rowv=NA,
          col=rev(redblue(100)),
          ColSideColors=foo,
          labCol="",
          breaks=seq(-2,2,length.out=101),
          keysize = .9,
          density.info=c("none"),
          cexRow=1,
          main="Key Players - BCL3, SOCS3, and PIM1")
dev.off()
##'-----------------------------------------------------------------------------------------#


##'Heatmap - Genes from Differential Expression
##'-----------------------------------------------------------------------------------------#
png("Heatmap_DE.png", width=4096, height=3096, units="px", res=300)
signature_exp                 <- exprs(raw_data_det)[rownames(raw_data_det) %in% filtered_in$ID,
                                                     grep("RA", pData(raw_data_det)$Update_Collapse_2)]
gene_names                    <- as.vector(anno_df[match(rownames(signature_exp),
                                                         anno_df$ID),]$symbol)
gene_names[is.na(gene_names)] <- "No Annotation"
rownames(signature_exp)       <- gene_names
foo                           <- as.vector(pData(raw_data_det)$Update_Collapse_2)
foo                           <- foo[grep("RA", foo)]
foo[foo == "RA"]              <- "black"
foo[foo == "NRA"]             <- "grey"
foo[foo == "UA"]              <- "white"
foo[foo == "OA"]              <- "grey"

par(cex.main=1.5)
heatmap.2(signature_exp,
          scale='row',
          dendrogram="column",
          trace='none',
          Rowv=NA,
          col=rev(redblue(100)),
          ColSideColors=foo,
          cexRow=1,
          labCol="",
          breaks=seq(-2,2,length.out=101),
          keysize = .9,
          density.info=c("none"),
          main="Combined Dataset DE Genes")
dev.off()
##'-----------------------------------------------------------------------------------------#



##'Utility - Extract Expression Data
##'-----------------------------------------------------------------------------------------#
anno_df$symbol <- as.vector(anno_df$symbol)
anno_df[match(signature[c(11,12)], anno_df$ID),]$symbol <- c("LOC731186", "GPRIN3")

signature_exp                 <- exprs(raw_data_det)[rownames(exprs(raw_data_det)) %in% signature,]
gene_names                    <- as.vector(anno_df[match(rownames(signature_exp),
                                                         anno_df$ID),]$symbol)
gene_names[is.na(gene_names)] <- "No Annotation"
rownames(signature_exp)       <- gene_names
# 
df_exp                        <- exprs(raw_data_det)[rownames(exprs(raw_data_det)) %in% filtered_in$ID,]
gene_names                    <- as.vector(anno_df[match(rownames(df_exp),
                                                         anno_df$ID),]$symbol)
gene_names[is.na(gene_names)] <- "No Annotation"
rownames(df_exp)              <- gene_names

write.csv(signature_exp, file="Data_Combination_Signature_Normalised_Expression.csv")
write.csv(df_exp, file="Data_Combination_DE_Normalised_Expression.csv")
write.csv(pData(raw_data_det), file="Data_Combination_Pheno_Data.csv")

##'-----------------------------------------------------------------------------------------#



##'Utility - Extract Expression Data
##'-----------------------------------------------------------------------------------------#
foo_mat <- signature_exp
colnames(foo_mat) <- foo
foo_ua            <- exprs(raw_data_det)[rownames(raw_data_det) %in% signature,
                                         grep("UA", pData(raw_data_det)$Update_Collapse_2)]
colnames(foo_ua)  <- rep("UA", ncol(foo_ua))
rownames(foo_ua)  <- rownames(foo_mat)
foo_out <- cbind(foo_mat[,colnames(foo_mat) == "RA"],
                 foo_mat[,colnames(foo_mat) == "NRA"],
                 foo_ua)
write.csv(foo_out, file="NewSignature_Data.csv")



foo_mat <- signature_exp
colnames(foo_mat) <- foo
foo_ua            <- exprs(raw_data_det)[rownames(raw_data_det) %in% filtered_in$ID,
                                         grep("UA", pData(raw_data_det)$Update_Collapse_2)]
colnames(foo_ua)  <- rep("UA", ncol(foo_ua))
rownames(foo_ua)  <- rownames(foo_mat)
foo_out <- cbind(foo_mat[,colnames(foo_mat) == "RA"],
                 foo_mat[,colnames(foo_mat) == "NRA"],
                 foo_ua)
write.csv(foo_out, file="NewDE_Data.csv")
##'-----------------------------------------------------------------------------------------#



##'Complex Heatmap
##'-----------------------------------------------------------------------------------------#
detach("package:ComplexHeatmap", unload=TRUE)
library(ComplexHeatmap)
library(circlize)


signature_exp                 <- exprs(raw_data_det)[rownames(raw_data_det) %in% signature,
                                                     grep("RA", pData(raw_data_det)$Update_Collapse_2)]

signature_exp                 <- exprs(raw_data_det)[rownames(raw_data_det) %in% signature[c(3,5,10)],
                                                     grep("RA", pData(raw_data_det)$Update_Collapse_2)]

signature_exp                 <- exprs(raw_data_det)[rownames(raw_data_det) %in% signature[c(3,10)],
                                                     grep("RA", pData(raw_data_det)$Update_Collapse_2)]


Rowv          <- colMeans(signature_exp)
d             <- dist(t(signature_exp))
h             <- hclust(d)
ddr           <- as.dendrogram(h)
ddr           <- reorder(ddr, Rowv)
ddrt          <- rev(ddr)


rm            <- rowMeans(signature_exp)
signature_exp <- sweep(signature_exp, 1, rm)
sx            <- apply(signature_exp, 1, sd)
signature_exp <- sweep(signature_exp, 1, sx, "/")

gene_names                    <- as.vector(anno_df[match(rownames(signature_exp),
                                                         anno_df$ID),]$symbol)
gene_names[is.na(gene_names)] <- "No Annotation"
rownames(signature_exp)       <- gene_names
foo                           <- as.vector(pData(raw_data_det)$Update_Collapse_2)
foo                           <- foo[grep("RA", foo)]

# colnames(signature_exp)       <- rep("", ncol(signature_exp))
df <- data.frame(Sample = foo)
ha = HeatmapAnnotation(df  = df, 
                       col = list(Sample = c("RA"  =  "black", 
                                             "NRA" =  "grey")))


ht2 = Heatmap(signature_exp, 
              name                     = "Row Z-Score", 
              column_title             = "", 
              bottom_annotation        = ha,
              col                      = colorRamp2(c(-2, 0, 2), 
                                                    c("blue", "white", "red")),
              show_row_hclust          = F,
              cluster_rows             = F,
              show_column_hclust       = T,
              show_column_names        = F,
              cluster_columns          = ddrt,
              column_hclust_height     = unit(1.5, "cm"),
              row_names_gp             = gpar(fontsize = 22)
)

png("Heatmap_Focus_V4.png", 
    width=3096, 
    height=2096, 
    units="px", 
    res=300)
draw(ht2, 
     show_heatmap_legend    = T, 
     show_annotation_legend = T,
     legend_grid_width      = unit(8, "mm"), 
     legend_title_gp        = gpar(fontsize = 18)
)
dev.off()
##'-----------------------------------------------------------------------------------------#



