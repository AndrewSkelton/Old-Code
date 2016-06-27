#-------------------------------------------------------------------------------------------#
#  Author      : Andrew J Skelton                                                           |
#  Language    : R Statistical Programming Language                                         |
#  Study       : Confirmation Study - Early Onset RA Study                                  |
#  Data Owner  : Newcastle University - Dr. Arthur Pratt                                    |
#  Description : Custom Plots                                                               |
#-------------------------------------------------------------------------------------------#



##'Set the Working Directory, load essential packages, and import files
##'-----------------------------------------------------------------------------------------#
library(gplots)
library(ggplot2)
library(reshape2)
library(scales)
##'-----------------------------------------------------------------------------------------#


##'Heatmap - Signature
##'-----------------------------------------------------------------------------------------#
png("Confirmation_Study_Heatmap_Signature.png", width=4096, height=3096, units="px", res=300)
anno_df$symbol <- as.vector(anno_df$symbol)
anno_df[match(signature[c(11,12)], anno_df$ID),]$symbol <- c("LOC731186", "GPRIN3")

#Useful
#nuID2IlluminaID(signature)

signature_exp                 <- batchCorrected_data[rownames(batchCorrected_data) %in% signature,]
gene_names                    <- as.vector(anno_df[match(rownames(signature_exp),
                                                         anno_df$ID),]$symbol)
gene_names[is.na(gene_names)] <- "No Annotation"
rownames(signature_exp)       <- gene_names
foo                           <- as.vector(pData(raw_data_det)$Baseline)
foo                           <- foo[grep("RA", foo)]
foo[foo == "RA"]              <- "black"
foo[foo == "NRA"]             <- "grey"
foo[foo == "UA"]              <- "white"
foo[foo == "OA"]              <- "grey"
par(cex.main=1.5)
heatmap.2(signature_exp,
          scale='row',
          dendrogram="none",
          trace='none',
          Rowv=NA,
          col=rev(redblue(100)),
          ColSideColors=foo,
          cexRow=1,
          labCol="",
          breaks=seq(-2,2,length.out=101),
          keysize = .9,
          density.info=c("none"),
          main="Signature Genes")
dev.off()
##'-----------------------------------------------------------------------------------------#



##'Signature Write
##'-----------------------------------------------------------------------------------------#
anno_df$symbol <- as.vector(anno_df$symbol)
anno_df[match(signature[c(11,12)], anno_df$ID),]$symbol <- c("LOC731186", "GPRIN3")

signature_exp                 <- batchCorrected_data[rownames(batchCorrected_data) %in% signature,]
gene_names                    <- as.vector(anno_df[match(rownames(signature_exp),
                                                         anno_df$ID),]$symbol)
gene_names[is.na(gene_names)] <- "No Annotation"
rownames(signature_exp)       <- gene_names
# 
df_exp                        <- batchCorrected_data[rownames(batchCorrected_data) %in% filtered_in$ID,]
gene_names                    <- as.vector(anno_df[match(rownames(df_exp),
                                                         anno_df$ID),]$symbol)
gene_names[is.na(gene_names)] <- "No Annotation"
rownames(df_exp)       <- gene_names

write.csv(signature_exp, file="Confirmation_Study_Signature_Normalised_Expression.csv")
write.csv(df_exp, file="Confirmation_Study_DE_Normalised_Expression.csv")
write.csv(pData(lumi.Q), file="Confirmation_Study_Pheno_Data.csv")

##'-----------------------------------------------------------------------------------------#



##'Volcano
##'-----------------------------------------------------------------------------------------#
png("Volcano_Confirmation_Study_Signature_andDE.png", width=4096, height=3096, units="px", res=300)
afc                <- 1.2
pval               <- 0.05
filtered_in        <- topTable(fit2, coef=1, p.value=pval, lfc=log2(afc),
                               number=Inf, adjust.method="none")
unfiltered_in      <- topTable(fit2, coef=1, number=Inf,
                               adjust.method="none")
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
  ggtitle("Confirmation Study: NRA-RA, UA Exclude\nP < 0.05 | FC > log2(1.2) | BH Uncorrected") +
  theme_bw()
dev.off()
##'-----------------------------------------------------------------------------------------#















