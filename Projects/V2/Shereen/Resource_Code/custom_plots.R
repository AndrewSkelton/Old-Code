#-------------------------------------------------------------------------------------------#
#  Author      : Andrew J Skelton                                                           |
#  Language    : R Statistical Programming Language                                         |
#  Study       : Globin Clear vs Full RNA                                                   |
#  Data Owner  : Newcastle University - Prof. Fai Ng, Shereen Al-Ali                        |
#  Description : Custom Plots for use in illustrating key analyses of the study             |
#-------------------------------------------------------------------------------------------#


##'PCA Code utilising the ggplot2 library
##'-----------------------------------------------------------------------------------------#
pca         <- prcomp(t(normalised_data))
d           <- as.data.frame(pca$x)
d           <- cbind(d, pData(raw_data_det))
d$pairs     <- as.factor(d$pairs)

gg1 <- ggplot(d, aes(x=PC1, y=PC2, shape=class)) +
       geom_point(size=3) +
       geom_text(label=d$pairs, size=4, vjust=1.2, hjust=-0.2) +
       theme_bw() +
       theme(axis.title.x = element_text(size=20),
             axis.title.y = element_text(size=20))

gg2 <- ggplot(d, aes(x=PC1, y=PC2, shape=class)) +
       geom_point(aes(colour=treatment), size=3) +
       geom_text(label=d$pairs, size=4, vjust=1.2, hjust=-0.2) +
       scale_colour_grey(start = 0, end = .9) +
       theme_bw() +
       theme(axis.title.x = element_text(size=20),
             axis.title.y = element_text(size=20))
##'-----------------------------------------------------------------------------------------#




##'Venn Diagrams of Detection Threshold numbers
##'-----------------------------------------------------------------------------------------#
library(VennDiagram)
normal_globin <- rownames(normalised_data)
normal_rna    <- rownames(normalised_data)
normal_raw    <- rownames(raw_data)

grid.newpage()
venn.plot <- draw.triple.venn(length(normal_globin),
                              length(normal_rna),
                              length(normal_raw),
                              length(intersect(normal_globin, normal_rna)),
                              length(intersect(normal_rna, normal_raw)),
                              length(intersect(normal_globin, normal_raw)),
                              length(intersect(intersect(normal_globin, normal_raw), normal_rna)),
                              c("Globin Clear", "Full RNA", "Total Probes on HT-12v4 Array"),
                              scaled=F)
grid.draw(venn.plot)

grid.newpage()
venn.plot <- draw.pairwise.venn(length(normal_globin),
                                length(normal_rna),
                                length(intersect(normal_globin, normal_rna)),
                                scaled=T,
                                c("Globin Clear", "Full RNA"))
grid.draw(venn.plot)
##'-----------------------------------------------------------------------------------------#




##'Custom Volcano Plots to show the difference between Full RNA and Globin Clear Preps
##'-----------------------------------------------------------------------------------------#
library(VennDiagram)
unfiltered_fullrna    <- topTable(fit2, coef="SS - control", number=Inf,
                                  adjust.method="BH")
unfiltered_fullrna_in <- unfiltered_fullrna
unfiltered_globin     <- topTable(fit2, coef="SS - control", number=Inf,
                                  adjust.method="BH")
unfiltered_globin_in  <- unfiltered_globin

afc  <- 1.2
pval <- 0.05

unfiltered_fullrna_in$pass  <- unfiltered_fullrna_in$adj.P.Val < pval & abs(unfiltered_fullrna_in$logFC) > log2(afc)
unfiltered_fullrna_in$class <- "Full_RNA"
unfiltered_globin_in$pass   <- unfiltered_globin_in$adj.P.Val < pval & abs(unfiltered_globin_in$logFC) > log2(afc)
unfiltered_globin_in$class  <- "Globin_Clear"
df                          <- rbind(unfiltered_fullrna_in, unfiltered_globin_in)

ggplot(df, aes(x=logFC, y=-log(adj.P.Val, 10))) +
      geom_point(aes(colour=pass, shape=pass), show_guide=F) +
      scale_colour_manual(values=c(alpha('grey', 0.5), 'black')) +
      geom_hline(yintercept=-log(pval, 10), colour="black", linetype=2) +
      geom_vline(xintercept=-log2(afc), colour="black", linetype=2) +
      geom_vline(xintercept=log2(afc), colour="black", linetype=2) +
      facet_grid(. ~ class) +
      theme_bw()

ggplot(unfiltered_fullrna_in, aes(x=logFC, y=-log(adj.P.Val, 10))) +
      geom_point(aes(colour=pass, shape=pass), show_guide=F) +
        scale_colour_manual(values=c(alpha('grey', 0.5), 'black')) +
      geom_hline(yintercept=-log(pval, 10), colour="black", linetype=2) +
      geom_vline(xintercept=-log2(afc), colour="black", linetype=2) +
      geom_vline(xintercept=log2(afc), colour="black", linetype=2) +
      ggtitle("Full RNA Differential Expression: PSS - Control") +
      theme_bw()

ggplot(unfiltered_globin_in, aes(x=logFC, y=-log(adj.P.Val, 10))) +
      geom_point(aes(colour=pass, shape=pass), show_guide=F) +
      scale_colour_manual(values=c(alpha('grey', 0.5), 'black')) +
      geom_hline(yintercept=-log(pval, 10), colour="black", linetype=2) +
      geom_vline(xintercept=-log2(afc), colour="black", linetype=2) +
      geom_vline(xintercept=log2(afc), colour="black", linetype=2) +
      ggtitle("Globin Clear Differential Expression: PSS - Control") +
      theme_bw()
##'-----------------------------------------------------------------------------------------#




##'Custom Heatmaps visualisng 25 Genes of Interest and DE Genes
##'-----------------------------------------------------------------------------------------#
#pdf("Custom_Heatmaps_Followup_Array.pdf", paper="a4r")
genes_of_interest       <- read.table("../Shereen_Sept_2014/25_genes.txt", sep='\t', header=T)
genes_of_interest       <- merge(genes_of_interest, anno_df, by.x="probe.ID", by.y="ID")
mat_in                  <- normalised_data[match(genes_of_interest$probe_list, rownames(normalised_data)),]
mat_in                  <- mat_in[,grep("SS|lymphoma", pData(lumi.Q)$Treatment)]
rownames(mat_in)        <- genes_of_interest$Gene.symbol
pheno_in                <- pData(lumi.Q)[match(colnames(mat_in), rownames(pData(lumi.Q))),]
foo                     <- pheno_in$Treatment
foo[foo == "SS"]        <- "black"
foo[foo == "lymphoma"]  <- "grey"

heatmap.2(mat_in,
          scale='row',
          dendrogram="column",
          trace='none',
          Rowv=NA,
          # col=my_palette,
          col=gray.colors(100),
          ColSideColors=foo,
          cexRow=.6,
          labCol="",
          keysize = .9,
          density.info=c("none"),
          main="25 Genes of Interest")


mat_in                  <- normalised_data[match(gene_list$probe_list, rownames(normalised_data)),]
mat_in                  <- mat_in[,grep("SS|lymphoma", pData(lumi.Q)$Treatment)]
rownames(mat_in)        <- gene_list$symbol
pheno_in                <- pData(lumi.Q)[match(colnames(mat_in), rownames(pData(lumi.Q))),]
foo                     <- pheno_in$Treatment
foo[foo == "SS"]        <- "black"
foo[foo == "lymphoma"]  <- "grey"

heatmap.2(mat_in,
          scale='row',
          dendrogram="column",
          trace='none',
          Rowv=NA,
          # col=my_palette,
          col=gray.colors(100),
          ColSideColors=foo,
          cexRow=.6,
          labCol="",
          keysize = .9,
          density.info=c("none"),
          main="DE Genes - pSS Vs Lymphoma")
#dev.off()
##'-----------------------------------------------------------------------------------------#
