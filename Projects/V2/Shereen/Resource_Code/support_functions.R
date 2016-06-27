#-------------------------------------------------------------------------------------------#
#  Author      : Andrew J Skelton                                                           |
#  Language    : R Statistical Programming Language                                         |
#  Study       : Globin Clear vs Full RNA                                                   |
#  Data Owner  : Newcastle University - Prof. Fai Ng, Shereen Al-Ali                        |
#  Description : Support Functions for use in analysing Illumina HT-12 v4 Microarray Data   |
#-------------------------------------------------------------------------------------------#




##'PCA Plots Utilising the ggplot2 library
##'-----------------------------------------------------------------------------------------#
ggpca <- function(matrix_in, pheno, labels=1, colours=2, Third=T, objRet=F) {
  pca         <- prcomp(t(matrix_in))
  d           <- as.data.frame(pca$x)
  d$coA       <- pheno[,labels]
  d$coB       <- as.factor(pheno[,colours])

  gg          <- ggplot(d, aes(x=PC1, y=PC2))
  if(Third == T) {
    # 3D Illusion (depth of PC3 based on size of point)
    gg        <- gg + geom_point(aes(colour=coB, size=PC3))
  } else {
    gg        <- gg + geom_point(aes(colour=coB), size=3)
  }
  gg          <- gg + theme_bw() +
                      theme(legend.text=element_text(size = 12),
                            text=element_text(size=12)) +
                            theme(panel.border=element_rect(size=1, colour = "black")) +
                            geom_text(label=d$coA, size=4, vjust=1.2, hjust=-0.2) +
                            theme(axis.text.x=element_text(vjust=0.5, size=10),
                                  axis.text.y=element_text(size=10),
                                  axis.title.y=element_text(size=6))
  if(objRet == T) {
    return(gg)
  } else {
    print(gg)
  }
}
##'-----------------------------------------------------------------------------------------#




##'Dendrogram Utilising the ggplot2 library
##'-----------------------------------------------------------------------------------------#
gg_dendro = function(matrix_in, pheno, pos) {
  require(ggdendro)
  d          <- dist(t(matrix_in), method='euclidian')
  h          <- hclust(d, method='complete')
  dd         <- dendro_data(h)
  labs       <- label(dd)
  labs$group <- pheno[match(labs$label, rownames(pheno)),pos]

  g          <- ggdendrogram(h, leaf_labels=F, labels=F, rotate=T) +
                             scale_x_discrete(labels = c()) +
                             scale_y_continuous(expand = c(0.1, 0)) +
                             geom_text(data=labs, aes(label=label, x=x, y=y,
                                       colour=as.factor(group)), hjust=1.0, size=2.0) +
                             theme(legend.text=element_text(size=6),
                                   legend.title=element_text(size=6))

  return(g)
}
##'-----------------------------------------------------------------------------------------#




##'Volcano Plot Utilising the ggplot2 library
##'-----------------------------------------------------------------------------------------#
gg_volcano <- function(unfiltered_toptable, afc=2, pval=0.05, title=NA, report=FALSE) {
  unfiltered_toptable$pass <- unfiltered_toptable$adj.P.Val < pval &
                              abs(unfiltered_toptable$logFC) > log2(afc)
  g                        <- ggplot(unfiltered_toptable, aes(x=logFC, y=-log(adj.P.Val, 10)))
  g                        <- g + geom_point(aes(colour=pass), show_guide=F) +
                                  scale_colour_manual(values=c(alpha('black', 0.5), 'red')) +
                                  geom_hline(yintercept=-log(pval, 10), colour="red", linetype=2) +
                                  geom_vline(xintercept=-log2(afc), colour="red", linetype=2) +
                                  geom_vline(xintercept=log2(afc), colour="red", linetype=2) +
                                  theme_bw()
  if (!is.na(title)) g     <- g + ggtitle(title)
  if (report) {
    g                      <- g + theme(plot.title=element_text(size=6),
                                        axis.title=element_text(size=6))
  }
  return(g)
}
##'-----------------------------------------------------------------------------------------#




##'Helper Function to write high and lower resolution images based on a ggplot2 object
##'-----------------------------------------------------------------------------------------#
image_deploy <- function(ggobj, name="Default") {
  png(paste(name, "_HighRes.png", sep=""), width=12.53, height=6.98, units="in", res=600)
  print(ggobj)
  dev.off()
  png(paste(name, "_LowRes.png", sep=""), width=1400)
  print(ggobj)
  dev.off()
}
##'-----------------------------------------------------------------------------------------#




##'Helper Function to run a hypergeometric test that looks for over-represented GO Terms
##'based on a set of differentially expressed genes. The function tests for BP, MF, and CC
##'-----------------------------------------------------------------------------------------#
goTerms <- function(contrast, fit2, pval, fc) {
  base <- getwd()
  dir.create(file.path(base, "GO"), showWarnings = FALSE)
  setwd(file.path(base, "GO"))

  library("GOstats")
  sig_genes = topTable(fit2, coef=contrast, number=Inf, p.value=pval, lfc=log2(fc))
  if(nrow(sig_genes) > 4) {
    sig_probes = as.character(sig_genes$probe_list)
    entrez = unique(unlist(lookUp(nuIDs[sig_probes], "lumiHumanAll.db", "ENTREZID")))
    entrez = entrez[!is.na(entrez)]
    entrez = as.character(entrez)
    entrez_universe = unique(unlist(lookUp(nuIDs, "lumiHumanAll.db", "ENTREZID")))
    entrez_universe = entrez_universe[!is.na(entrez_universe)]
    entrez_universe = as.character(entrez_universe)
    bp_params = new("GOHyperGParams", geneIds=entrez, universeGeneIds=entrez_universe, annotation="lumiHumanAll.db",
                    ontology="BP", pvalueCutoff=0.05, conditional=FALSE, testDirection="over")
    bp_hyperg = hyperGTest(bp_params)
    bp_fdr = p.adjust(pvalues(bp_hyperg), method="fdr")
    bp_table = summary(bp_hyperg)
    bp_table$fdr = bp_fdr[bp_table$GOBPID]
    sig_bp_table = bp_table[bp_table$fdr < 0.05,]
    if (nrow(sig_bp_table) > 0) {
      write.csv(sig_bp_table, file=paste("GOTerms_BP_", contrast, ".csv", sep=""))
    }

    mf_params = new("GOHyperGParams", geneIds=entrez, universeGeneIds=entrez_universe, annotation="lumiHumanAll.db",
                    ontology="MF", pvalueCutoff=0.05, conditional=FALSE, testDirection="over")
    mf_hyperg = hyperGTest(mf_params)
    mf_fdr = p.adjust(pvalues(mf_hyperg), method="fdr")
    mf_table = summary(mf_hyperg)
    mf_table$fdr = mf_fdr[mf_table$GOMFID]
    sig_mf_table = mf_table[mf_table$fdr < 0.05,]
    if (nrow(sig_mf_table) > 0) {
      write.csv(sig_mf_table, file=paste("GOTerms_MF_", contrast, ".csv", sep=""))
    }

    cc_params = new("GOHyperGParams", geneIds=entrez, universeGeneIds=entrez_universe, annotation="lumiHumanAll.db",
                    ontology="CC", pvalueCutoff=0.05, conditional=FALSE, testDirection="over")
    cc_hyperg = hyperGTest(cc_params)
    cc_fdr = p.adjust(pvalues(cc_hyperg), method="fdr")
    cc_table = summary(cc_hyperg)
    cc_table$fdr = cc_fdr[cc_table$GOCCID]
    sig_cc_table = cc_table[cc_table$fdr < 0.05,]
    if (nrow(sig_cc_table) > 0) {
      write.csv(sig_cc_table, file=paste("GOTerms_CC_", contrast, ".csv", sep=""))
    }
  }
  setwd(base)
}
##'-----------------------------------------------------------------------------------------#




##'UPDATED FUNCTION'
##'Helper Function to run a hypergeometric test that looks for over-represented GO Terms
##'based on a set of differentially expressed genes. The function tests for BP, MF, and CC
##'-----------------------------------------------------------------------------------------#
goTerms <- function(contrast, fit2, pval, fc) {
  library(KEGGREST)
  library(org.Hs.eg.db)
  library(Category)
  library(pathview)
  library(GO.db)
  library(GOstats)

  gene_list_unfil   <- topTable(fit2, coef=contrast, number=Inf, adjust.method="BH")
  gene_list         <- topTable(fit2, coef=contrast, p.value=pval,
                                lfc=log2(fc), number=Inf, adjust.method = "BH")
  significant.genes <- as.vector(entrez_map[match(as.vector(gene_list$ID),
                                                  as.vector(entrez_map$nuID)),]$EntrezID)
  all.geneIDs       <- unique(entrez_map$EntrezID)

  hyperg            <- Category:::.doHyperGInternal
  hyperg.test       <-function(pathway.genes, significant.genes, all.genes, over=TRUE) {
    white.balls.drawn     <- length(intersect(significant.genes, pathway.genes))
    white.balls.in.urn    <- length(pathway.genes)
    total.balls.in.urn    <- length(all.genes)
    black.balls.in.urn    <- total.balls.in.urn - white.balls.in.urn
    balls.pulled.from.urn <- length(significant.genes)
    hyperg(white.balls.in.urn,    black.balls.in.urn,
           balls.pulled.from.urn, white.balls.drawn, over)
  }

  test_type               <- c("BP", "MF", "CC")
  for(i in 1:length(test_type)) {
    bp_params               <- new("GOHyperGParams",
                                   geneIds=unique(significant.genes),
                                   universeGeneIds=unique(as.vector(all.geneIDs)),
                                   annotation="lumiHumanAll.db",
                                   ontology=test_type[i],
                                   pvalueCutoff=0.05,
                                   conditional=FALSE,
                                   testDirection="over")
    bp_hyperg               <- hyperGTest(bp_params)
    bp_table                <- summary(bp_hyperg)
    bp_table$FDR            <- p.adjust(bp_table$Pvalue, method="fdr")
    bp_table                <- bp_table[bp_table$FDR < pval,]
    bp_table                <- bp_table[,c(1,2,8,3:7)]
    write.csv(bp_table, file=paste0(test_type[i], "_GO_Terms.csv"))
  }

}
##'-----------------------------------------------------------------------------------------#




##'Probe Profile Utilising the ggplot2 library
##'-----------------------------------------------------------------------------------------#
probeProfile <- function(matrix_in, pheno, contrast) {
  require(reshape2)
  require(scales)
  df           <- melt(matrix_in)
  df$Treatment <- pheno[df$Var2,]$treatment
  gg           <- ggplot(data=df, aes(x=Var2, y=value, colour=as.factor(Var1), group=Var1)) +
                         theme_bw() +
                         theme(axis.text.x=element_text(angle=90, vjust=0.5, size=10),
                               axis.text.y=element_text(size=11)) +
                         geom_point(size=2) + geom_line(size=1) +
                         ggtitle(paste(contrast, "")) +
                         facet_grid(Var1 ~ Treatment, space="free", scales="free_x") +
                         geom_line(stat="hline", yintercept="mean",
                                   colour="black", size=0.5) +
                         theme(strip.text.y=element_text(size = 11, angle = -90))
  return(gg)
}
##'-----------------------------------------------------------------------------------------#




##'Probe Profile Utilising the ggplot2 library
##'-----------------------------------------------------------------------------------------#
plot_profile <- function(expression_matrix, treatments=colnames(expression_matrix),
                         sep=FALSE, ret_obj=FALSE) {
  require(ggplot2)
  require(reshape2)
  require(scales)

  for_profile           <- as.data.frame(expression_matrix)
  for_profile$probe     <- row.names(for_profile)
  melted                <- melt(for_profile, id='probe')
  melted$facet_split    <- rep(treatments, each=nrow(for_profile))
  melted                <- melted[with(melted, order(probe)), ]

  melted$value_in_first <- rep(melted[melted$probe==melted$probe &
                               melted$variable==colnames(expression_matrix)[1],]$value,
                               each=ncol(expression_matrix))

  g                     <- ggplot(data=melted, aes(x=as.factor(variable), y=value)) +
                                  theme_bw() +
                                  theme(axis.text.x=element_text(angle=90,
                                                                 vjust=0.5,
                                                                 size=10)) +
                                  geom_point(aes(group=probe, colour=value_in_first)) +
                                  geom_line(aes(group=probe, colour=value_in_first)) +
                                  scale_x_discrete(breaks=levels(factor(melted$variable))) +
                                  scale_colour_gradient2(low='blue',
                                                         high='red',
                                                         mid='yellow',
                                                         midpoint=(max(melted$value_in_first) +
                                                                   min(melted$value_in_first))/2)

  if(sep == TRUE){
    g                  <- g + facet_grid(. ~ facet_split, space="free", scales="free_x")
  }

  if(ret_obj == TRUE) {
    return(g)
  } else {
    print(g)
  }
}
##'-----------------------------------------------------------------------------------------#




##'Intensity boxplots Utilising the ggplot2 library
##'-----------------------------------------------------------------------------------------#
norm_boxplot <- function(raw_data, normalised_data, pheno, pos) {
  group_factor <- factor(pheno[,pos], levels=unique(pheno[,pos]))
  rd           <- melt(raw_data)
  rd$value     <- log2(rd$value)
  rd$treat     <- 'raw'
  rd$group     <- rep(group_factor, each=nrow(normalised_data))
  nd           <- melt(normalised_data)
  nd$treat     <- 'norm'
  nd$group     <- rep(group_factor, each=nrow(normalised_data))
  ad           <- rbind(rd, nd)
  ad$treat     <- factor(ad$treat, levels = c("raw", "norm"))

  gg           <- ggplot(data=ad, aes(x=Var2, y=value)) +
                        geom_boxplot(aes(fill=group), outlier.size=0.5, size=0.2) +
                        facet_grid(. ~ treat) +
                        scale_x_discrete(name="") +
                        scale_y_continuous(name="Amplitude") +
                        theme_bw() +
                        theme(axis.text.x=element_text(angle=90, vjust=0.5, size=6),
                              axis.text.y=element_text(size=6),
                              axis.title.y=element_text(size=6),
                              legend.text=element_text(size=6),
                              legend.position="top",
                              legend.title=element_text(size=6))
  return(gg)
}
##'-----------------------------------------------------------------------------------------#




##'Helper Function to look for over-represented pathways in KEGG
##'-----------------------------------------------------------------------------------------#
keggGen <- function(contrast, fit2, pval, fc) {
  base <- getwd()
  dir.create(file.path(base, "KEGG"), showWarnings = FALSE)
  setwd(file.path(base, "KEGG"))

  sig_genes = topTable(fit2, coef=contrast, number=Inf, p.value=pval, lfc=log2(fc))
  sig_probes = as.character(sig_genes$probe_list)
  entrez = unique(unlist(lookUp(nuIDs[sig_probes], "lumiHumanAll.db", "ENTREZID")))
  entrez = entrez[!is.na(entrez)]
  entrez = as.character(entrez)
  entrez_universe = unique(unlist(lookUp(nuIDs, "lumiHumanAll.db", "ENTREZID")))
  entrez_universe = entrez_universe[!is.na(entrez_universe)]
  entrez_universe = as.character(entrez_universe)
  kegg_params = new("KEGGHyperGParams", geneIds=entrez, universeGeneIds=entrez_universe, annotation="lumiHumanAll.db",
                    pvalueCutoff=0.05, testDirection="over")
  kegg_hyperg = hyperGTest(kegg_params)
  kegg_fdr = p.adjust(pvalues(kegg_hyperg), method="fdr")
  kegg_table = summary(kegg_hyperg)
  kegg_table$fdr = kegg_fdr[kegg_table$KEGGID]
  sig_kegg_table = kegg_table[kegg_table$fdr < 0.05,]

  write.csv(sig_kegg_table, file=paste("KEGG_Table_", contrast, ".csv", sep=""))

  genes_in        <- sig_genes$logFC
  fc_entrez       <- unlist(lookUp(as.character(sig_genes$ID),
                                "lumiHumanAll.db", "ENTREZID"))
  names(genes_in) <- fc_entrez
  genes_in        <- genes_in[na.omit(names(genes_in))]
  for(j in 1:length(sig_kegg_table$KEGGID))
  {
    pv <- pathview(gene.data=genes_in,
                   pathway.id=sig_kegg_table$KEGGID[j],
                   species="hsa",
#                    limit=list(gene=c(round(min(sig_genes$logFC)),
#                                      round(max(sig_genes$logFC))), cpd=1),
                  limit=list(gene=c(-4,4), cpd=1),
                   low = list(gene = "blue", cpd="blue"),
                   mid = list(gene = "light grey", cpd="light grey"),
                   high = list(gene = "red", cpd="red")
                   )
  }
  dir.create(file.path(base, "KEGG/TMP"), showWarnings = FALSE)
  file.copy(list.files(".", "*path*"), "./TMP")
  file.remove(list.files(".", "*.xml"))
  file.remove(list.files(".", "*.png"))
  setwd(file.path(base, "KEGG/TMP"))
  file.copy(list.files(".", "*path*"), "../")
  setwd(base)
  unlink("./KEGG/TMP", T, T)
}
##'-----------------------------------------------------------------------------------------#



##'UPDATED FUNCTION
##'Helper Function to look for over-represented pathways in KEGG
##'-----------------------------------------------------------------------------------------#
pathways          <- keggList("pathway", "hsa")
human.pathways    <- sub("path:", "", names(pathways))
foo               <- c()
for(i in 1:length(pathways)) {
  foo             <- c(foo, setNames(keggGet(human.pathways[i]), human.pathways[i]))
}
pathways.all      <- foo

genes.by.pathway  <- lapply(pathways.all, function(pathways.all) {
  pathways.all$GENE[c(TRUE, FALSE)]
})
all.geneIDs       <- unique(entrez_map$EntrezID)
significant.genes <- as.vector(entrez_map[match(as.vector(gene_list$ID),
                                                as.vector(entrez_map$nuID)),]$EntrezID)

hyperg <- Category:::.doHyperGInternal
hyperg.test <-function(pathway.genes, significant.genes, all.genes, over=TRUE) {
    white.balls.drawn     <- length(intersect(significant.genes, pathway.genes))
    white.balls.in.urn    <- length(pathway.genes)
    total.balls.in.urn    <- length(all.genes)
    black.balls.in.urn    <- total.balls.in.urn - white.balls.in.urn
    balls.pulled.from.urn <- length(significant.genes)
    hyperg(white.balls.in.urn,    black.balls.in.urn,
           balls.pulled.from.urn, white.balls.drawn, over)
}

pVals.by.pathway       <- t(sapply(genes.by.pathway, hyperg.test,
                                   significant.genes, all.geneIDs))
table_out              <- as.data.frame(pVals.by.pathway)
table_out$FDR          <- p.adjust(table_out$p, method="fdr")
table_out$size         <- 0
table_out$overlap      <- 0
for(i in 1:length(rownames(table_out))) {
  table_out$size[i]    <- (length(pathways.all[[grep(rownames(table_out)[i],
                                                     names(pathways.all))]]$GENE) / 2)
  table_out$overlap[i] <- length(intersect(significant.genes,pathways.all[[grep(rownames(table_out)[i],
                                                               names(pathways.all))]]$GENE))
}
table_out              <- table_out[table_out$FDR < 0.05,]

gene_data              <- gene_list$logFC
names(gene_data)       <- entrez_map[match(gene_list$ID, entrez_map$nuID),]$EntrezID

base                   <- getwd()
dir.create(file.path(base, "KEGG"), showWarnings=F)
setwd(file.path(base, "KEGG"))

for(i in 1:nrow(table_out)) {
  pathview(gene.data=gene_data,
           pathway.id=rownames(table_out)[i],
           species="hsa",
           out.suffix="pathway",
           kegg.native=T,
           limit=list(gene=c(floor(min(gene_list$logFC)),
                             ceiling(max(gene_list$logFC))), cpd=1))
}

dir.create(file.path(base, "KEGG/TMP"), showWarnings=F)
write.csv(as.matrix(table_out), file="Pathway_Analysis_Table.csv")
file.copy(list.files(".", "*path*"), "./TMP")
file.remove(list.files(".", "*.xml"))
file.remove(list.files(".", "*.png"))
setwd(file.path(base, "KEGG/TMP"))
file.copy(list.files(".", "*path*"), "../")
setwd(base)
unlink("./KEGG/TMP", T, T)
##'-----------------------------------------------------------------------------------------#




##'Unity Function to tie everything together
##'-----------------------------------------------------------------------------------------#
downstreamBundle <- function(contrast, fit2, pval, fc, lumiSubset, expression_subset, folderID) {
  cat("Generating Directory Structure...")
  baseRoot <- getwd()
  dir.create(file.path(baseRoot, folderID), showWarnings = FALSE)
  setwd(file.path(baseRoot, folderID))
  cat("Done!\n")

  cat("Generating Gene Lists...")
  gene_list_unfiltered <- topTable(fit2, coef = contrast, number=Inf, adjust.method="BH")
  gene_list            <- topTable(fit2, coef = contrast, p.value=pval, lfc = log2(fc),
                                   number=Inf, adjust.method = "BH")

  base <- getwd()
  dir.create(file.path(base, "DE_Gene_Lists"), showWarnings = FALSE)
  setwd(file.path(base, "DE_Gene_Lists"))
  write.csv(gene_list_unfiltered, file=paste("UnfilteredGeneList_", folderID, ".csv", sep=""))
  write.csv(gene_list, file=paste("DEGeneList_", folderID, ".csv", sep=""))
  setwd(base)
  cat("Done!\n")

  cat("Making Volcano Plot...")
  ggVol <- gg_volcano(gene_list_unfiltered, fc, pval, contrast, F)
  dir.create(file.path(base, "Volcano_Plot"), showWarnings = FALSE)
  setwd(file.path(base, "Volcano_Plot"))
  image_deploy(ggVol, paste("VolcanoPlot_", contrast, sep=""))
  setwd(base)
  cat("Done!\n")

  if(nrow(gene_list) > 0) {
    cat("Finding GO Terms...")
    goTerms(contrast, fit2, pval, fc)
    cat("Done!\n")

    cat("Generating Probe Profile...")
    if(nrow(gene_list) > 10) {
      pp <- probeProfile(expression_subset[rownames(expression_subset) %in% gene_list$probe_list[1:10],],
                         pData(lumiSubset), contrast)
    } else if (nrow(gene_list) > 1) {
      pp <- probeProfile(expression_subset[rownames(expression_subset) %in% gene_list$probe_list,],
                         pData(lumiSubset), contrast)
    }
    if(nrow(gene_list)>1) {
      dir.create(file.path(base, "Probe_Profile"), showWarnings = FALSE)
      setwd(file.path(base, "Probe_Profile"))
      image_deploy(pp, paste("ProbeProfile_", contrast, sep=""))
    }
    setwd(base)
    cat("Done!\n")

    cat("Finding KEGG Data...")
    keggGen(contrast, fit2, pval, fc)
    cat("Done!\n")
  }

  setwd(baseRoot)
}
##'----------------------------------------------------------------------------------------------------------#
