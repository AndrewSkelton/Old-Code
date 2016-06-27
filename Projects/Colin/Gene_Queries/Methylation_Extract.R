#-------------------------------------------------------------------------------------------#
#  Author      : Andrew J Skelton                                                           |
#  Language    : R                                                                          |
#  Study       : MCF2L Data from 450K                                                       |
#  Data Owner  : Louise Reynard                                                             |
#  Description : Visualise the methylation of MCF2L - Requires 450K data						        |
#-------------------------------------------------------------------------------------------#



##'Set the Working Directory, load essential packages, and import files
##'-----------------------------------------------------------------------------------------#
library(ggplot2)
library(biomaRt)
library(gridExtra)
library(biovizBase)
library(GenomicFeatures)
library(ggbio)


biocLite(c("biovizBase", "GenomicFeatures", "ggbio"))
setwd("~/Colin/Gene_Queries/")
cpg_of_interest <- "MCF2L" #"ALDH1A2"

##'-----------------------------------------------------------------------------------------#



##'Build TxDb
##'-----------------------------------------------------------------------------------------#
transcriptdb <- makeTxDbFromGFF(file     = "/data/genomes/GRCh37_Ensembl_Homo_Sapiens/genome.gtf", 
                                format   = 'gtf', 
                                organism = 'Homo sapiens')
##'-----------------------------------------------------------------------------------------#



##'Extract Probe Data
##'-----------------------------------------------------------------------------------------#
MCF2L_Anno      <- annotation_other[grep(paste(paste0(cpg_of_interest, ";"),
                                               paste0(cpg_of_interest, "$"),
                                               sep="|"), 
                                         annotation_other$UCSC_RefGene_Name),]
MCF2L_Probes    <- hm450.hg19[rownames(hm450.hg19) %in% rownames(MCF2L_Anno),]
##'-----------------------------------------------------------------------------------------#



##'Query GRCh37
##'-----------------------------------------------------------------------------------------#
ensembl    <- useMart("ENSEMBL_MART_ENSEMBL",
                      dataset="hsapiens_gene_ensembl",
                      host="grch37.ensembl.org",
                      path="/biomart/martservice")
annotation <- getBM(attributes=c("ensembl_gene_id",
                                 "external_gene_name",
                                 "description"),
                    filters="external_gene_name",
                    values=cpg_of_interest,
                    ensembl)
##'-----------------------------------------------------------------------------------------#



##'Get Betas
##'-----------------------------------------------------------------------------------------#
MCF2L_Beta      <- m2beta(exprs(lumi.norm[rownames(lumi.norm) %in% rownames(MCF2L_Probes),]))
##'-----------------------------------------------------------------------------------------#



##'Create Data Structures
##'-----------------------------------------------------------------------------------------#
df_in           <- melt(MCF2L_Beta)
df_in           <- cbind(df_in, 
                         pData(lumi.norm)[match(df_in$Var2,
                                                rownames(pData(lumi.norm))),][,c(1,6)])
rownames(df_in) <- 1:nrow(df_in)
# df_in           <- cbind(df_in, pData(lumi.norm)[,c(1,6)])
anno_base_in    <- hm450.hg19[match(df_in$Var1, rownames(hm450.hg19)),]
anno_other_in   <- MCF2L_Anno[match(df_in$Var1, rownames(MCF2L_Anno)),]
df_in           <- cbind(df_in, anno_base_in[,c(1:3)])
# df_in           <- cbind(df_in, anno_other_in[,c(5,12:13)])
colnames(df_in) <- c("CpG", "Sample", "Beta", "Sample", 
                     "Type", "Chr", "Start", "Stop")
df_in           <- df_in[with(df_in, order(Start)),]
##'-----------------------------------------------------------------------------------------#



##'Plot Magic
##'-----------------------------------------------------------------------------------------#
pdf(paste0(cpg_of_interest, "_Methylation.pdf"), paper="a4r")
for(i in seq(1, length(unique(df_in$CpG)), 8)) {
  
  df_sub       <- df_in[df_in$CpG %in% as.vector(unique(df_in$CpG)[c(i:(i+7))]),]
  df_sub       <- df_sub[with(df_sub, order(Start)),]
  sort_df      <- data.frame(CpG   = unique(df_sub$CpG),
                             Start = unique(df_sub$Start))
  
  gg_a <- ggplot(data=df_sub, aes(x=Start, y=Beta, group=Type, colour=Type)) +
    scale_y_continuous(limits = c(0, 1)) +
    scale_x_continuous(limits = c(min(df_sub$Start), 
                                  max(df_sub$Start))) +
    geom_point(alpha = 0.4) +
    # geom_jitter() +
    # geom_boxplot() +
    theme_bw() +
    scale_fill_brewer(palette="Set1") + 
    scale_color_brewer(palette="Set1") +
    stat_summary(fun.y = mean, 
                 geom="line") + 
    ggtitle(cpg_of_interest)
  
  df_sub$CpG <- factor(df_sub$CpG, 
                         levels=sort_df$CpG)
  
  gg_b <- ggplot(data=df_sub, aes(x=CpG, y=Beta)) +
    scale_y_continuous(limits = c(0, 1)) +
    geom_boxplot(aes(fill = Type),
                 notch=T) +
    geom_point(aes(shape=Type, 
                   fill=Type),
               position=position_dodge(width=.75), 
               pch=21, 
               size=2) +
    theme_bw() +
    scale_fill_brewer(palette="Set1") +
    theme(axis.text.x = element_text(angle = 90, 
                                     hjust = 1))
  
  genename      <- annotation$ensembl_gene_id
  for_intersect <- transcriptsBy(transcriptdb, 
                                 by = "gene")[genename]
  trans         <- crunch(transcriptdb, 
                          which=for_intersect)
  gr1           <- split(trans, 
                         trans$tx_name)
  p1            <- autoplot(gr1) + 
                   theme_bw() +
                   geom_vline(xintercept=c(min(df_sub$Start),
                                           max(df_sub$Start)), 
                   colour="red")
  
  p2 <- p1@ggplot
  multiplot(p2, gg_a)
  
  
  
  
  
  fixed(p1)     <- T
  fixed(gg_a)   <- T
  fixed(gg_b)   <- T
  labeled(p1)   <- F
  labeled(gg_a) <- F
  labeled(gg_b) <- F
  foo         <- list(A1=p1, A2=gg_a, A3=gg_b)
  foo         <- list(p1, gg_a)
  tracks(foo)

  plot_grid(gg_a, gg_b, p2, labels=c("A","B","C"))
  multiplot(gg_a, gg_b, cols=1)
  
}
dev.off()
##'-----------------------------------------------------------------------------------------#



##'Save Outputs
##'-----------------------------------------------------------------------------------------#
write.csv(MCF2L_Beta, file=paste0(cpg_of_interest, "_Betas.csv"))
write.csv(pData(lumi.norm), file=paste0(cpg_of_interest, "_Pheno_Data.csv"))
write.csv(anno_base_in, file="Annotation_1.csv")
write.csv(anno_other_in, file="Annotation_2.csv")
##'-----------------------------------------------------------------------------------------#


cpg         <- c("cg13979708","cg19254793","cg20913747","cg18551225")
norm_betas  <- m2beta(exprs(lumi.norm[rownames(lumi.norm) %in% cpg,]))
anno        <- hm450.hg19[match(cpg, rownames(hm450.hg19)),]


plot_transcript_expression = function(norm_trans, norm_gene, gene_id, transcriptdb, pdata, sample=NA, ...) {
  
  
  
  trans1 = trans
  trans2 = trans
  rm(trans)
  trans1$celltype = levels(pdata)[1]
  trans1$expression = apply(norm_trans[match(trans1$tx_name,rownames(norm_trans)),which(pdata==levels(pdata)[1])], 1, mean) 
  
  trans2$celltype = levels(pdata)[2]
  trans2$expression = apply(norm_trans[match(trans2$tx_name,rownames(norm_trans)),which(pdata==levels(pdata)[2])], 1, mean)
  
  gr1 = split(trans1, trans1$tx_name)
  gr2 = split(trans2, trans2$tx_name)
  
  all_exp = c(trans1$expression, trans2$expression)
  
  p1 = autoplot(gr1, aes(fill=expression), colour=alpha('black', 0.2)) +
    scale_fill_gradient(low="blue", high="red", limits=c(0,max(all_exp))) + theme_clear()
  p2 = autoplot(gr2, aes(fill=expression), colour=alpha('black', 0.2)) +
    scale_fill_gradient(low="blue", high="red", limits=c(0,max(all_exp))) + theme_clear()
  track_col = gg_color_hue(2)
  track_col = c("white", track_col, "white")
  t = tracks(undiff=p1, diff=p2, title=paste0(gene_id, " (", genename, ")"),
             label.bg.fill=track_col)
  return(t)
}
