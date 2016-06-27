#-------------------------------------------------------------------------------------------#
#  Author      : Andrew J Skelton                                                           |
#  Language    : R Statistical Programming Language                                         |
#  Study       : miRNA Sequencing of Pooled Bovine Cartilage Samples                        |
#  Data Owner  : Newcastle University - Prof. David Young, (Cardiff Collaboration           |
#  Description : Useful Functions                                                           |
#-------------------------------------------------------------------------------------------#



##'Image Deploy
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



##'ggPCA
##'-----------------------------------------------------------------------------------------#
ggpca <- function(matrix_in, pheno, labels=1, colours=2, Third=T, objRet=F) {
  pca         <- prcomp(t(matrix_in))
  d           <- as.data.frame(pca$x)
  d$coA       <- pheno[,labels]
  d$coB       <- as.factor(pheno[,colours])
  
  gg <- ggplot(d, aes(x=PC1, y=PC2)) 
  if(Third == T) {
    gg <- gg + geom_point(aes(colour=coB, size=PC3)) 
  } else {
    gg <- gg + geom_point(aes(colour=coB), size=3) 
  }
  gg <- gg + theme_bw() +
    theme(legend.text = element_text(size = 12), text = element_text(size=12)) + 
    theme(panel.border = element_rect(size=1, colour = "black"))+ 
    geom_text(label=d$coA, size=4, vjust=1.2, hjust=-0.2) +
    theme(axis.text.x=element_text(vjust=0.5, size=10), 
          axis.text.y=element_text(size=10), axis.title.y=element_text(size=6)) 
  #                scale_x_continuous(limits=c(-30,50))
  if(objRet == T) {
    return(gg)
  } else {
    print(gg)
  }
}
##'-----------------------------------------------------------------------------------------#



##'Plot a Pairing
##'-----------------------------------------------------------------------------------------#
plot_paired <- function(mirna, counts, pheno, comparison="") {
  subset           <- counts[mirna,]
  subset           <- melt(subset)
  subset$Pairing   <- pheno$Prefix
  subset$TimePoint <- pheno$TimePoint
  subset$Pressure  <- pheno$Pressure
  
  cbPalette <- c("#999999", "#E69F00", "#56B4E9", 
                 "#009E73", "#F0E442", "#0072B2", 
                 "#D55E00", "#CC79A7")
  
  gg <- ggplot(subset, aes(x=Pairing, 
                           y=value, 
                           fill=Pressure)) +
    geom_bar(stat="identity", 
             position="dodge") +
    ggtitle(paste0(comparison,
                   " - ",
                   mirna)) +
    labs(y="Normalised Counts") +
    scale_fill_brewer(palette="Spectral") +
    theme_bw()
  return(gg)
}
##'-----------------------------------------------------------------------------------------#



##'plot_paired_comparison
##'-----------------------------------------------------------------------------------------#
plot_paired_comp <- function(mirna_in, 
                             hairpin_counts, 
                             mature_counts, 
                             pheno) {
  hairpin_df           <- hairpin_counts[grep(paste0(mirna_in,"$"), 
                                              rownames(hairpin_counts)),]
  hairpin_df           <- melt(hairpin_df)
  hairpin_df$mirna     <- mirna_in
  hairpin_df$origin    <- "hairpin"
  hairpin_df$sample    <- rownames(hairpin_df)
  hairpin_df           <- hairpin_df[,c(2,4,1,3)]
  
  mature_df            <- mature_counts[grep(paste(paste0(gsub("[r]",
                                                               "R",
                                                               mirna_in),
                                                          "$"),
                                                   paste0(gsub("[r]",
                                                               "R",
                                                               mirna_in),
                                                          "-"),
                                                   sep="|"),
                                             rownames(mature_counts)),]
  mature_df            <- melt(mature_df)
  mature_df$origin     <- "mature"
  colnames(mature_df)  <- colnames(hairpin_df)
  
  df_in                <- rbind(hairpin_df, 
                                mature_df)
  rownames(df_in)      <- 1:nrow(df_in)
  df_in$pairing        <- pheno[match(df_in$sample, 
                                      pheno$Sample_ID),]$Prefix
  df_in$pressure       <- pheno[match(df_in$sample, 
                                      pheno$Sample_ID),]$Pressure
  
  cbPalette <- c("#999999", "#E69F00", "#56B4E9", 
                 "#009E73", "#F0E442", "#0072B2", 
                 "#D55E00", "#CC79A7")
  
  gg <- ggplot(df_in, 
               aes(x=pairing, 
                   y=value, 
                   fill=pressure)) +
    geom_bar(stat="identity", 
             position="dodge") +
    scale_fill_brewer(palette="Spectral") +
    facet_grid(origin ~ mirna) +
    theme_bw()
  return(gg)
}
##'-----------------------------------------------------------------------------------------#




