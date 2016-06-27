#-------------------------------------------------------------------------------------------#
#  Author      : Andrew J Skelton                                                           |
#  Language    : R Statistical Programming Language                                         |
#  Study       : si Knockout and Interlukin Treatment Time Series                           |
#  Data Owner  : Newcastle University - Prof. David Young                                   |
#  Description : Custom Plots                                                               |
#-------------------------------------------------------------------------------------------#



##'Set the Working Directory, load essential packages, and import files
##'-----------------------------------------------------------------------------------------#
library(gplots)
library(ggplot2)
library(reshape2)
library(scales)
##'-----------------------------------------------------------------------------------------#


##'PCA - Normalised
##'-----------------------------------------------------------------------------------------#
png("PCA_Normalised_Data_2.png", width=4096, height=3096, units="px", res=300)
pca          <- prcomp(t(normalised_data))
d            <- as.data.frame(pca$x)
d            <- cbind(d, pData(lumi.Q))

ggplot(d, aes(x=PC1, y=PC2)) + #, shape=Sample_ID
  geom_point(aes(colour=Treatment), size=3) +
  geom_text(label=d$Sample_ID, size=4, vjust=1.2, hjust=-0.2) +
  theme_bw() +
  ggtitle("PCA of Normalised Microarray Data") +
  theme(axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15))
dev.off()
##'-----------------------------------------------------------------------------------------#


##'ggVolcano Function
##'Parameters:
##'           Unfiltered Toptable from Limma contrast test
##'           Actual Fold Change Cutoff
##'           Ajusted P Value Cutoff
##'           Contrast In
##'Returns:
##'           A ggplot2 Object
##'-----------------------------------------------------------------------------------------#
ggvolcano <- function(unfiltered_toptable = NULL, 
                      actual_fc_cutoff    = 1.2, 
                      adj_p_cutoff        = 0.01,
                      contrast            = "") {
  
  unfiltered_toptable$Threshold <- unfiltered_toptable$adj.P.Val  < adj_p_cutoff & 
                                   abs(unfiltered_toptable$logFC) > log2(actual_fc_cutoff)
  gg      <-   ggplot(unfiltered_toptable, 
                      aes(x=logFC, 
                          y=-log(adj.P.Val, 10))) +
                      geom_point(aes(colour=Threshold, 
                                     shape=Threshold), 
                                 show_guide=F) +
                      scale_colour_manual(values=c(alpha('grey', 0.5), 
                                                   'red')) +
                      geom_hline(yintercept=-log(adj_p_cutoff, 10), 
                                 colour="black", 
                                 linetype=2) +
                      geom_vline(xintercept=-log2(actual_fc_cutoff), 
                                 colour="black", 
                                 linetype=2) +
                      geom_vline(xintercept=log2(actual_fc_cutoff), 
                                 colour="black", 
                                 linetype=2) +
                      ggtitle(contrast) +
                      theme_bw()
  return(gg)
}
##'-----------------------------------------------------------------------------------------#






