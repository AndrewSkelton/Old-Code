#-------------------------------------------------------------------------------------------#
#  Author      : Andrew J Skelton                                                           |
#  Language    : R Statistical Programming Language                                         |
#  Study       : Interferon Study of a TET2 Mutation in pediatric Patients - Rare Samples   |
#  Data Owner  : Newcastle University - Prof. Sophie Hambleton                              |
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



##'Volcano Facet of Different Interferon Treatments
##'-----------------------------------------------------------------------------------------#
df_in <- c()
treat <- c("alpha", "alpha", "beta", "beta", "gamma", "gamma")
type  <- c("Stock", "Patient", "Stock", "Patient", "Stock", "Patient")
point <- 1
for(i in c(5,7,10,8,6,9)) {
  tmp        <- topTable(fit2, 
                         coef=comparisons[i], 
                         number=Inf, 
                         adjust.method="BH")
  tmp$Source <- type[point]
  tmp$Treat  <- treat[point]
  point      <- point + 1
  df_in      <- rbind(df_in, tmp)
}

p_cut_off    <- 0.01
fold_change  <- 2
mtc          <- 'BH'

df_in$cutoff <- df_in$adj.P.Val  < p_cut_off & 
                abs(df_in$logFC) > log2(fold_change)

gg           <- ggplot(df_in, aes(x=logFC, y=-log(adj.P.Val, 10))) + 
                  theme_bw() +
                  geom_point(aes(colour=cutoff), 
                             show_guide=F) + 
                  scale_colour_manual(values=c(alpha('black', 0.5), 
                                               'red')) +
                  geom_hline(yintercept=-log(p_cut_off, 10), 
                             colour="red", 
                             linetype=2) +
                  geom_vline(xintercept=-log2(fold_change), 
                             colour="red", 
                             linetype=2) +
                  geom_vline(xintercept= log2(fold_change), 
                             colour="red", 
                             linetype=2) +
                  facet_grid(Treat ~ Source, 
                             space="free", 
                             labeller=label_parsed) +
                  theme(strip.background = element_blank(), 
                        strip.text       = element_blank(),
                        axis.text        = element_text(size   = 18),
                        strip.text.y     = element_text(size   = 18, 
                                                        colour = "black", 
                                                        angle  = 0),
                        strip.text.x     = element_text(size   = 16, 
                                                        colour = "black"),
                        text             = element_text(size=22))

png("Multi_Volcano_Custom_V3.png", width=3096, height=4096, units="px", res=300)
print(gg)
dev.off()
##'-----------------------------------------------------------------------------------------#















