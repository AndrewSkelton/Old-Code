#-------------------------------------------------------------------------------------------#
#  Author      : Andrew J Skelton                                                           |
#  Language    : R Statistical Programming Language                                         |
#  Optimised   : Sjogren's Syndome | Exiqon miRNA Microarray (2 Colour)                     |
#  Data Owner  : Newcastle University - Fai | Jessica Tarn                                  |
#  Description : Analysis Pipeline for Normalisation of a 2 Colour array using background   |
#                Correction, Within array and between array normalisation. Contrasts and    |
#                Gene Lists are generated using the Limma Package.                          |
#  Plots       : Boxplot | Chip Pseudo Images | Cy3/Cy5 Intensity Scatter | MA Plots |      |
#                Volcano Plots | Density Plots |                                            |
#-------------------------------------------------------------------------------------------#

source("http://bioconductor.org/biocLite.R")
biocLite()
biocLite(c("lumi", "gplots", "ggplot2", "limma", "ExiMiR"))
library(lumi)
library(gplots)
library(ggplot2)
library(limma)
library(ExiMiR)

setwd("~/Documents/Bioinformatics/Customers/Jess/May_2014_Exiqon_miRNA_Microarray/Raw_Data/Array_1/")
targets        <- readTargets(path = "./" )
RGList         <- read.maimages(targets[,c("Cy3", "Cy5")], source = "imagene", path = "./" )
RGList$genes   <- readGAL(path = "./")
RGList$printer <- getLayout(RGList$genes)

##'===========================================================================================================
##'
##' RAW INTENSITIES PLOTS
##'
pdf("Raw_Intensities.pdf")

x_lab_write <- function(x_labs)
{
  axis(1, at=seq(1, length(x_labs), by=1), labels = FALSE)
  text(seq(1+.25, length(x_labs)+.75, by=1), par("usr")[3]-.4, 
       labels = x_labs, srt = -90, pos = 1, xpd = TRUE, cex=0.4)
}

# treatment <- c("Co", "Co", "Co", "Co", "Co", "Co", "Ly", "Ly", "Ly", "Ly", "Ly", "Ly", "Ly", "Ly", "Ly", 
#                "Ly", "Ly", "Ly", "Hi", "Hi", "Hi", "Hi", "Hi", "Hi", "Hi", "Hi", "Hi", "Hi", "Hi", "Hi")

treatment <- c("Co", "Co", "Co", "Co", "Co", "Co", "Co", "Co", "Co", "Co", "Co", "Co", "Ly", "Ly", "Ly", "Ly", 
               "Ly", "Ly", "Ly", "Ly", "Ly", "Ly", "Ly", "Ly", "Hi", "Hi", "Hi", "Hi", "Hi", "Hi", "Hi", "Hi", "Hi", "Hi", "Hi", "Hi")

treatment_s <- treatment
treatment_s[grep('Hi|Lo', treatment_s)] <- 'Sj'

x_labs    <- c(paste(treatment, "_Cy3", sep=""), paste(treatment, "_Cy5", sep=""))

boxplot(data.frame(cbind(log2(RGList$G),log2(RGList$R))), cex=.2,
        main="Red & Green Channel Intensities", xaxt="n",
        col=rep(c("green","red"),each=length(treatment)))
x_lab_write(x_labs)

boxplot(data.frame(cbind(log2(RGList$Gb),log2(RGList$Rb))), cex=.2,
        main="Red & Green Background Channel Intensities", xaxt="n",
        col=rep(c("green","red"),each=length(treatment)))
x_lab_write(x_labs)

for(i in 1:length(colnames(RGList$Gb)))
{
  par(mfrow = c(2, 3))
  imageplot(log2(RGList$R[,i]),RGList$printer,  cex.main=.5, 
            main=paste("Array", i, "(", treatment[i] ,
                       ") \nChip Pseudo-Image - Red Channel", sep=" "))
  imageplot(log2(RGList$G[,i]),RGList$printer,  cex.main=.5, 
            main=paste("Array", i, "(", treatment[i] ,
                       ") \nChip Pseudo-Image - Green Channel", sep=" "))
  imageplot(log2(RGList$Rb[,i]),RGList$printer, cex.main=.5, 
            main=paste("Array", i, "(", treatment[i] ,
                       ") \nChip Pseudo-Image - Red Background Channel", sep=" "))
  imageplot(log2(RGList$Gb[,i]),RGList$printer, cex.main=.5, 
            main=paste("Array", i, "(", treatment[i] ,
                       ") \nChip Pseudo-Image - Green Background Channel", sep=" "))
  
  plot(log2(RGList$Rb[,i]),log2(RGList$R[,i]), cex=.2, cex.main=.5,
       main=paste("Scatter Plot of Array", i, "(", treatment[i], ") \nChip Pseudo-Image - Red Channel Signal Vs Background"))
  lines(c(-9,99),c(-9,99),col=2)
  
  plot(log2(RGList$Gb[,i]),log2(RGList$G[,i]), cex=.2, cex.main=.5,
       main=paste("Scatter Plot of Array", i, "(", treatment[i], ") \nChip Pseudo-Image - Green Channel Signal Vs Background"))
  lines(c(-9,99),c(-9,99),col=2)
}

dev.off()
##'

##'===========================================================================================================
##'
##' BACKGROUND CORRECTION
##'
RG         <- read.maimages(targets[,c("Cy3", "Cy5")], source = "imagene", path = "./")
RG$genes   <- readGAL(path = "./")
RG$printer <- getLayout(RGList$genes)
RG         <- backgroundCorrect(RGList, printer=RGList$printer, method="normexp", verbose=T, offset=50)
##'

##'===========================================================================================================
##'
##' BACKGROUND CORRECTED PLOTS
##'
pdf("BG_Corr_Intensities.pdf")

boxplot(log2(data.frame( RGList_BG_Corr$R,RGList_BG_Corr$G))[,as.vector(rbind(1:30,31:60))],
        col=rep(c("red","green"),30), xaxt="n", cex=.2, main="Background Corrected Intensities")
axis(1, at=seq(1, length(x_labs), by=2), labels=treatment, cex=.5)

plotDensities(RGList_BG_Corr)

par(mfrow = c(2, 2))
for(i in 1:length(colnames(RGList_BG_Corr$R)))
{
  plotMA(RGList_BG_Corr[,i], 
         main=paste("MA Plot of Array", i, "(", treatment[i], ")"))
}

dev.off()
##'

##'===========================================================================================================
##'
##' NORMALISATION
##'
MA=normalizeWithinArrays(RG,method="loess") 
MA=normalizeWithinArrays(MA,method="median")
normali
##'

##'===========================================================================================================
##'
##' NORMALISED PLOTS
##'
pdf("Normalised_Intensities.pdf")

boxplot(data.frame(MA$M), xaxt="n", cex=.2, main="Normalised Intensities")
axis(1, at=seq(1, length(x_labs), by=2), labels=treatment, cex=.5)

plotDensities(MA)

par(mfrow = c(2, 2))
for(i in 1:length(colnames(MA)))
{
  plotMA(MA[,i], 
         main=paste("MA Plot of Normalised Array", i, "(", treatment[i], ")"))
}

par(mfrow = c(2, 2))
for(i in 1:length(colnames(RGList_BG_Corr$G)))
{
  imageplot(log2(RGList_BG_Corr$R[,i]),RGList_BG_Corr$printer,  cex.main=.5, 
            main=paste("Normalised Array", i, "(", treatment[i] ,
                       ") \nChip Pseudo-Image - Red Channel", sep=" "))
  imageplot(log2(RGList_BG_Corr$G[,i]),RGList_BG_Corr$printer,  cex.main=.5, 
            main=paste("Normalised Array", i, "(", treatment[i] ,
                       ") \nChip Pseudo-Image - Green Channel", sep=" "))
}

dev.off()
##'

##'===========================================================================================================
##'
##' DIFFERENTIAL EXPRESSION
##'
treatments       <- c("Co", "Ly", "Hi")
array_names      <- treatment
# treatments       <- c("Co", "Ly", "Sj")
# array_names      <- treatment_s
design           <- model.matrix(~0 + factor(array_names, levels = treatments))
colnames(design) <- treatments
num_parameters   <- ncol(design)
fit              <- lmFit(MA$M, design)

cont_mat <- makeContrasts(Co-Ly, Co-Hi, Hi-Ly, levels=treatments)
#cont_mat <- makeContrasts(Co-Sj, Ly-Sj,  levels=treatments)
fit2 <- contrasts.fit(fit, contrasts=cont_mat)
fit2 <- eBayes(fit2)
##'

##'===========================================================================================================
##'
##' TABLES
##'
comparisons <- c("Co - Ly", "Co - Hi", "Hi - Ly")
#comparisons <- c("Co - Sj", "Ly - Sj")
p_cut_off   <- 0.01
fold_change <- 2
i           <- 1

fit2$genes <- MA$genes
for(i in 1:length(comparisons))#
{#
  gene_list            <- topTable(fit2, coef = comparisons[i], p.value = p_cut_off, lfc = log2(fold_change),
                                   nrow(fit2$genes), adjust.method = "BH")
  print(length(gene_list$logFC))
}
  
interesting <- c()
for(i in 1:length(comparisons))#
{#
  foo <- topTable(fit2, coef = comparisons[i], adjust.method="BH", number = nrow(fit2$genes))
  fooA <- foo[grep('miR-335$|miR-21$', foo$Name),]
  fooA$Comparison <- comparisons[i]
  interesting <- rbind(interesting, fooA)
}
write.csv(interesting, "MiR_of_interest.csv")

library(ggplot2)

pdf("volcano_plots.pdf")#
for(i in 1:length(comparisons))#
{#
  topTable(fit2, coef = comparisons[i], p.value = p_cut_off, lfc = log2(fold_change), adjust.method="BH")
  gene_list_unfiltered <- topTable(fit2, coef = comparisons[i], number = nrow(fit2$genes), adjust.method="BH")
  gene_list            <- topTable(fit2, coef = comparisons[i], p.value = p_cut_off, lfc = log2(fold_change),
                                   nrow(fit2$genes), adjust.method = "BH")
  
  comparison_split <- strsplit(comparisons[i], split=" ")
  #controls      <- batchCorrected_data[as.character(gene_list$probe_list),
  #                                     grep(comparison_split[[1]][1], array_names)]
  #conditions    <- batchCorrected_data[as.character(gene_list $probe_list),
  #                                     grep(comparison_split[[1]][3], array_names)]
  
  #gene_list$ControlMean   <- apply(controls,   1, mean)
  #gene_list$ConditionMean <- apply(conditions, 1, mean)
  
  if(length(gene_list)>0)
  {
    write.csv(gene_list, paste("gene_list_", comparison_split[[1]][1], "-", comparison_split[[1]][3], ".csv",
                               sep="",collapse=""))
  }
  
  #}#
  
  ##'VOLCANO_PLOT----------------------------------------------------------------------------------------------#
  
  gene_list_unfiltered$threshold = as.factor(abs(gene_list_unfiltered$logFC) > log2(fold_change) &
                                               gene_list_unfiltered$adj.P.Val < p_cut_off)
  g = ggplot(data=gene_list_unfiltered, aes(x=gene_list_unfiltered$logFC,
                                            y=-log10(gene_list_unfiltered$adj.P.Val),
                                            colour=gene_list_unfiltered$threshold)) +
    geom_point(alpha=0.4, size=1.75) +
    theme(legend.position = "none") +
    #xlim(c(-4, 4)) + ylim(c(0, 4.5)) +
    xlab("log2 fold change") + ylab("-log10 Adj.P-value") +
    ggtitle(paste("Volcano Plot of Fold Change vs Adj.P-Value", comparisons[i], sep=" "))
  print(g)
  
  
}#
dev.off()#
##'----------------------------------------------------------------------------------------------------------#
##'
