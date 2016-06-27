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

##'LIBRARIES--------------------------------------------------------------------------------#
source("http://bioconductor.org/biocLite.R")

biocLite(c("lumi", "gplots", "ggplot2", "limma", "ExiMiR", "MmPalateMiRNA", "targetscan.Hs.eg.db",
           "microRNA"))
biocLite()

library(microRNA)
library(targetscan.Hs.eg.db)
library(lumi)
library(gplots)
library(ggplot2)
library(limma)
library(affycoretools)
library(MmPalateMiRNA)
library(vsn)

setwd("~/Documents/Bioinformatics/Customers/Jess/May_2014_Exiqon_miRNA_Microarray/Raw_Data/Array_1/")

##'-----------------------------------------------------------------------------------------#

##'EXPERIMENTAL SETUP-----------------------------------------------------------------------#

treatment <- c("Co", "Ly", "Lo", "Lo", "Lo", "Ly", "Hi", "Lo", "Hi", "Hi", "Lo", "Hi", "Hi", 
               "Co", "Lo", "Ly", "Ly", "Ly", "Co", "Lo", "Hi", "Hi", "Co", "Ly", "Lo", "Co", 
               "Ly", "Co", "Ly", "Co", "Hi", "Co", "Co", "Lo", "Hi", "Ly")

#Array 1
treatment <- c("Co", "Co", "Co", "Co", "Co", "Co", "Co", "Co", "Co", "Co", "Co", "Co", "Ly", 
               "Ly", "Ly", "Ly", "Ly", "Ly", "Ly", "Ly", "Ly", "Ly", "Ly", "Ly", "Hi", "Hi",
               "Hi", "Hi", "Hi", "Hi", "Hi", "Hi", "Hi", "Hi", "Hi", "Hi")

##'-----------------------------------------------------------------------------------------#

##'READ IN RAW DATA-------------------------------------------------------------------------#

targets        <- readTargets(path = "./" )
RGList         <- read.maimages(targets[,c("Cy3", "Cy5")], source = "imagene", path = "./" )
RGList$genes   <- readGAL(path = "./")
RGList$printer <- getLayout(RGList$genes)

probes_in <- unique(substr(RGList$genes$Name, 1, 4))[-c(1,3,11)]

foo <- RGList$genes$Name
foo[grep('hsa*', foo)] <- "Human"
foo[grep(paste(probes_in[-c(1,3)], collapse="|"), foo)] <- "Other"
foo[grep('spike|Hy3', foo)] <- "Spike"
foo[foo == ""] <- "Empty"
RGList$genes$probe.type <- as.factor(foo)

RGList$genes$'Gene ID' <- RGList$genes$ID

##'CREATE BACKUP RAW DATA
RGList.bck    <- RGList
treatment.bck <- treatment

##'REMOVE POTENTIAL OUTLIER - 19715600
RGList <- RGList.bck
# RGList    <- RGList.bck[,-(grep('19715600', colnames(RGList.bck)))]
# treatment <- treatment.bck[-(grep('19715600', colnames(RGList.bck)))]
# RGList    <- RGList.bck[,-(grep('19715600|610', colnames(RGList.bck)))]
# treatment <- treatment.bck[-(grep('19715600|610', colnames(RGList.bck)))]
# RGList    <- RGList.bck[,-(grep('19715600|610|584', colnames(RGList.bck)))]
# treatment <- treatment.bck[-(grep('19715600|610|584', colnames(RGList.bck)))]

RGList    <- RGList.bck[,-(grep('14311357', colnames(RGList.bck)))]
treatment <- treatment.bck[-(grep('14311357', colnames(RGList.bck)))]
split   <- strsplit(colnames(RGList),'_')
short_colnames <- unique(grep('143', unlist(split), value=T))

##'-----------------------------------------------------------------------------------------#

##'RAW DATA VISUALISATION-------------------------------------------------------------------#

res <- levelplot(RGList[, c(1:36)],
                 channel="G", group="probe.type",
                 subset=c("Human", "Other", "Spike"),
                 scales = list(rot=c(45, 45)))

jpeg(filename="raw_profile_green.jpeg", height=779, width=1400, quality=100)
pp <- plot_profile(log2(RGList$G))
print(pp)
dev.off()

jpeg(filename="raw_density.jpeg", height=779, width=1400, quality=100)
res <- densityplot(RGList, channel="G", group="probe.type",
                   subset = c("Human", "Other", "Spike"))
print(res)
dev.off()

jpeg(filename="raw_log2_boxplot.jpeg", height=779, width=1400, quality=100)
boxplot(log2(cbind(RGList$G, RGList$R)), cex=.2,
        main="Raw log2 Red & Green Channel Intensities", xaxt="n",
        col=rep(c("green","red"),each=length(treatment)))
dev.off()


# densityplot(RGList, channel="G")
# densityplot(RGList, channel="R")

# plot_profile(log2(RGList$G), treatment)
# plot_profile(RGList$R, treatment)

plotPCA(log2(RGList$G), groups=as.numeric(as.factor(treatment)), 
        groupnames=levels(as.factor(treatment)), addtext=treatment, pcs=c(1,2,3), 
        plot3d=TRUE)

plotPCA(log2(RGList$G), groups=as.numeric(as.factor(treatment)), 
        groupnames=levels(as.factor(treatment)), addtext=colnames(RGList), 
        pcs=c(1,2), plot3d=FALSE, main="")

##'-----------------------------------------------------------------------------------------#

##'OUTLIER DETECTION------------------------------------------------------------------------#

outliers <- checkOutliers(RGList)
# boxplot(as.data.frame(t(RGList$R[outliers$Rout,])))

RGList$R <- fixOutliers(RGList$R, outliers$Rout, RGList$genes$Name)
RGList$G <- fixOutliers(RGList$G, outliers$Gout, RGList$genes$Name)

mvs <- checkMVs(RGList)
RGList$Rb <- fixMVs(RGList$Rb, mvs$Rb.na, RGList$genes$Name)
RGList$Gb <- fixMVs(RGList$Gb, mvs$Gb.na, RGList$genes$Name)

jpeg(filename="profile_plot_log2G_PostProbeCorrection.jpeg", height=779, width=1400, quality=100)
tmp <- RGList
colnames(tmp) <- short_colnames
pp <- plot_profile(log2(tmp$G))
print(pp)
dev.off()

##'-----------------------------------------------------------------------------------------#

##'FILTERING--------------------------------------------------------------------------------#
##'NOTE: RGList$genes$Gene ID is a required field.
##'RGList$genes$'Gene ID' <- RGList$genes$ID

reducedSet <- filterArray(RGList,  keep=probes_in, 
                          frac=1.2, number=floor(ncol(RGList$G)*.2), reps=4)



# reducedSet.spike <- filterArray(RGList,  keep=probes_in[-3], 
#                           frac=1.2, number=floor(ncol(RGList$G)*.2), reps=4)

jpeg(filename="profile_plot_log2G_PreNorm.jpeg", height=779, width=1400, quality=100)
tmp <- reducedSet
colnames(tmp) <- short_colnames
pp <- plot_profile(log2(tmp$G))
print(pp)
dev.off()

# plot_profile(reducedSet$R, treatment)
# 
# plotPCA(log2(reducedSet$G), groups=as.numeric(as.factor(treatment)), 
#         groupnames=levels(as.factor(treatment)), addtext=treatment, pcs=c(1,2,3), 
#         plot3d=TRUE)
# 
# plotPCA(log2(reducedSet$G), groups=as.numeric(as.factor(treatment)), 
#         groupnames=levels(as.factor(treatment)), addtext=colnames(RGList), pcs=c(1,2), 
#         plot3d=FALSE, main="")

##'-----------------------------------------------------------------------------------------#

##'NORMALISATION----------------------------------------------------------------------------#

##'None
ndata.none        <- normalizeWithinArrays(reducedSet, method="none")
##'Median
ndata.median      <- normalizeWithinArrays(reducedSet, method="median")
##'Loess
ndata.loess       <- normalizeWithinArrays(reducedSet, method="loess")
##'Quantile
ndata.quantile    <- normalizeBetweenArrays(reducedSet, method="quantile")
##'VSN
ndata.vsn.limma   <- normalizeVSN(reducedSet)
##'Spike-in VSN
idx.control       <- grep('spike', reducedSet$genes$Name)
spikein.fit       <- vsn2(reducedSet[idx.control,], lts.quantile=1, 
                          backgroundsubtract=TRUE)
ndata.spikein.vsn <- predict(spikein.fit, newdata=reducedSet)
spike             <- new("RGList", as.list(assayData(ndata.spikein.vsn)))
spike$genes       <- reducedSet$genes 
spike.norm        <- normalizeWithinArrays(spike, method="none")


ndata.all         <- list(ndata.none, ndata.median, ndata.loess,
                          ndata.quantile, ndata.vsn.limma, spike.norm)
#                   ndata.spikein.vsn)
names(ndata.all)  <- c("None", "Median", "Loess", "Quantile", "VSN",
                      "Spike-in VSN")

##'SPIKE-IN NORMALISATION TO MAList
for(i in 1:length(names(ndata.all))) {
  dp <- densityplot(list(ndata.all[[i]]), channel="G", group="probe.type",
              subset = c("Human", "Other", "Spike"))
  print(dp)
  cat(names(ndata.all)[i])
  cat ("Press [enter] to continue")
  line <- readline()
}
##'-----------------------------------------------------------------------------------------#

##'NORMALISATION VISUALISATION--------------------------------------------------------------#
dplot <- densityplot(ndata.all, channel="G", group="probe.type",
                     subset = c("Human", "Other", "Spike"),
                     par.strip.text=list(cex=0.75))

res <- MADvsMedianPlot(ndata.all, channel="R", group="probe.type",
                       subset=c("Human", "Other", "Spike"))

lp <- levelplot(ndata.all, channel="G", order=c(1:36),
                scales = list(rot=c(45, 45)))


##'-----------------------------------------------------------------------------------------#

##'PCA NORMALISED REFERENCE CHANEL----------------------------------------------------------#

pca_mirna <- function(x, file_name) {
  x <- RG.MA(x)
  y <- x$G
  rownames(y) <- x$genes$Name
  std_dev <- apply(y,1,sd)
  std_dev <- sort(std_dev, TRUE)
  std_dev <- as.matrix(std_dev)
  std_dev <- as.matrix(std_dev[grep("HIDDEN|Empty|spike|Hy3", rownames(std_dev), invert=T),])
  std_dev <- as.matrix(std_dev[grep('hsa', rownames(std_dev)),])
  filter  <- unique(rownames(std_dev))[1:50]
  final   <- y[rownames(y) %in% filter,]
  final   <- final[!duplicated(rownames(final)), ]

  split   <- strsplit(colnames(x),'_')
  colnames(final) <- unique(grep('143', unlist(split), value=T))
  
#   jpeg(filename="normalised_profile.jpeg", height=779, width=1400, quality=100)
  jpeg(filename="profile.jpeg", height=779, width=1400, quality=100)  
  g <- plot_profile(log2(na.omit(final)), sep=F)
  print(g)
  dev.off()

#   jpeg(filename="normalised_pca.jpeg", height=779, width=1400, quality=100)
  jpeg(filename="pca.jpeg", height=389, width=700, quality=100)
  plotPCA(log2(na.omit(final)), groups=as.numeric(as.factor(treatment)), 
          groupnames=levels(as.factor(treatment)), addtext=colnames(final), pcs=c(1,2), 
          plot3d=F, main="")
  dev.off()

#   jpeg(filename="normalised_heatmap.jpeg", height=779, width=1400, quality=100)
  jpeg(filename="heatmap.jpeg", height=389, width=700, quality=100)
  heatmap.2(log2(na.omit(final)), scale="row", cexRow=0.5, col=redgreen(75), key=TRUE, 
            symkey=FALSE, density.info="none", trace="none")
  dev.off()
  
#   jpeg(filename="normalised_density.jpeg", height=779, width=1400, quality=100)
  jpeg(filename="density.jpeg", height=779, width=1400, quality=100)
#   x <- x[x$genes$Name %in% filter,]
  den <- densityplot(x, channel="G", group="probe.type",
              subset = c("Human", "Other", "Spike"),
              par.strip.text=list(cex=0.75))
  print(den)
  dev.off()

  system("convert +append pca.jpeg heatmap.jpeg merge.jpeg")
  system("convert -append merge.jpeg profile.jpeg merge1.jpeg")
  system(paste("convert -append merge1.jpeg density.jpeg ", file_name, ".jpeg", sep=""))
  system("rm pca.jpeg profile.jpeg merge.jpeg heatmap.jpeg density.jpeg")
}

##'None
pca_mirna(ndata.all[[1]], "normalisation_none")
##'Median
pca_mirna(ndata.all[[2]], "normalisation_median")
##'Loess
pca_mirna(ndata.all[[3]], "normalisation_loess")
##'Quantile
pca_mirna(ndata.all[[4]], "normalisation_quantile")
 ##'VSN
pca_mirna(ndata.all[[5]], "normalisation_vsn")
##'Spike-in VSN
pca_mirna(spike.norm, "normalisation_spike_filtered")

plot_profile(log2(na.omit(spike$G)), treatment, F)

densityplot(list(spike.norm), channel="G", group="probe.type",
            subset = c("Human", "Other", "Spike"))

##'-----------------------------------------------------------------------------------------#

##'LIMMA------------------------------------------------------------------------------------#
##'
treatments       <- c("Co", "Ly", "Hi")#, "Lo")#
# treatment_s        <- treatment
# treatment_s[grep('Hi|Lo', treatment)] <- 'Sj'
array_names      <- treatment
# treatments       <- c("Co", "Ly", "Sj")
# array_names      <- treatment_s
design           <- model.matrix(~0 + factor(array_names, levels = treatments))
colnames(design) <- treatments
num_parameters   <- ncol(design)

cont_mat <- makeContrasts(Co-Ly, Co-Hi, Hi-Ly, levels=treatments)
# cont_mat <- makeContrasts(Co-Sj, Ly-Sj,  levels=treatments)
# cont_mat <- makeContrasts(Co-Ly, Co-Hi, Co-Lo, Hi-Ly, Ly-Lo, Hi-Lo, levels=treatments)

## Order data by probes
idx      <- order(ndata.quantile$genes$Name)
ndata    <- ndata.quantile[idx,]
idx.rm   <- which(ndata$genes$probe.type=="Spike")
ndata    <- ndata[-idx.rm,]
## compute correlations between same probes on each chip
corfit   <- duplicateCorrelation(ndata, design, ndups=4,
                                    weights=ndata$weights)

fit      <- lmFit(ndata, design, ndups=4, correlation=corfit$consensus)
fit2     <- contrasts.fit(fit, contrasts=cont_mat)
fit2     <- eBayes(fit2)

##'-----------------------------------------------------------------------------------------#

##'DIFFERENTIAL EXPRESSION------------------------------------------------------------------#
##'

# comparisons <- c("Co - Ly", "Co - Hi", "Co - Lo", "Hi - Ly", "Hi - Lo", "Ly - Lo")
comparisons <- c("Co - Ly", "Co - Hi", "Hi - Ly")
# comparisons <- c("Co - Sj", "Ly - Sj")
p_cut_off   <- 0.05
fold_change <- 2
i           <- 3

data_in <- list()

# for(i in 1:length(comparisons)) {
  gene_list            <- topTable(fit2, coef = comparisons[i], p.value = p_cut_off, lfc = log2(fold_change),
                                   nrow(fit2$p.value), adjust.method = "BH")
  
  write(gene_list$Name, file=paste(comparisons[i], ".txt"))
  gene_list_unfiltered <- topTable(fit2, coef = comparisons[i], number = nrow(fit2$p.value), adjust.method="BH")
  data_in[[3]] <- gene_list_unfiltered
  cat(comparisons[i])
  print(length(gene_list$logFC))
  
#   gene_list_unfiltered$threshold = as.factor(abs(gene_list_unfiltered$logFC) > log2(fold_change) &
#                                                gene_list_unfiltered$adj.P.Val < p_cut_off)
#   g = ggplot(data=gene_list_unfiltered, aes(x=gene_list_unfiltered$logFC,
#                                             y=-log10(gene_list_unfiltered$adj.P.Val),
#                                             colour=gene_list_unfiltered$threshold)) +
#     geom_point(alpha=0.4, size=1.75) +
#     theme(legend.position = "none") +
#     #xlim(c(-4, 4)) + ylim(c(0, 4.5)) +
#     xlab("log2 fold change") + ylab("-log10 Adj.P-value") +
#     ggtitle(paste("Volcano Plot of Fold Change vs Adj.P-Value", comparisons[i], sep=" "))
#   print(g)

  print(gg_volcano(gene_list_unfiltered))
  
  cat ("Press [enter] to continue")
  line <- readline()
}
##'-----------------------------------------------------------------------------------------#

res <- decideTests(fit2)
vennDiagram(res, include=c("up","down"),
               counts.col=c("red","green"), cex=1.25)




if (require("multtest")) {
  rawp <- fit2$F.p.value
  BH <- mt.rawp2adjp(rawp, proc = "BH")
  sum(BH$adjp[,"BH"] < 0.05)
}
contr.helmert(3)

contrast.helmert <- cont_mat
fitc.helmert <- contrasts.fit(fit, contrast.helmert)
fitc.helmert <- eBayes(fitc.helmert)
Fstats <- topTable(fit2, coef="Co - Ly", number=nrow(ndata)/4, adjust="BH", sort.by="P")
sum(Fstats$adj.P.Val < 0.05)



## average probes for plotting
avedata <- avedups(ndata, ndups=4, spacing=1)
sigFgenes <- Fstats$Gene.ID[which(Fstats$adj.P.Val < 0.05)]
## ssgenes <- unique(c(sig13v12$Gene, sig14v12$Gene, sig14v13$Gene))
## sum(sigFgenes%in%ssgenes) ## all 79 show up in individual lists
mat <- as.matrix(avedata[match(sigFgenes, avedata$genes$Gene), ])
colnames(mat) <- treatment
rownames(mat) <- sigFgenes


treatment_num <- treatment
treatment_num[treatment_num == "Co"] <- 1
treatment_num[treatment_num == "Ly"] <- 2
treatment_num[treatment_num == "Hi"] <- 3
# treatment_num[treatment_num == "Lo"] <- 4
require("clValid")
## This creates averages over replicates for each day
aveExpr <- t(apply(mat, 1, function(x) tapply(x, treatment_num, mean)))
clRes <- clValid(aveExpr, 6, clMethod=c("hierarchical", "diana", "sota","kmeans"),
                 validation=c("internal"), metric="correlation")
summary(clRes)


exists("clRes")
clusters <- cutree(clRes@clusterObjs$hierarchical, 6)
## Scales the average expression values
aveExpr <- t(scale(t(aveExpr)))
colnames(aveExpr) <- unique(treatment)
jpeg(filename="clustering_cont1.jpeg", height=779, width=1400, quality=100)
clustPlot(clusters, aveExpr, 3, 2)
dev.off()

miRNames <- function(ids, geneNames, geneIDs) {
  ids.names <- geneNames[which(geneIDs%in%ids)]
  ids.mmu.names <- lapply(strsplit(ids.names, " "), function(x) grep("hsa", x, value=TRUE))
  ids.mmu.names <- unlist(lapply(ids.mmu.names, function(x) gsub("\\;", "", x)))
  ## Convert to standard nomenclature ...
  miRs <- tolower(ids.mmu.names)
  miRs <<- gsub("mir", "miR", miRs)
} 


paste(miRs1, collapse=", ")


## Pull out the hsa specific names
ids1 <- names(clusters[which(clusters==2)])
miRs1 <- miRNames(ids1, avedata$genes$Name, avedata$genes$"Gene ID")

require("targetscan.Hs.eg.db")
## These are "miRNAAnnDbBimap" - Bi-mappings
res01 <- miRs1%in%ls(targetscan.Hs.egMIRNA) ## only two false
miRs1[!res01] ## "hsa-miR-5684"
miRs1 <- miRs1[res01]
# miRs1 <- miRs1[-15]
miRs1 <- unique(miRs1)

miRs1.list <- mget(miRs1, targetscan.Hs.egMIRNA)
miRs1.fams <- mget(miRs1, targetscan.Hs.egMIRBASE2FAMILY)
miRs1.targets <- mget(as.character(miRs1.fams), revmap(targetscan.Hs.egTARGETS))
targets.tscan <- unique(unlist(miRs1.targets))
length(targets.tscan) #321 Targets | 285 Targets


require("microRNA")
## Try identifying targets for miRNAs in cluster one ...
data(hsTargets)
targets.miRB <- hsTargets$target[which(hsTargets$name%in%miRs1)]
targets.miRB <- unique(targets.miRB)
length(targets.miRB)
## target ensembl IDs
## head(targets.miRB) ## Ensembl IDs
## head(targets.tscan) ## Entrez gene idenifiers
## Convert to common naming for use with GSEA, use Entrez Gene IDs

require("org.Hs.eg.db")
## org.Mm.egENSEMBLTRANS Map Ensembl transcript acession numbers with
## Entrez Gene identifiers
idx.miRB <- as.character(targets.miRB)%in%ls(revmap(org.Hs.egENSEMBLTRANS))
targets.miRB.list <- as.character(targets.miRB)[idx.miRB]
targets.miRB.entrez <- unlist(mget(targets.miRB.list, revmap(org.Hs.egENSEMBLTRANS)))
targets.intsect <- intersect(targets.tscan, targets.miRB.entrez)
length(targets.intsect)







