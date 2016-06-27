#-------------------------------------------------------------------------------------------#
#  Author      : Andrew J Skelton                                                           |
#  Language    : R Statistical Programming Language                                         |
#  Study       : Dr. Louise Reynard - MSc Project                                           |
#  Data Owner  : Newcastle University - Prof. Fai Ng, Shereen Al-Ali                        |
#  Description : 450K Visualisations                                                        |
#-------------------------------------------------------------------------------------------#



##'Set Directory and load packages
##'-----------------------------------------------------------------------------------------#
setwd("~/Documents/Bioinformatics/Customers/Louise/")

source('http://bioconductor.org/biocLite.R')
biocLite("IlluminaHumanMethylation450kanno.ilmn12.hg19")

library(readr)
library(ggplot2)
library(biomaRt)
library(gridExtra)
library(biovizBase)
library(GenomicFeatures)
library(ggbio)
library(FDb.InfiniumMethylation.hg19)
library(reshape2)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
##'-----------------------------------------------------------------------------------------#



##'Load Known Annotation (hg19)
##'-----------------------------------------------------------------------------------------#
hm450.hg19          <- as.data.frame(getPlatform(platform = 'HM450', 
                                                 genome   = 'hg19'))
annotation_env      <- IlluminaHumanMethylation450kanno.ilmn12.hg19@data
annotation_other    <- as.data.frame(annotation_env$Other)
##'-----------------------------------------------------------------------------------------#



##'Load Known GTF (hg19)
##'-----------------------------------------------------------------------------------------#
gtf_path     <- "./genome.gtf"
transcriptdb <- makeTxDbFromGFF(file     = gtf_path, 
                                format   = 'gtf', 
                                organism = 'Homo sapiens')
##'-----------------------------------------------------------------------------------------#



##'Load CpG Data (from pyrosequencing)
##'-----------------------------------------------------------------------------------------#
cpg_data_in  <- read_delim("blood.txt", "\t") 
# cpg_data_in  <- read_delim("synovium.txt", "\t") 
# cpg_data_in  <- read_delim("fat_pad.txt", "\t") 

foo <- c()
for(i in 6:12) {
  tmp              <- cpg_data_in[,c(1:5,i)]
  tmp$Probe        <- colnames(cpg_data_in)[i]
  tmp$Tissue       <- "Blood"
  colnames(tmp)[6] <- "Perc_Methyl" 
  foo <- rbind(foo, tmp)
}

# blood_data <- foo
# syno_data  <- foo
# fp_data    <- foo

df_in             <- na.omit(rbind(blood_data[,c(2,6:8)], 
                                   fp_data[,c(3,6:8)],
                                   syno_data[,c(2,5:7)]))
df_in             <- df_in[df_in$Perc_Methyl != "-",]
df_in$Perc_Methyl <- as.numeric(df_in$Perc_Methyl)
colnames(df_in)[1]<- "Genotype"
df_in             <- df_in[df_in$Probe != "mean",]
df_in$Genotype    <- factor(df_in$Genotype, 
                            levels=c("AA", "GG", "AG"))
df_in$Probe       <- factor(df_in$Probe, 
                            levels=c("cg13979708","cg19254793","cg20913747",
                                     "cg18551225","+7bp","+11bp"))

cpg_in     <- unique(grep("cg", 
                          blood_data$Probe, 
                          value=T))

cpg_in_an  <- hm450.hg19[cpg_in,]
annotation_other[cpg_in,]
##'-----------------------------------------------------------------------------------------#


##'Query GRCh37
##'-----------------------------------------------------------------------------------------#
ensembl    <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                      dataset = "hsapiens_gene_ensembl",
                      host    = "grch37.ensembl.org",
                      path    = "/biomart/martservice")
annotation <- getBM(attributes=c("ensembl_gene_id",
                                 "external_gene_name",
                                 "description"),
                    filters   =  "external_gene_name",
                    values    =  "RUNX2",
                    mart      = ensembl)
##'-----------------------------------------------------------------------------------------#

chr_subset    <- as.numeric(unique(cpg_in_an$seqnames))
start_subset  <- min(cpg_in_an$start)-1000
stop_subset   <- max(cpg_in_an$start)+1000000
trans         <- crunch(transcriptdb, 
                        which=GRanges(seqnames=chr_subset, 
                                      IRanges(start_subset, 
                                              stop_subset)))
gr1           <- split(trans, 
                       trans$tx_name)
p1            <- autoplot(gr1) + 
                  theme_bw() +
                  geom_vline(xintercept=c(min(cpg_in_an$start),
                                          max(cpg_in_an$start)), 
                             colour="red")
# p1            <- p1@ggplot
p.ideo        <- Ideogram(genome = "hg19", 
                          subchr=paste0("chr",
                                        chr_subset)) +
                 xlim(GRanges(paste0("chr",
                                     chr_subset), 
                              IRanges(start_subset, 
                                      stop_subset)))
# p.ideo        <- p.ideo@ggplot

# grid.arrange(arrangeGrob(p.ideo@ggplot,p1@ggplot), 
             # p.ideo@ggplot, ncol=1)
# ga <- grid.arrange(p.ideo@ggplot,p1@ggplot, heights=(c(2,6)), ncol=1)

gg <- ggplot(df_in, 
             aes(x=Probe, 
                 y=Perc_Methyl)) +
        geom_boxplot(aes(fill      = Genotype,
                         colour    = Genotype),
                     notch         = F,
                     outlier.shape = NA,
                     alpha         = 0.4) +
        geom_point(aes(shape       = Genotype, 
                       fill        = Genotype),
                   position        = position_dodge(width=.75), 
                   pch             = 21, 
                   size            = 2,
                   alpha           = 0.4) +
        stat_summary(fun.y         = median, 
                     geom          = "line",
                     position      = position_dodge(width=0.75),
                     aes(colour    = Genotype, 
                         group     = Genotype)) +
        theme_bw() + 
        xlab("") +
        ylab("% Methylation") +
        theme(axis.text.x = element_text(angle = -90, 
                                         hjust = 1)) +
        facet_grid(Tissue ~ .)

png("Louise_visualisation_4_mean.png", width=6096, height=6096, units="px", res=300)
grid.arrange(arrangeGrob(p.ideo@ggplot,
                          p1@ggplot, 
                          heights=(c(1,8)), 
                          ncol=1), 
             gg, ncol=2)
dev.off()

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
ensembl    <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                      dataset = "hsapiens_gene_ensembl",
                      host    = "grch37.ensembl.org",
                      path    = "/biomart/martservice")
annotation <- getBM(attributes=c("ensembl_gene_id",
                                 "external_gene_name",
                                 "description"),
                    filters   =  "external_gene_name",
                    values    =  "RUNX2",
                    mart      = ensembl)
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







