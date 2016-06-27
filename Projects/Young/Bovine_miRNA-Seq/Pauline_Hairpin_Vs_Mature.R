#-------------------------------------------------------------------------------------------#
#  Author      : Andrew J Skelton                                                           |
#  Language    : R Statistical Programming Language                                         |
#  Study       : miRNA Sequencing of Pooled Bovine Cartilage Samples                        |
#  Data Owner  : Newcastle University - Prof. David Young, (Cardiff Collaboration           |
#  Description : Mature Sequences Vs Hairpin Sequences                                      |
#-------------------------------------------------------------------------------------------#



##'Set the Working Directory, load essential packages, and import files
##'-----------------------------------------------------------------------------------------#
setwd("~/Paulina/Resequencing/April/Alignment/mature/cow/")

source('http://bioconductor.org/biocLite.R')
biocLite()

library(stringr)
library(Rsamtools)
library(ShortRead)
library(DESeq2)
library(ggplot2)
##'-----------------------------------------------------------------------------------------#



##'Mature - Load Files, Count, build Count Structure
##'-----------------------------------------------------------------------------------------#
setwd("~/Paulina/Resequencing/April/Alignment/mature/cow/")
matures         <- readDNAStringSet("../../../bowtie2Idx/DNA/mature/cow/cow_mature_dna.fa", 
                                    use.names=TRUE)
names(matures)  <- str_match(names(matures),
                             "^\\S+") 
matureGR        <- GRanges(names(matures), 
                           IRanges(1,width(matures)), 
                           strand="*")
bamView         <- BamViews(list.files("./", 
                                       pattern="*.bam$"), 
                            bamRanges=matureGR)
bamcounts       <- countBam(bamView)
foo             <- as.data.frame(bamcounts)

out             <- data.frame(miRBase_ID=unique(foo$space))
for(i in unique(foo$group_name)) {
  tmp           <- foo[foo$group_name == i,]
  tmp           <- tmp[,c(3,8)]
  out           <- merge(out, 
                         tmp, 
                         by.x="miRBase_ID", 
                         by.y="space")
  colnames(out)[ncol(out)] <- i
}
rownames(out)   <- out$miRBase_ID
out_mature      <- out[,-1]
##'-----------------------------------------------------------------------------------------#



##'Hairpin - Load Files, Count, build Count Structure
##'-----------------------------------------------------------------------------------------#
setwd("~/Paulina/Resequencing/April/Alignment/hairpin/cow/")
hairpins        <- readDNAStringSet("../../../bowtie2Idx/DNA/hairpin/cow/cow_hairpin_dna.fa", 
                                    use.names=TRUE)
names(hairpins) <- str_match(names(hairpins),"^\\S+") 
hairpinGR       <- GRanges(names(hairpins), 
                           IRanges(1,width(hairpins)), strand="*")
bamView         <- BamViews(list.files("./", 
                                       pattern="*.bam$"), 
                            bamRanges=hairpinGR)
bamcounts       <- countBam(bamView)
foo             <- as.data.frame(bamcounts)

out             <- data.frame(miRBase_ID=unique(foo$space))
for(i in unique(foo$group_name)) {
  tmp           <- foo[foo$group_name == i,]
  tmp           <- tmp[,c(3,8)]
  out           <- merge(out, 
                         tmp, 
                         by.x="miRBase_ID", 
                         by.y="space")
  colnames(out)[ncol(out)] <- i
}
rownames(out)   <- out$miRBase_ID
out_hairpin     <- out[,-1]
##'-----------------------------------------------------------------------------------------#



##'DESeq2 - Setup
##'-----------------------------------------------------------------------------------------#
setwd("~/Paulina/Resequencing/April/DESeq2/")
sampleTable             <- read.table("pheno.txt", 
                                      header=T, 
                                      row.names=1, 
                                      sep="\t", 
                                      stringsAsFactors=T)
colnames(out)           <- rownames(sampleTable)
sampleTable_bck         <- sampleTable
sampleTable$UniqueClass <- paste0(sampleTable$Pressure, 
                                  sampleTable$TimePoint)
sampleTable$TimePoint   <- factor(sampleTable$TimePoint, 
                                  levels=c("2hr", "6hr", "24hr"))
sampleTable$Prefix      <- factor(sampleTable$Prefix, 
                                  labels=c("A","B","C","D","E","F","G","H","I"))

colnames(out_hairpin)   <- rownames(sampleTable)
colnames(out_mature)    <- rownames(sampleTable)

time_in                 <- "6hr"

out_hairpin_in          <- out_hairpin[, grep(time_in, 
                                              sampleTable$TimePoint)]
out_mature_in           <- out_mature[, grep(time_in, 
                                             sampleTable$TimePoint)]
sampleTable_in          <- sampleTable[sampleTable$TimePoint == time_in,]
##'-----------------------------------------------------------------------------------------#



##'DESeq2 - Mature - Run
##'-----------------------------------------------------------------------------------------#
dds2_hairpin           <- DESeqDataSetFromMatrix(countData = out_hairpin_in,
                                                 colData   = sampleTable_in,
                                                 design    = ~ Prefix + Pressure)
dds2_hairpin           <- DESeq(dds2_hairpin)

deseq_rld               <- rlog(dds2_hairpin)
deseq_vst               <- varianceStabilizingTransformation(dds2_hairpin)
deseq_hairpin_normCoun  <- counts(dds2_hairpin, normalized=T)
deseq_hairpin_rawCoun   <- counts(dds2_hairpin, normalized=F)
##'-----------------------------------------------------------------------------------------#



##'DESeq2 - Hairpin - Run
##'-----------------------------------------------------------------------------------------#
dds2_mature            <- DESeqDataSetFromMatrix(countData = out_mature_in,
                                                 colData   = sampleTable_in,
                                                 design    = ~ Prefix + Pressure)
dds2_mature           <- DESeq(dds2_mature)

deseq_rld               <- rlog(dds2_mature)
deseq_vst               <- varianceStabilizingTransformation(dds2_mature)
deseq_mature_normCoun   <- counts(dds2_mature, normalized=T)
deseq_mature_rawCoun    <- counts(dds2_mature, normalized=F)
##'-----------------------------------------------------------------------------------------#



##'Plot Hairpin Vs Mature Comparison 
##'-----------------------------------------------------------------------------------------#
mirna          <- "bta-mir-21"
hairpin_counts <- deseq_hairpin_normCoun
mature_counts  <- deseq_mature_normCoun

# pdf("bta-mir-21.pdf", paper="a4r")
plot_paired_comp(mirna, 
                 deseq_hairpin_normCoun, 
                 deseq_mature_normCoun, 
                 sampleTable_in)
# dev.off()
##'-----------------------------------------------------------------------------------------#



