#-------------------------------------------------------------------------------------------#
#  Author      : Andrew J Skelton                                                           |
#  Language    : R Statistical Programming Language                                         |
#  Study       : Methylation of Cartilage Study - Hip / Knee OA                             |
#  Data Owner  : Newcastle University - Prof. John Loughlin, Dr. Louise Reynard             |
#  Description : Combine 3 Illumina Infinium 450K Methylations Array Studies to increase    |
#                Overall Experimental Power.                                                |
#-------------------------------------------------------------------------------------------#



##'Set the Working Directory, load essential packages, and import files
##'-----------------------------------------------------------------------------------------#
setwd("~/Loughlin/project_master_branch/methylation/")

library(methyAnalysis)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(FDb.InfiniumMethylation.hg19)
# library(IlluminaHumanMethylation450k.db)
library(BSgenome.Hsapiens.UCSC.hg19)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(lumi)
library(limma)
library(ggplot2)
library(reshape2)
library(scales)
library(sva)
library(bumphunter)
library(rtracklayer)
library(GenomicFeatures)
library(readr)

pheno_in <- read.table("pheno/all_pheno_4.txt",
                       row.names=1,
                       header=T,
                       sep="\t",
                       stringsAsFactors=F)

pheno_ex <- pheno_in[ (pheno_in$KL_Score == "-") |
                      (pheno_in$Discrepancy_KL == 1) |
                      (pheno_in$Discrepancy_Other == 1) |
                      (pheno_in$Discrepancy_Gender == 1) |
                      (pheno_in$Gender == "ND") |
                      (pheno_in$QC_Exclude == 1), ]

pheno_in <- pheno_in[ (pheno_in$KL_Score != "-") &
                      (pheno_in$Discrepancy_KL != 1) &
                      (pheno_in$Discrepancy_Other != 1) &
                      (pheno_in$Discrepancy_Gender != 1) &
                      (pheno_in$Gender != "ND") &
                      (pheno_in$QC_Exclude != 1), ]

# foo      <- read_delim("Louise_Samples_in.txt", 
#                        delim="\t", 
#                        col_names=F)
# pheno_in <- pheno_in[pheno_in$Sample.ID %in% foo[,1],]
# pheno_in <- pheno_in[!(pheno_in$between_experiment_replicate == 1 & 
#                          pheno_in$QC_Exclude == 1),]
##'-----------------------------------------------------------------------------------------#



##'Read in Raw IDAT Files
##'-----------------------------------------------------------------------------------------#
raw.idat <- importMethyIDAT(pheno_in,
                           lib="FDb.InfiniumMethylation.hg19",
                           dataPath="raw_data/All/")
##'-----------------------------------------------------------------------------------------#



##'Create Annotation Data Structures
##'-----------------------------------------------------------------------------------------#
annotation_env      <- IlluminaHumanMethylation450kanno.ilmn12.hg19@data
annotation_snp      <- as.data.frame(annotation_env$SNPs.Illumina)
annotation_snp      <- annotation_snp[annotation_snp$Probe_SNPs != "" |
                                      annotation_snp$Probe_SNPs_10 != "",]
annotation_island   <- as.data.frame(annotation_env$Islands.UCSC)
annotation_location <- as.data.frame(annotation_env$Locations)
annotation_manifest <- as.data.frame(annotation_env$Manifest)
annotation_other    <- as.data.frame(annotation_env$Other)
annotation_snp132   <- as.data.frame(annotation_env$SNPs.132CommonSingle)
annotation_snp135   <- as.data.frame(annotation_env$SNPs.135CommonSingle)
annotation_snp137   <- as.data.frame(annotation_env$SNPs.137CommonSingle)
annotation_snpIll   <- as.data.frame(annotation_env$SNPs.Illumina)
annotation_snpIll   <- annotation_snpIll[annotation_snpIll$Probe_SNPs != "" |
                                         annotation_snpIll$Probe_SNPs_10 != "",]
annotation_snp132   <- annotation_snp132[!is.na(annotation_snp132$Probe_rs),]
hm450.hg19          <- as.data.frame(getPlatform(platform='HM450', genome='hg19'))
##'-----------------------------------------------------------------------------------------#



##'Identify Probes for Removal
##'-----------------------------------------------------------------------------------------#
cpg_remove          <- rownames(annotation_snpIll)
# cpg_remove          <- rownames(annotation_snp132)

cpg_remove          <- c(cpg_remove, rownames(hm450.hg19[hm450.hg19$seqnames == "chrX" |
                                                         hm450.hg19$seqnames == "chrY",]))
cpg_remove          <- c(cpg_remove, rownames(raw.idat)[grep("rs", rownames(raw.idat))])

raw.idat.in         <- raw.idat[-match(cpg_remove, rownames(raw.idat)),]
##'-----------------------------------------------------------------------------------------#



##'Normalise Data -
##'Background Correct (Separate Colours)
##'Colour Bias Correction
##'Quantile Normalise
##'-----------------------------------------------------------------------------------------#
lumiMethy.b.adj     <- lumiMethyB(raw.idat.in,
                                  method="bgAdjust2C",
                                  separateColor=TRUE)

lumiMethy.c.adj     <- lumiMethyC(lumiMethy.b.adj,
                                  method="quantile")

lumiMethy.c.n       <- lumiMethyN(lumiMethy.c.adj,
                                  method="quantile",
                                  separateColor=TRUE)
##'-----------------------------------------------------------------------------------------#



##'Probe Level Filtering
##'Detection P-Value
##'-----------------------------------------------------------------------------------------#
exprs_data          <- exprs(lumiMethy.c.n)
present_count       <- detectionCall(lumiMethy.c.n)
normalised_data     <- exprs_data[present_count > 0, ]
##'-----------------------------------------------------------------------------------------#



##'Technical Effect Correction
##'-----------------------------------------------------------------------------------------#
#Separate Latest set of arrays to account for technical effect
pData(lumiMethy.c.n)[pData(lumiMethy.c.n)$SENTRIX_BARCODE == 9985131037,]$BATCH <- 4

batches             <- pData(lumiMethy.c.n)$BATCH
pheno               <- data.frame(sample=c(1:ncol(normalised_data)),
                                  outcome=pData(lumiMethy.c.n)$SAMPLETYPE,
                                  batch=batches)
rownames(pheno)     <- colnames(normalised_data)
mod                 <- model.matrix(~as.factor(outcome), data=pheno)
batchCorrected_data <- ComBat(dat=normalised_data,
                              batch=batches,
                              mod=mod,
                              par.prior=T,
                              prior.plots=F)

lumi.norm           <- lumiMethy.c.n
lumi.norm           <- lumi.norm[present_count > 0,]
exprs(lumi.norm)    <- batchCorrected_data
##'-----------------------------------------------------------------------------------------#



##'Utility - Extract Data about CpGs
##'-----------------------------------------------------------------------------------------#
lookup <- c("cg02963878")
foo    <- m2beta(exprs(lumi.norm)[match(lookup,rownames(lumi.norm)),])
anno   <- cbind(annotation_island[match(lookup, rownames(annotation_island)),],
                annotation_location[match(lookup, rownames(annotation_location)),],
                annotation_other[match(lookup, rownames(annotation_other)),],
                annotation_manifest[match(lookup, rownames(annotation_manifest)),])
write.csv(foo, file="CpG_Query.csv")
write.csv(anno, file="CpG_Query_Annotation.csv")
write.csv(pData(lumi.norm), file="CpG_Query_Pheno.csv")
##'-----------------------------------------------------------------------------------------#



##'Differential Methylation - MethyAnalysis Methodology - M Values Only
##'-----------------------------------------------------------------------------------------#
lumi.norm.in           <- lumi.norm[, grep("Hip_Control|OA_Hip", pData(lumi.norm)$SAMPLETYPE)]
sampleType             <- pData(lumi.norm.in)$SAMPLETYPE
allResult              <- detectDMR.slideWin(lumi.norm.in, sampleType=sampleType)
allDMRInfo             <- identifySigDMR(allResult, p.adjust.method="fdr",
                                         pValueTh=0.05, fdrTh=0.05)
DMRInfo.ann            <- annotateDMRInfo(allDMRInfo, 'TxDb.Hsapiens.UCSC.hg19.knownGene')
DMRInfo.ann.df.R       <- as.data.frame(DMRInfo.ann$sigDMRInfo)
DMRInfo.ann.df.P       <- as.data.frame(DMRInfo.ann$sigDataInfo)

export.DMRInfo(DMRInfo.ann, lumi.norm.in, savePrefix='Hip_OA_Vs_Hip_Control_B')
head(DMRInfo.ann)

foo                            <- as.data.frame(DMRInfo.ann$sigDataInfo)
foo_expr                       <- as.data.frame(exprs(lumi.norm.in)[foo$PROBEID,])
foo_expr                       <- m2beta(foo_expr)
Beta_Mean_Hip_Control          <- as.numeric(rowMeans(foo_expr[,grep("Hip_Control", sampleType)]))
Beta_Mean_OA_Hip               <- as.numeric(rowMeans(foo_expr[,grep("OA_Hip", sampleType)]))
Beta_diff                      <- Beta_Mean_Hip_Control - Beta_Mean_OA_Hip
out                            <- data.frame(probe_id = rownames(foo),
                                             Beta_Mean_Hip_Control = Beta_Mean_Hip_Control,
                                             Beta_Mean_OA_Hip = Beta_Mean_OA_Hip,
                                             Beta_diff = Beta_diff)
df_out                         <- merge(foo, out, by.x="PROBEID", by.y="probe_id")

write.csv(df_out, file="DM_Hip_Control-OA_Hip.csv")
write.csv(df_out[,c(1, 19, 15:16, 7, 23:25)], file="DM_Hip_Control-OA_Hip_Streamline.csv")

foo <- MethyLumiM2GenoSet(lumi.norm.in,
                          lib="FDb.InfiniumMethylation.hg19")
export.methyGenoSet(foo,
                    file.format="bw",
                    exportValue="beta",
                    savePrefix="Hip_OA_Vs_Hip_Control")

pheno_hm           <- as.data.frame(pData(foo)[,6])
colnames(pheno_hm) <- "SampleType"
rownames(pheno_hm) <- rownames(pData(foo))

heatmapByChromosome(foo,
                    gene='8111',
                    genomicFeature='TxDb.Hsapiens.UCSC.hg19.knownGene',
                    includeGeneBody=TRUE,
                    phenoData=pheno_hm)
dev.off()
##'-----------------------------------------------------------------------------------------#



##'Differential Methylation - bumphunter
##'-----------------------------------------------------------------------------------------#
library(doParallel)
registerDoParallel(cores = 20)

hm450.hg19 <- getPlatform(platform='HM450', genome='hg19')
range_test <- hm450.hg19[rownames(lumi.norm.in)]
gr         <- sortSeqlevels(range_test)
gr         <- sort(gr)
range_test <- gr

chr        <- as.character(seqnames(range_test))
pos        <- start(range_test)

clusters   <- bumphunter::clusterMaker(chr,
                                       pos,
                                       assumeSorted=T,
                                       maxGap=300)
design.fu  <- model.matrix(~ pData(lumi.norm.in)$SAMPLETYPE)
bh_m       <- bumphunter::bumphunter(exprs(lumi.norm.in),
                                     design.fu,
                                     chr,
                                     pos,
                                     clusters,
                                     coef=2,
                                     B=100,
                                     verbose=TRUE,
                                     pickCutoff=TRUE)

cpg.cur <- read.delim("http://rafalab.jhsph.edu/CGI/model-based-cpg-islands-hg19.txt",
                      as.is=TRUE)

plotRegions(data.frame(chr="chr2", start=71785254,end=71795254), m2beta(exprs(lumi.norm.in)), chr, pos,
            exposure=pData(lumi.norm.in)$SAMPLETYPE, Genome=BSgenome.Hsapiens.UCSC.hg19,
            cpg.islands=cpg.cur , outfile='dmrTest.pdf')

# foo        <- as.data.frame(range_test)
# foo <- foo[foo$seqnames == "chr2" &
#            foo$start > 71770200 &
#            foo$end < 71792300,]
# range_test <- makeGRangesFromDataFrame(foo)
# mcols(range_test) <- foo[,6:ncol(foo)]
##'-----------------------------------------------------------------------------------------#



##'Differential Methylation - Limma - Grouped Design
##'-----------------------------------------------------------------------------------------#
treatments           <- unique(pData(lumi.norm)$SAMPLETYPE)
treatment_arrays     <- pData(lumi.norm)$SAMPLETYPE
design               <- model.matrix(~0 + factor(treatment_arrays,
                                                 levels=treatments))
colnames(design)     <- treatments
fit                  <- lmFit(m2beta(exprs(lumi.norm)), design)

cont_mat             <- makeContrasts(OA_Hip-Hip_Control,
                                      OA_Knee-Hip_Control,
                                      OA_Knee-OA_Hip,
                                      levels=treatments)
fit2                 <- contrasts.fit(fit,
                                      contrasts=cont_mat)
fit2                 <- eBayes(fit2)

contrast             <- c("OA_Hip - Hip_Control")
# contrast             <- c("OA_Knee - OA_Hip")
# contrast             <- c("OA_Knee - Hip_Control")

gene_list_unfiltered <- topTable(fit2,
                                 coef=contrast,
                                 number=Inf,
                                 adjust.method="BH",
                                 sort.by="P")
gene_list            <- topTable(fit2,
                                 coef=contrast,
                                 p.value=0.05,
                                 lfc=0.1,
                                 number=Inf,
                                 adjust.method="BH",
                                 sort.by="P")
nrow(gene_list)

file_name <- "OA_Knee_Vs_Hip_Control_P0.05_DB0.1"
lookup <- "cg00071274"
foo    <- m2beta(exprs(lumi.norm)[match(lookup,rownames(lumi.norm)),])
anno   <- cbind(annotation_island[match(lookup, rownames(annotation_island)),],
                annotation_location[match(lookup, rownames(annotation_location)),],
                annotation_other[match(lookup, rownames(annotation_other)),],
                annotation_manifest[match(lookup, rownames(annotation_manifest)),])
write.csv(foo, file=paste0(file_name, "_Betas.csv"))
write.csv(anno, file=paste0(file_name, "_Annotation.csv"))
write.csv(gene_list, file=paste0(file_name, "_DEData.csv"))

write.csv(pData(lumi.norm), file=paste0(file_name, "_Betas.csv"))
##'-----------------------------------------------------------------------------------------#


##'Differential Methylation - Limma - Regression by Phenotype - Age
##'-----------------------------------------------------------------------------------------#
#OA Hip
design       <- model.matrix(~pData(lumi.norm)[pData(lumi.norm)$SAMPLETYPE == "OA_Hip",]$AGE)
fit          <- lmFit((exprs(lumi.norm)[, grep("OA_Hip", pData(lumi.norm)$SAMPLETYPE)]), design)
fit_OA_Hip   <- eBayes(fit)
gene_list_OA_Hip   <- topTable(fit_OA_Hip,
                               coef=2,
                               p.value=0.001,
                               # lfc=0.03,
                               number=Inf,
                               adjust.method="none",
                               sort.by="P")

#OA Knee
design       <- model.matrix(~pData(lumi.norm)[pData(lumi.norm)$SAMPLETYPE == "OA_Knee",]$AGE)
fit          <- lmFit((exprs(lumi.norm)[, grep("OA_Knee", pData(lumi.norm)$SAMPLETYPE)]), design)
fit_OA_Knee  <- eBayes(fit)
gene_list_OA_Knee  <- topTable(fit_OA_Knee,
                               coef=2,
                               p.value=0.001,
                               # lfc=0.03,
                               number=Inf,
                               adjust.method="none",
                               sort.by="M")

#Hip Control
design       <- model.matrix(~pData(lumi.norm)[pData(lumi.norm)$SAMPLETYPE == "Hip_Control",]$AGE)
fit          <- lmFit((exprs(lumi.norm)[, grep("Hip_Control", pData(lumi.norm)$SAMPLETYPE)]), design)
fit_Hip_Ctrl <- eBayes(fit)
gene_list_Hip_Ctrl <- topTable(fit_Hip_Ctrl,
                               coef=2,
                               p.value=0.001,
                               # lfc=0.03,
                               number=Inf,
                               adjust.method="none",
                               sort.by="P")

file_name <- "Age_OA_Hip_P0.001"
lookup <- rownames(gene_list_OA_Hip)
foo    <- m2beta(exprs(lumi.norm)[match(lookup,rownames(lumi.norm)),])
anno   <- cbind(annotation_island[match(lookup, rownames(annotation_island)),],
                annotation_location[match(lookup, rownames(annotation_location)),],
                annotation_other[match(lookup, rownames(annotation_other)),],
                annotation_manifest[match(lookup, rownames(annotation_manifest)),])
write.csv(foo, file=paste0(file_name, "_Betas.csv"))
write.csv(anno, file=paste0(file_name, "_Annotation.csv"))
write.csv(gene_list_OA_Hip, file=paste0(file_name, "_DEData.csv"))
##'-----------------------------------------------------------------------------------------#



##'Differential Methylation - Limma - Identify CpGs with Opposite Slope Directions in 
##'                                   Different Tissues
##'Requires the Above Section to be run first... 
##'-----------------------------------------------------------------------------------------#
OA_Hip_in   <- data.frame(CpG=rownames(gene_list_OA_Hip),   logFC=gene_list_OA_Hip$logFC)
OA_Knee_in  <- data.frame(CpG=rownames(gene_list_OA_Knee),  logFC=gene_list_OA_Knee$logFC)
Hip_Ctrl_in <- data.frame(CpG=rownames(gene_list_Hip_Ctrl), logFC=gene_list_Hip_Ctrl$logFC)
  
length(intersect(OA_Knee_in$CpG, OA_Hip_in$CpG))  
df_compare_OA           <- merge(OA_Hip_in, OA_Knee_in, by="CpG")
colnames(df_compare_OA) <- c("CpG", "OA_Hip", "OA_Knee")
df_compare_OA[(df_compare_OA$OA_Hip > 0 & df_compare_OA$OA_Knee < 0) |
              (df_compare_OA$OA_Hip < 0 & df_compare_OA$OA_Knee > 0),]


length(intersect(Hip_Ctrl_in$CpG, OA_Hip_in$CpG)) 
df_compare_HipC_OAHip           <- merge(OA_Hip_in, Hip_Ctrl_in, by="CpG")
colnames(df_compare_HipC_OAHip) <- c("CpG", "OA_Hip", "Hip_Control")
df_compare_OA[(df_compare_HipC_OAHip$OA_Hip > 0 & df_compare_HipC_OAHip$Hip_Control < 0) |
                (df_compare_HipC_OAHip$OA_Hip < 0 & df_compare_HipC_OAHip$Hip_Control > 0),]

length(intersect(Hip_Ctrl_in$CpG, OA_Knee_in$CpG)) 
df_compare_HipC_OAKnee           <- merge(OA_Knee_in, Hip_Ctrl_in, by="CpG")
colnames(df_compare_HipC_OAKnee) <- c("CpG", "OA_Knee", "Hip_Control")
df_compare_OA[(df_compare_HipC_OAKnee$OA_Knee > 0 & df_compare_HipC_OAKnee$Hip_Control < 0) |
                (df_compare_HipC_OAKnee$OA_Knee < 0 & df_compare_HipC_OAKnee$Hip_Control > 0),]

##'-----------------------------------------------------------------------------------------#



##'Differential Methylation - Limma - Plot Regression
##'-----------------------------------------------------------------------------------------#
cpg_in  <- c("cg00071274")
cpg_in  <- rownames(gene_list_OA_Hip[1:5,])
df      <- melt(m2beta(exprs(lumi.norm)[cpg_in,]))

if(length(cpg_in) > 1) {
  colnames(df) <- c("CpG", "Sample", "Beta")
} else {
  df$CpG          <- cpg_in
  df$Sample       <- rownames(df)
  rownames(df)    <- c(1:nrow(df))
  df              <- df[,c(2,3,1)]
  colnames(df)[3] <- "Beta"
}

df$Age  <- pData(lumi.norm)[match(df$Sample, rownames(pData(lumi.norm))),]$AGE
df$Type <- pData(lumi.norm)[match(df$Sample, rownames(pData(lumi.norm))),]$SAMPLETYPE

gg      <- ggplot(df, aes(x=Age, y=Beta, colour=Type)) + 
             geom_point() +
             geom_smooth(method=lm, se=T, aes(group=1)) +
             # scale_y_continuous(limits=c(0,1)) +
             theme_bw()

if(length(cpg_in) > 1) {
  gg <- gg + 
        facet_grid(CpG ~ Type)
} else {
  gg <- gg + 
        facet_grid(. ~ Type) +
        ggtitle(cpg_in)
}
print(gg)
##'-----------------------------------------------------------------------------------------#



##'Plot Range
##'-----------------------------------------------------------------------------------------#

##'-----------------------------------------------------------------------------------------#



##'Differential Methylation - Limma - Grouped Design - Gender Stratification
##'-----------------------------------------------------------------------------------------#
treatments           <- unique(pData(lumi.norm)$GENDER)
treatment_arrays     <- pData(lumi.norm)[pData(lumi.norm)$SAMPLETYPE == "Hip_Control",]$GENDER
design               <- model.matrix(~0 + factor(treatment_arrays,
                                                 levels=treatments))
colnames(design)     <- treatments
fit                  <- lmFit(m2beta(exprs(lumi.norm[, pData(lumi.norm)$SAMPLETYPE == "Hip_Control"])), design)

cont_mat             <- makeContrasts(Male-Female,
                                      levels=treatments)
fit2                 <- contrasts.fit(fit,
                                      contrasts=cont_mat)
fit2                 <- eBayes(fit2)

contrast             <- c("Male - Female")

gene_list_unfiltered <- topTable(fit2,
                                 coef=contrast,
                                 number=Inf,
                                 adjust.method="BH",
                                 sort.by="P")
gene_list            <- topTable(fit2,
                                 coef=contrast,
                                 p.value=0.01,
                                 lfc=0.1,
                                 number=Inf,
                                 adjust.method="BH",
                                 sort.by="P")
nrow(gene_list)



file_name <- "Hip_Control_Male_Vs_Female_P0.01_DB0.1"
lookup <- rownames(gene_list)
foo    <- m2beta(exprs(lumi.norm)[match(lookup,rownames(lumi.norm)),])
anno   <- cbind(annotation_island[match(lookup, rownames(annotation_island)),],
                annotation_location[match(lookup, rownames(annotation_location)),],
                annotation_other[match(lookup, rownames(annotation_other)),],
                annotation_manifest[match(lookup, rownames(annotation_manifest)),])
write.csv(foo, file=paste0(file_name, "_Betas.csv"))
write.csv(anno, file=paste0(file_name, "_Annotation.csv"))
write.csv(gene_list, file=paste0(file_name, "_DEData.csv"))

write.csv(pData(lumi.norm), file=paste0(file_name, "_Betas.csv"))
##'-----------------------------------------------------------------------------------------#



##'Chromatin State
##'-----------------------------------------------------------------------------------------#
hm450.hg19          <- getPlatform(platform="HM450", genome="hg19")
InfiniumMethylation <- features(FDb.InfiniumMethylation.hg19)
session             <- browserSession()
genome(session)     <- "hg19"
# trackNames(session)
query               <- ucscTableQuery(session, "wgEncodeBroadHmm")

tableName(query)    <- "wgEncodeBroadHmmGm12878HMM"
Gm12878             <- track(query, asRangedData = FALSE)
tableName(query)    <- "wgEncodeBroadHmmH1hescHMM"
H1hesc              <- track(query, asRangedData = FALSE)

Gm12878_split       <- split(Gm12878, Gm12878$name)
H1hesc_split        <- split(H1hesc, H1hesc$name)

df_chrom            <- c()
for(i in 1:length(Gm12878_split)) {
  tmp      <- as.data.frame(subsetByOverlaps(InfiniumMethylation, Gm12878_split[[i]]))
  val      <- unique(elementMetadata(Gm12878_split[[i]])$name)
  tmp      <- data.frame(cpg_id=rownames(tmp), type=val)
  df_chrom <- rbind(df_chrom, tmp)
}
head(df_chrom)

df_chrom           <- df_chrom[-na.omit(match(cpg_remove, df_chrom$cpg_id)),]
matched_exprs      <- m2beta(na.omit(exprs(lumi.norm)[match(df_chrom$cpg_id, rownames(lumi.norm)),]))
matched_exprs_melt <- melt(matched_exprs)

foo           <- cbind(matched_exprs_melt, df_chrom[match(matched_exprs_melt$Var1, df_chrom$cpg_id),])
colnames(foo) <- c("CpG", "Sample", "Beta", "CpG_Map", "Chromatin_State")
foo$Type      <- pData(lumi.norm)[match(foo$Sample, rownames(pData(lumi.norm))),]$SAMPLETYPE

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", 
                "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
                "#999999", "#E69F00", "#56B4E9", "#009E73", 
                "#F0E442", "#0072B2", "#D55E00")

png("Chromatin_State.png", width=4096, height=6096, units="px", res=300)
ggplot(foo, aes(x=Chromatin_State, y=Beta, fill=Chromatin_State)) +
  geom_violin() +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90, hjust=1)) +
  scale_fill_manual(values=cbbPalette) +
  facet_grid(Type ~ ., scales="free_x")
dev.off()
#'-----------------------------------------------------------------------------------------#



