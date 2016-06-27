#-------------------------------------------------------------------------------------------#
#  Author      : Andrew J Skelton                                                           |
#  Language    : R                                                                          |
#  Study       : Quantification of SRX544934 from ENA                                       |
#  Data Owner  : ENA - At the request of Colin Shepard                                      |
#  Description : Visualise the counts of MCF2L        															        |
#-------------------------------------------------------------------------------------------#



##'Set the Working Directory, load essential packages, and import files
##'-----------------------------------------------------------------------------------------#
source("http://bioconductor.org/biocLite.R")
biocLite()

library(EBSeq)
library(ggplot2)
library(biomaRt)

setwd("~/Young/Pipeline_Analysis/salmon/grch38_gene/genes/")

pheno_in <- read.table("../../../scripts/Pheno.txt",
                       sep="\t",
                       header=T,
                       row.names=1)
##'-----------------------------------------------------------------------------------------#


##'Utility Function - get_counts
##'                 - Description : Reads in a Salmon Output file and converts it to a
##'                                 useful dataframe structure
##'                 - Parameters  : file - path and filename of the Salmon .sf file
##'2 - TPM
##'3 - FPKM
##'-----------------------------------------------------------------------------------------#
get_counts = function(file) {
  dat               <- read.table(file, skip=10, sep="\t", stringsAsFactors=F, row.names=1)
  raw_counts        <- dat[,3]
  names(raw_counts) <- rownames(dat)
  return(raw_counts)
}
##'-----------------------------------------------------------------------------------------#



##'Salmon - Gene Level - Get Normalised Counts
##'-----------------------------------------------------------------------------------------#
gene_files_in       <- list.files("./",
                                  pattern="*.sf",
                                  full.names=T)
gene_counts         <- lapply(gene_files_in,
                              get_counts)
gene_count_table    <- matrix(nrow=length(gene_counts[[1]]),
                              ncol=length(gene_files_in))

for(i in 1:length(gene_counts)) {
  gene_count_table[,i] <- gene_counts[[i]]
}

rownames(gene_count_table) <- names(gene_counts[[1]])
colnames(gene_count_table) <- pheno_in$Sample #c("SRR1289651", "SRR364315", "SRR364316")
# Sizes                      <- MedianNorm(gene_count_table)
# Salmon_Gene_Counts         <- GetNormalizedMat(gene_count_table, Sizes)
##'-----------------------------------------------------------------------------------------#



##'Salmon - Transcript Level - Get Normalised Counts
##'-----------------------------------------------------------------------------------------#
transcript_files_in    <- list.files("../transcript/",
                                     pattern="*.sf",
                                     full.names=T)
transcript_counts      <- lapply(transcript_files_in,
                                 get_counts)
transcript_count_table <- matrix(nrow=length(transcript_counts[[1]]),
                                 ncol=length(transcript_files_in))

for(i in 1:length(transcript_counts)) {
  transcript_count_table[,i] <- transcript_counts[[i]]
}

rownames(transcript_count_table) <- names(transcript_counts[[1]])
colnames(transcript_count_table) <- pheno_in$Sample #c("SRR1289651", "SRR364315", "SRR364316")
# Sizes                            <- MedianNorm(transcript_count_table)
# Salmon_Trans_Counts              <- GetNormalizedMat(transcript_count_table, Sizes)
##'-----------------------------------------------------------------------------------------#



##'Save Outputs to keep from processing again
##'-----------------------------------------------------------------------------------------#
# save(pheno_in, Salmon_Gene_Counts, Salmon_Trans_Counts, file="Salmon_Data.Rdata")
##'-----------------------------------------------------------------------------------------#


##'Salmon - Transcript Level - Get Normalised Counts
##'-----------------------------------------------------------------------------------------#
Sal_Gene_melt            <- melt(as.matrix(gene_count_table))
colnames(Sal_Gene_melt)  <- c("BioID", "SampleID", "Count")
Sal_Gene_melt$Tool       <- "Salmon"

Sal_Trans_melt           <- melt(as.matrix(transcript_count_table))
colnames(Sal_Trans_melt) <- c("BioID", "SampleID", "Count")
Sal_Trans_melt$Tool      <- "Salmon"

counts_master            <- rbind(Sal_Gene_melt, Sal_Trans_melt)
##'-----------------------------------------------------------------------------------------#



##'Get Annotation for all Ensemble GRCh38 Known Transcripts
##'"MCF2L"
##'-----------------------------------------------------------------------------------------#
transcript_names <- "ALDH1A2"
ensembl          <- useMart("ensembl",
                            dataset="hsapiens_gene_ensembl")
annotation       <- getBM(attributes=c("ensembl_transcript_id",
                                       "ensembl_gene_id",
                                       "external_gene_name",
                                       "description",
                                       "chromosome_name",
                                       "start_position",
                                       "end_position",
                                       "strand"),
                          filters="external_gene_name",
                          values=transcript_names,
                          ensembl)

range_in <- annotation[1,]
range_in <- paste0(range_in$chromosome_name,
                   ":", 
                   (range_in$start_position-1000000),
                   "-",
                   (range_in$end_position+1000000))
anno_range       <- getBM(attributes = c("ensembl_transcript_id",
                                         "ensembl_gene_id",
                                         "external_gene_name",
                                         "description",
                                         "chromosome_name",
                                         "start_position",
                                         "end_position",
                                         "strand"),
                          filters    = "chromosomal_region",
                          values     = range_in, 
                          mart       = ensembl)
##'-----------------------------------------------------------------------------------------#



##'Extract Just the MCF2L Data
##'-----------------------------------------------------------------------------------------#
IDs_to_find      <- c(unique(annotation$ensembl_transcript_id), 
                      unique(annotation$ensembl_gene_id))

df_in            <- counts_master[counts_master$BioID %in% IDs_to_find,]

df_in$Tissue     <- pheno_in[match(df_in$SampleID, pheno_in$Sample),]$Treatment
df_in$Type       <- "Gene"
df_in[grep("ENST", df_in$BioID),]$Type <- "Transcript"
##'-----------------------------------------------------------------------------------------#



##'Visualise Data
##'-----------------------------------------------------------------------------------------#
gg <- ggplot(df_in, aes(x=BioID, y=Count)) +
        geom_bar(stat="identity") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 60,
                                         vjust = 1,
                                         hjust = 1)) +
        ggtitle("ALDH1A2") +
        facet_grid(SampleID ~ .) +
        labs(x = "",
             y = "TPM")

png("ALDH1A2_All.png",
    width=3096, 
    height=4096, 
    units="px", 
    res=300)
print(gg)
dev.off()



for(i in 1:length(unique(df_in$SampleID))) {
  png(paste0(unique(df_in$SampleID)[i], ".png"),
      width=3096, 
      height=4096, 
      units="px", 
      res=300)
  
  gg <- ggplot(df_in[df_in$SampleID == unique(df_in$SampleID)[i],], 
               aes(x=BioID, 
                   y=Count)) +
    geom_bar(stat="identity") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 60,
                                     vjust = 1,
                                     hjust = 1)) +
    ggtitle("ALDH1A2") +
    facet_grid(SampleID ~ .) +
    labs(x = "",
         y = "TPM")
  print(gg)
  dev.off()
}

##'-----------------------------------------------------------------------------------------#



##'Save Data
##'-----------------------------------------------------------------------------------------#
write.csv(df_in, file="ALDH1A2_All.csv")
##'-----------------------------------------------------------------------------------------#



##'Visualise Range
##'-----------------------------------------------------------------------------------------#
IDs_to_find      <- c(unique(anno_range$ensembl_transcript_id), 
                      unique(anno_range$ensembl_gene_id))
df_in            <- counts_master[match(IDs_to_find, counts_master$BioID),]

anno_range       <- anno_range[anno_range$ensembl_transcript_id %in% df_in$BioID |
                               anno_range$ensembl_gene_id %in% df_in$BioID,]

for(i in 1:length(unique(anno_range$external_gene_name))) { 
  anno_in   <- anno_range[anno_range$external_gene_name == unique(anno_range$external_gene_name)[i],]
  IDs_in    <- c(anno_in$ensembl_gene_id[1], 
                 anno_in$ensembl_transcript_id)
  df_sub    <- df_in[df_in$BioID %in% IDs_in,]
  
  gg <- ggplot(df_sub, aes(x=BioID, y=Count)) +
    geom_bar(stat="identity") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 60,
                                     vjust = 1,
                                     hjust = 1)) +
    ggtitle(paste0(unique(anno_in$external_gene_name),
                   " ", 
                   anno_in$chromosome_name[1],
                   ":",
                   anno_in$start_position[1],
                   "-",
                   anno_in$end_position[1])) +
    facet_grid(SampleID ~ .) +
    labs(x = "",
         y = "TPM")
  
  if(max(df_sub$Count) < 1) {
    gg <- gg + scale_y_continuous(limits = c(0, 1))
  }
  
  png(paste0(unique(anno_in$external_gene_name), 
             ".png"),
      width=3096, 
      height=4096, 
      units="px", 
      res=300)
  print(gg)
  dev.off()
}


write.csv(df_in, file="MegaBaseSpan_of_ALDHIA2_TPM.csv")
write.csv(anno_range, file="MegaBaseSpan_of_ALDHIA2_Annotation.csv")

##'-----------------------------------------------------------------------------------------#



##'Visualise Data - Replicates 
##'-----------------------------------------------------------------------------------------#
gg <- ggplot(df_in, aes(x=SampleID, y=Count, colour=BioID, group=BioID)) +
  geom_point(aes(shape=Type, size=4)) +
  geom_line() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60,
                                   vjust = 1,
                                   hjust = 1)) +
  ggtitle(transcript_names) +
  facet_grid(Type ~ Tissue, scales="free_x") +
  labs(x = "",
       y = "FPKM")

png("~/Colin/ALDH1A2_All.png",
    width=3096, 
    height=4096, 
    units="px", 
    res=300)
print(gg)
dev.off()

write.csv(df_in, file="~/Colin/ALDH1A2_All.csv")
##'-----------------------------------------------------------------------------------------#



##'Visualise Data - Replicates Range
##'-----------------------------------------------------------------------------------------#
IDs_to_find      <- c(unique(anno_range$ensembl_transcript_id), 
                      unique(anno_range$ensembl_gene_id))
df_in            <- counts_master[counts_master$BioID %in% IDs_to_find,]
df_in$Tissue     <- pheno_in[match(df_in$SampleID, pheno_in$Sample),]$Treatment
df_in$Type       <- "Gene"
df_in[grep("ENST", df_in$BioID),]$Type <- "Transcript"

anno_range       <- anno_range[anno_range$ensembl_transcript_id %in% df_in$BioID |
                                 anno_range$ensembl_gene_id %in% df_in$BioID,]

for(i in 1:length(unique(anno_range$external_gene_name))) { 
  anno_in   <- anno_range[anno_range$external_gene_name == unique(anno_range$external_gene_name)[i],]
  IDs_in    <- c(anno_in$ensembl_gene_id[1], 
                 anno_in$ensembl_transcript_id)
  df_sub    <- df_in[df_in$BioID %in% IDs_in,]
  
  gg <- ggplot(df_sub, aes(x=SampleID, y=Count, colour=BioID, group=BioID)) +
    geom_point(aes(shape=Type, size=4)) +
    geom_line() +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 60,
                                     vjust = 1,
                                     hjust = 1)) +
    facet_grid(Type ~ Tissue, scales="free_x") +
    labs(x = "",
         y = "FPKM") +
    ggtitle(paste0(unique(anno_in$external_gene_name),
                   " ", 
                   anno_in$chromosome_name[1],
                   ":",
                   anno_in$start_position[1],
                   "-",
                   anno_in$end_position[1])) 
  
  png(paste0("~/Colin/FPKM/Megabase/Plots/",
             unique(anno_in$external_gene_name), 
             ".png"),
      width=3096, 
      height=4096, 
      units="px", 
      res=300)
  print(gg)
  dev.off()
}


write.csv(df_in, file="~/Colin/FPKM/Megabase/MegaBaseSpan_of_ALDHIA2_TPM_NOF_OA.csv")
write.csv(anno_range, file="~/Colin/FPKM/Megabase/MegaBaseSpan_of_ALDHIA2_Annotation_NOF_OA.csv")
##'-----------------------------------------------------------------------------------------#









