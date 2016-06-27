#-------------------------------------------------------------------------------------------#
#  Author      : Andrew J Skelton                                                           |
#  Language    : R Statistical Programming Language                                         |
#  Study       : Methylation of Cartilage Study - Hip / Knee OA                             |
#  Data Owner  : Newcastle University - Prof. John Loughlin, Dr. Louise Reynard             |
#  Description : Extract CpGs from Gene Short Names                                         |
#-------------------------------------------------------------------------------------------#



##'Find CpGs of Gene Short Names
##'-----------------------------------------------------------------------------------------#
genes_of_interest <- c("ABCD4", "BRMS1")

library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
annotation_env    <- IlluminaHumanMethylation450kanno.ilmn12.hg19@data
annotation_other  <- as.data.frame(annotation_env$Other)

cpg_sites         <- annotation_other[grep(paste(c(paste0(genes_of_interest, ";"),
                                                   paste0(genes_of_interest, "$")),
                                                 collapse="|"),
                                           annotation_other$UCSC_RefGene_Name),]

cpg_probes        <- rownames(cpg_sites)
##'-----------------------------------------------------------------------------------------#