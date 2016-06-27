setwd("~/Young/Pipeline_Analysis/salmon/grch38/alignment")

source("http://bioconductor.org/biocLite.R")
biocLite("EBSeq")

library("EBSeq")
library("biomaRt")
library("ggplot2")

get_counts = function(file) {
  dat               <- read.table(file, skip=10, sep="\t", stringsAsFactors=F, row.names=1)
  raw_counts        <- dat[,4]
  names(raw_counts) <- rownames(dat)
  return(raw_counts)
}

get_results = function(norm, cond, res, ann) {
  if (length(cond) != ncol(norm)) {
    stop("Number of conditions does not equal number of samples")
  }
  cond1             <- which(cond==levels(cond)[1])
  cond2             <- which(cond==levels(cond)[2])
  results           <- data.frame(row.names=rownames(res), transcript_id=row.names(res), 
                                  gene_id=NA, baseMean=NA, cond1Mean=NA, cond2Mean=NA, log2FC=NA, 
                                  PPEE=NA, PPDE=NA)
  results$gene_id   <- ann[match(row.names(res), ann$ensembl_transcript_id), "ensembl_gene_id"]
  results$baseMean  <- apply(norm[row.names(res),], 1, mean)
  results$cond1Mean <- apply(norm[row.names(res),cond1], 1, mean)
  results$cond2Mean <- apply(norm[row.names(res),cond2], 1, mean)
  results$log2FC    <- apply(norm[row.names(res),], 1, function(x) log2(mean(x[cond2]))-log2(mean(x[cond1])))
  results$PPEE      <- res[row.names(res),"PPEE"]
  results$PPDE      <- res[row.names(res),"PPDE"]
  return(results)
}

files_in    <- list.files("./", pattern="*.sf")
counts      <- lapply(files_in, get_counts)
count_table <- matrix(nrow=length(counts[[1]]), ncol=length(files_in))

i = 1
for (c in counts) {
  count_table[,i] = c
  i = i + 1
}
rownames(count_table) <- names(counts[[1]])
colnames(count_table) <- gsub(".sf", "", files_in)

pheno <- data.frame(
  row.names = c("F1_S1","F1_S2","F1_S3","F1_S4","F1_S5","F1_S6","F1_S7","F1_S8",
                "F2_S1","F2_S2","F2_S3","F2_S4","F2_S5","F2_S6","F2_S7","F2_S8"),
  condition = c("OA","NOF","OA","OA","NOF","OA","NOF","OA","OA","NOF","OA","OA",
                "NOF","NOF","OA","OA"),
  libType   = c(rep("paired-end", 16)), 
  countName = list.files("./", pattern="*.sf", full.names=T)
)

transcript_names <- rownames(count_table)
ensembl          <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
annotation       <- getBM(attributes=c("ensembl_transcript_id", "ensembl_gene_id", "external_gene_name", "description"),  #"ensembl_gene_id
                          filters="ensembl_transcript_id", 
                          values=transcript_names, 
                          ensembl)

annotation            <- annotation[match(transcript_names, annotation$ensembl_transcript_id),]
transcript_gene_names <- annotation$ensembl_gene_id

transcript_sizes      <- MedianNorm(count_table)
ng_list               <- GetNg(transcript_names, transcript_gene_names)
iso_ng_trun           <- ng_list$IsoformNgTrun

iso_eb_out            <- EBTest(Data=data.matrix(count_table), NgVector=iso_ng_trun,
                                Conditions=pheno$condition, sizeFactors=transcript_sizes, maxround=5)
iso_pp                <- GetPPMat(iso_eb_out)
iso_de                <- rownames(iso_pp)[which(iso_pp[,"PPDE"]>=0.95)]

norm_data             <- GetNormalizedMat(data.matrix(count_table), transcript_sizes)

#results for filtered set (only those tested)
results_table         <- get_results(norm_data, pheno$condition, iso_pp, annotation)
results_table$abs.fc  <- abs(results_table$log2FC)

results_merged        <- merge(results_table, annotation, by.x="transcript_id", by.y="ensembl_transcript_id")

results_table_out     <- results_merged[,c(1,10:12,3:6,9,7:8)]
write.csv(results_table_out, file="NOF_OA_Salmon_EBSeq_Filtered.csv")

sorted_res            <- results_table[with(results_table, order(-abs(log2FC))),]
head(sorted_res, 100)
sorted_res['ENST00000394299',]


