setwd("~/Young/Pipeline_Analysis/cuffdiff/grch38_No_Novel_v2/")

gene_map                     <- read.table("tss_groups.fpkm_tracking", sep="\t", header=T, stringsAsFactors=F)
gene_exp_diff                <- read.table("gene_exp.diff", sep="\t", header=T, stringsAsFactors=F)
tail(gene_map)
# gene_exp_out <- gene_exp_diff[(gene_exp_diff$value_1 > 5 | gene_exp_diff$value_2 > 5) & 
#                                abs(gene_exp_diff$log2.fold_change) > 1.5 & 
#                                gene_exp_diff$q_value < 0.05,]
gene_exp                     <- gene_exp_diff
gene_exp$gene                <- gsub(",", ";", gene_exp$gene)
gene_exp$abs.fc              <- abs(gene_exp$log2.fold_change.)
gene_exp$min.val.5           <- "No"
gene_exp[gene_exp$value_1 > 5 | gene_exp$value_2 > 5,]$min.val.5 <- "Yes"
gene_exp_out                 <- gene_exp[,c(1:5,8,6,9,16,10,15,12:14)]
gene_exp_out                 <- cbind(gene_exp[,1:5], gene_exp[,], gene_exp[,15], gene_exp[,11:14])

write.table(gene_exp_out, file="cuffdiff_genes_out.txt", sep="\t", quote=F, col.names=T, row.names=F)
# write.csv(gene_exp_out, file="cuffdiff_genes_FPKM.csv", quote=F)

isoform_map                    <- read.table("isoforms.fpkm_tracking", sep="\t", header=T, stringsAsFactors=F)
isoform_exp_diff               <- read.table("isoform_exp.diff", sep="\t", header=T, stringsAsFactors=F)
# isoform_exp_out                <- isoform_exp_diff[(isoform_exp_diff$value_1 > 5 | isoform_exp_diff$value_2 > 5) & 
#                                                     abs(isoform_exp_diff$log2.fold_change) > 1.5 & 
#                                                     isoform_exp_diff$q_value < 0.05,]
isoform_exp_out                <- isoform_exp_diff
colnames(isoform_exp_out)[1]   <- colnames(isoform_map)[1]
isoform_exp_out_anno           <- merge(isoform_exp_out, isoform_map[,c(1:5, 8)], by=colnames(isoform_exp_out)[1])
isoform_exp_out_anno$abs.fc    <- abs(isoform_exp_out_anno$log2.fold_change.)
isoform_exp_out_anno$min.val.5 <- "No"
isoform_exp_out_anno[isoform_exp_out_anno$value_1 > 5 | isoform_exp_out_anno$value_2 > 5,]$min.val.5 <- "Yes"
isoform_exp_out_anno_out       <- isoform_exp_out_anno[,c(1,2,16,3,4,19,15,5,8,6,9,21,10,20,12,13,14)]

isoform_exp_out_anno_out[isoform_exp_out_anno_out$q_value < 0.05,]
write.table(isoform_exp_out_anno_out, file="cuffdiff_isoforms_out.txt", sep="\t", quote=F,
            col.names=T, row.names=F)
nrow(isoform_map[isoform_map$class_code == 'j' &
                 (isoform_map$NOF_FPKM > 5 | isoform_map$OA_FPKM > 5),])

# write.csv(isoform_exp_out_anno_out, file="cuffdiff_isoforms_FPKM.csv", quote=F)

novel_isoforms   <- isoform_exp_out_anno_out[isoform_exp_out_anno_out$gene == "-",]
novel_isoforms_j <- isoform_exp_out_anno_out[isoform_exp_out_anno_out$class_code == "j",]
# write.csv(novel_isoforms_j, file="cuffdiff_isoforms_FPKM-5_log2FC-1.5_qval-0.05_class-j.csv")


foo  <- isoform_exp_out_anno_out[grep("MCF2L$", isoform_exp_out_anno_out$gene),]
fooB <- foo[foo$class_code == "=",]
# write.csv(fooB, file="../../../NOF_OA_MCF2L_Colin_Tuxedo.csv")

foo  <- gene_exp_diff[grep("MCF2L$", gene_exp_diff$gene),]
# write.csv(foo, file="../../../NOF_OA_GeneLevel_MCF2L_Colin_Tuxedo.csv")


transcript_names <- isoform_exp_out_anno_out$nearest_ref_id
ensembl          <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
annotation       <- getBM(attributes=c("ensembl_transcript_id", "ensembl_gene_id", "external_gene_name", "description", "transcript_biotype"),  #"ensembl_gene_id
                          filters="ensembl_transcript_id", 
                          values=transcript_names, 
                          ensembl)
listAttributes(ensembl)

foo <- merge(isoform_exp_out_anno_out, annotation, by.x="nearest_ref_id", by.y="ensembl_transcript_id")
head(foo)
foo_out <- foo[,c(1:4, 20:21, 5:17)]
foo_out$description <- gsub(",", ";", foo_out$description)
head(foo_out)


write.csv(foo_out, file="cuffdiff_isoform_FPKM_No_Novel.csv", quote=F)
