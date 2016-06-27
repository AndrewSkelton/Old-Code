genenames = c("ACAN", "SOX9", "COL2A1")
illuminanames <- c("ILMN_2326509", "ILMN_2147517", "ILMN_1712678", "ILMN_2224103")
chrIn <- c("15:88803443:88875354:1")

mart_ob <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")

IdMultimapper <-function (x, type="gene_name", mart_object) {
  if(type == "gene_name") {
    return(getBM(attributes=c("refseq_mrna","hgnc_symbol","ensembl_gene_id", "illumina_humanht_12_v4", "ensembl_transcript_id", "description"), 
                 filters="hgnc_symbol", values=x, mart=mart_object))
  } else if(type == "illumina") {
    return(getBM(attributes=c("refseq_mrna","hgnc_symbol","ensembl_gene_id", "illumina_humanht_12_v4", "ensembl_transcript_id", "description"), 
                 filters="illumina_humanht_12_v4", values=x, mart=mart_object))
  } else if(type == "ensembl_transcript") {
    return(getBM(attributes=c("refseq_mrna","hgnc_symbol","ensembl_gene_id", "illumina_humanht_12_v4", "ensembl_transcript_id", "description"), 
                 filters="ensembl_gene_id", values=x, mart=mart_object))
  } else if(type == "chromosomal") {
    return(getBM(attributes=c("refseq_mrna","hgnc_symbol","ensembl_gene_id", "illumina_humanht_12_v4", "ensembl_transcript_id", "description"), 
                 filters="chromosomal_region", values=x, mart=mart_object))
  } else {
    message("Error: bad type argument")
    return()
  }
}

ID2TXF <- function(x, type, mart, hoRse.environment) {
  mapping  <- IdMultimapper(x, type, mart)
  vec_in   <- mapping$refseq_mrna
  transcripts_GR <- transcripts(hoRse.environment$refGene_txdb)[
    mcols(transcripts(hoRse.environment$refGene_txdb))$tx_name %in% vec_in]
  splitTXF     <- split(hoRse.environment$TXF_Ruddy_2014, 
                        names(hoRse.environment$TXF_Ruddy_2014))
  
  overlapTable <- as.data.frame(findOverlaps(splitTXF, transcripts_GR ))
  overlapTable$txf         <- names(splitTXF[overlapTable$queryHits])
  overlapTable$inputID     <- mcols(GR_genename_Andrew[overlapTable$subjectHits])$tx_name
  overlapTable$Gene        <- mapvalues(overlapTable$inputID, from=input_data$refseq_mrna, to=input_data$hgnc_symbol)
  overlapTable$Illumina    <- mapvalues(overlapTable$inputID, from=input_data$refseq_mrna, to=input_data$illumina_humanht_12_v4)
  overlapTable$Ensembl     <- mapvalues(overlapTable$inputID, from=input_data$refseq_mrna, to=input_data$ensembl_transcript_id)
  overlapTable$Description <- mapvalues(overlapTable$inputID, from=input_data$refseq_mrna, to=input_data$description)
  return(overlapTable)
}


# refGene_txdb       <- makeTranscriptDbFromUCSC(genome="hg19", tablename="refGene")

# input_data         <- IdMultimapper(illuminanames, "illumina")
input_data         <- IdMultimapper(chrIn, "chromosomal", mart_ob)
vec_in             <- input_data$refseq_mrna

GR_genename_Andrew <- transcripts(refGene_txdb)[mcols(transcripts(refGene_txdb))$tx_name %in% vec_in]

splitTXF     <- split(hoRse.env$TXF_Ruddy_2014, names(hoRse.env$TXF_Ruddy_2014))
overlapTable <- as.data.frame(findOverlaps(splitTXF, GR_genename_Andrew ))

overlapTable$txf         <- names(splitTXF[overlapTable$queryHits])
overlapTable$inputID     <- mcols(GR_genename_Andrew[overlapTable$subjectHits])$tx_name
overlapTable$Gene        <- mapvalues(overlapTable$inputID, from=input_data$refseq_mrna, to=input_data$hgnc_symbol)
overlapTable$Illumina    <- mapvalues(overlapTable$inputID, from=input_data$refseq_mrna, to=input_data$illumina_humanht_12_v4)
overlapTable$Ensembl     <- mapvalues(overlapTable$inputID, from=input_data$refseq_mrna, to=input_data$ensembl_transcript_id)
overlapTable$Description <- mapvalues(overlapTable$inputID, from=input_data$refseq_mrna, to=input_data$description)

output <- overlapTable[,c(3:(ncol(overlapTable)))]
output
output[output$Illumina == input_data$illumina_humanht_12_v4[1],]





