#-----------------------------------------------------------------------------------------------------------------------------#
#  Author      : Andrew & Ruddy                                                                                               |
#  Language    : R Statistical Programming Language                                                                           |
#  Optimised   :                                                                                                              |
#  Data Owner  : Ruddy                                                                                                        |
#  Description : hoRse Package Design                                                                                         |
#-----------------------------------------------------------------------------------------------------------------------------#

hoRse_setup <- function() {
  ##LOAD_PACKAGES_&_PREPARE---------------------------------------------------------------------------------------------------#
  message("Loading hoRse Dependences")
  pkgTest <- function(x)
  {
    if (!require(x,character.only = TRUE))
    {
      source("http://bioconductor.org/biocLite.R")
      biocLite(x)
      require(x)
    }
  }
  pkgTest("rtracklayer")
  pkgTest("FDb.InfiniumMethylation.hg19")
  pkgTest("GenomicFeatures")
  pkgTest("GenomicRanges")
  pkgTest("lumi")
  pkgTest("ggplot2")
  
  message("Creating hoRse Environment")
  hoRse.env <- new.env()
  ##'-------------------------------------------------------------------------------------------------------------------------#
  
  ##PRE PROCESSING -----------------------------------------------------------------------------------------------------------#
  message("Fetching Human Genome (hg19)")
  hoRse.env$InfiniumMethylation         <- features(FDb.InfiniumMethylation.hg19)
  message("Fetching hg19 Methylation Metadata")
  hoRse.env$met                         <- metadata(FDb.InfiniumMethylation.hg19)
  message("Finalising Reference Dataset")
  genome(hoRse.env$InfiniumMethylation) <- hoRse.env$met[which(hoRse.env$met[,"name"]=="Genome"),"value"]
  hoRse.env$InfiniumMethylation         <- sort(hoRse.env$InfiniumMethylation)
  
  ##'Tests
#   show(hoRse.env$InfiniumMethylation)
#   length(hoRse.env$InfiniumMethylation$probeType)
  
  ##TXF REFERENCE ------------------------------------------------------------------------------------------------------------#
  message("Dowloading Additional Resources (This might take a few minutes, grab a coffee...)")
  hoRse.env$fdb                         <-makeFeatureDbFromUCSC(genome="hg19",
                                                                track="Txn Factor ChIP",
                                                                tablename="wgEncodeRegTfbsClusteredV3")
  hoRse.env$TXF_Ruddy_2014 <- features(hoRse.env$fdb)
  hoRse.env$TXF_Ruddy_2014 <- sort(hoRse.env$TXF_Ruddy_2014)

  ##'Tests
#   show(hoRse.env$TXF_Ruddy_2014)

  ##TXF SUBSET ---------------------------------------------------------------------------------------------------------------#
  message("Subsetting TXF Data......")
  hoRse.env$Array <- subsetByOverlaps(hoRse.env$TXF_Ruddy_2014,hoRse.env$InfiniumMethylation)
  table(names(hoRse.env$Array))
  hoRse.env$All_161_TXFs_old<-as.data.frame(table(names(hoRse.env$Array)))
  hoRse.env$All_161_TXFs_old
  
  hoRse.env$cov                <- countOverlaps(hoRse.env$TXF_Ruddy_2014, hoRse.env$InfiniumMethylation)
  hoRse.env$cov_matrix         <- as.matrix(hoRse.env$cov)
  hoRse.env$split_cov_matrix   <- split(hoRse.env$cov_matrix, rownames(hoRse.env$cov_matrix))
  hoRse.env$Number_binding_TXF <- as.matrix(lapply(hoRse.env$split_cov_matrix,sum))
  
  hoRse.env$All_161_TXFs       <- hoRse.env$Number_binding_TXF
  hoRse.env$All_161_TXFs       <- as.data.frame(hoRse.env$All_161_TXFs)
  hoRse.env$All_161_TXFsv2     <- hoRse.env$All_161_TXFs_old
  hoRse.env$All_161_TXFsv2[,2] <- hoRse.env$All_161_TXFs[1]
  hoRse.env$All_161_TXFs       <- hoRse.env$All_161_TXFsv2
  message("hoRse.env Preliminary Setup Complete!")
  ##'-------------------------------------------------------------------------------------------------------------------------#
  return(hoRse.env)
}
##'TEST FUNCTION--------------#
# foo <- hoRse_setup()
##'---------------------------#
##'---------------------------------------------------------------------------------------------------------------------------#



##COUNT DATA EXTRACT ---------------------------------------------------------------------------------------------------------#
TFnames4cpg_v2 <- function (hoRse.environment, cpg_vec) {
  cov_matrix       <- as.matrix( countOverlaps(hoRse.environment$TXF_Ruddy_2014, 
                                               hoRse.environment$InfiniumMethylation[cpg_vec]))
  split_cov_matrix <- split(cov_matrix, rownames(cov_matrix))
  txf_counts       <- as.matrix(lapply(split_cov_matrix, sum))
  return(data.frame(txf=names(unlist(txf_counts[,1])), count=as.vector(unlist(txf_counts[,1]))))
}
##'TEST FUNCTION--------------#
# df <- TFnames4cpg_v2(foo, Adipose_tissue_hypercpgs)
##'---------------------------#
##'---------------------------------------------------------------------------------------------------------------------------#


#___________________Inverse function to know for a given name and cpglist what cpgs bind to this TXF chonsen

#the arguments of the funtion are:X= the list of cpgs. Y= the name of the TXF between "" ( for example("ZEB1"))
cpgs4TXF <- function(X,Y, env_in) {
  names(subsetByOverlaps(env_in$InfiniumMethylation[ X],env_in$TXF_Ruddy_2014[names(env_in$TXF_Ruddy_2014)==Y]))
}
# cpgs4TXF(OA_cpglist_hyper,"ZEB1")
cpgs4TXF(hyper_in[[2]], "ZEB1", hoRse.env)


##PREPARE COUNT DATA ---------------------------------------------------------------------------------------------------------#
hoRse_counts <- function(names_vec, hypo_list, hyper_list, hoRse.environment) {
  df <- c()
  for(i in 1:length(names_vec)) {
    df_tmp    <- c()
    hypo_tmp  <- TFnames4cpg_v2(hoRse.environment, hypo_list[[i]])
    hypo_tmp$class     <- "hypo"
    hypo_tmp$data      <- names_vec[i]
#     hypo_tmp$Perc      <- (as.numeric(hypo_tmp$count) * 100)/sum(as.numeric(hypo_tmp$count))
    hypo_tmp$Perc      <- (as.numeric(hypo_tmp$count) * 100)/length(hypo_list[[i]])
    hypo_tmp$Not_Bound <- length(hypo_list[[i]])-hypo_tmp$Perc 
    
    hyper_tmp <- TFnames4cpg_v2(hoRse.environment, hyper_list[[i]])
    hyper_tmp$class     <- "hyper"
    hyper_tmp$data      <- names_vec[i]
#     hyper_tmp$Perc      <- (as.numeric(hyper_tmp$count) * 100)/sum(as.numeric(hyper_tmp$count))
    hyper_tmp$Perc      <- (as.numeric(hyper_tmp$count) * 100)/length(hyper_list[[i]])
    hyper_tmp$Not_Bound <- length(hyper_list[[i]])-hyper_tmp$Perc 
      
    RRT <- (hypo_tmp$Perc / hyper_tmp$Perc)    
    hyper_tmp$RRT <- RRT
    hypo_tmp$RRT  <- RRT
    df <- rbind(df, hyper_tmp, hypo_tmp)
  } 
  
  df                 <- merge(df, hoRse.environment$All_161_TXFs, by.x="txf", by.y="Var1")
  df$Array_Perc      <- (as.numeric(df$Freq) * 100)/length(hoRse.environment$InfiniumMethylation$probeType) 
  df$Freq <- as.numeric(unlist(df$Freq))
  df$Array_Not_Bound <- length(hoRse.environment$InfiniumMethylation$probeType)-as.numeric(df$Freq)
  df$RRTlog          <- log2(df$RRT)
  return(df)
}
##'TEST FUNCTION--------------#
names_in <- c("Adipose", "hPSC", "OB")
hypo_in  <- list(Adipose_tissue_hypocpgs, hPSC_tissue_hypocpgs, OB_array2_cpglist_hypo)
hyper_in <- list(Adipose_tissue_hypercpgs, hPSC_tissue_hypercpgs, OB_array2_cpglist_hyper)
counts   <- hoRse_counts(names_in, hypo_in, hyper_in, hoRse.env)
# 
names_in <- c("hPSC")
hypo_in  <- list(hPSC_tissue_hypocpgs)
hyper_in <- list(hPSC_tissue_hypercpgs)
counts   <- hoRse_counts(names_in, hypo_in, hyper_in, hoRse.env)
# 
# names_in <- c("OB")
# hypo_in  <- list(OB_array2_cpglist_hypo)
# hyper_in <- list(OB_array2_cpglist_hyper)
# counts   <- hoRse_counts(names_in, hypo_in, hyper_in, foo)
##'---------------------------#
##'---------------------------------------------------------------------------------------------------------------------------#

##PLOT COUNT DATA ------------------------------------------------------------------------------------------------------------#
hoRse_plot_counts <- function(df, txf_of_interest, metric="count", display="both",
                              scale=F) {
  require(ggplot2)
  if(!missing(txf_of_interest)) {
    df <- df[grep(paste(txf_of_interest, collapse="|"), df$txf),] 
  }
  
  if(display == "hyper") {
    df <- df[df$class == "hyper",] 
  } else if(display == "hypo") {
    df <- df[df$class == "hypo",]
  }
  
  if(metric == 'count') {
    g <- ggplot(df, aes(x=factor(txf), y=count, fill=class) )  
  } else if(metric == 'perc') {
    g <- ggplot(df, aes(x=factor(txf), y=Perc, fill=class) )
    if(scale == T) { g <- g + scale_y_continuous(limits=c(0, 100)) }
  } else if(metric == 'binding') {
    g <- ggplot(df, aes(x=factor(txf), y=Not_Bound, fill=class) )  
  } else if(metric == 'freq') {
    g <- ggplot(df, aes(x=factor(txf), y=Freq, fill=class) )  
  } else if(metric == 'array_perc') {
    g <- ggplot(df, aes(x=factor(txf), y=Array_Perc, fill=class) )  
  } else if(metric == 'rrt') {
    g <- ggplot(df, aes(x=factor(txf), y=RRT, fill=class) )  
  } else if(metric == 'rrtlog') {
    g <- ggplot(df, aes(x=factor(txf), y=RRTlog, fill=class) )  
  } else { break } 
  
  g <- g + theme_bw() +
#           coord_flip() + #guides(fill=FALSE) + 
          theme(axis.text.y=element_text(size=7)) +             
          facet_grid( ~ data, space="free", scales="free_y") +
          labs(y=metric, x="Transcription Factor")
          
  if(display != "both") {
    g <- g + geom_bar(stat = "identity", position=position_dodge())
  } else {
    g <- g + geom_bar(stat = "identity")
  }
  return(g)
}
##'TEST FUNCTION--------------#
# hoRse_plot_counts(counts)
# hoRse_plot_counts(counts, metric="perc", scale=F, display='hyper', txf_of_interest=c("POLR2A"))
# hoRse_plot_counts(counts, metric="perc", txf_of_interest=c("ZEB1", "EZH2", "CTBP2", "SUZ12"))
# hoRse_plot_counts(counts, metric="binding")
# hoRse_plot_counts(counts, metric="freq")
# hoRse_plot_counts(counts, metric="array_perc")
# hoRse_plot_counts(counts, metric="rrt")
# hoRse_plot_counts(counts, c("ZNF274", "RUNX3", "IRF4", "HMGN3", "EBF1"), metric="rrt")
# hoRse_plot_counts(counts, c("ZNF274", "RUNX3", "IRF4", "HMGN3", "EBF1"), metric="count", display="hyper")
# hoRse_plot_counts(counts, c("ZNF274", "RUNX3", "IRF4", "HMGN3", "EBF1"), "perc", display="hyper")
# hoRse_plot_counts(counts, c("ZNF274", "RUNX3", "IRF4", "HMGN3", "EBF1"), "binding")
# hoRse_plot_counts(counts, c("ZNF274", "RUNX3", "IRF4", "HMGN3", "EBF1"), "freq")

hoRse_plot_counts(counts, metric="perc", display="both", txf_of_interest=c("EZH2","SUZ12"))
hoRse_plot_counts(counts, c("ZNF274", "RUNX3", "IRF4", "HMGN3", "EBF1"), metric="rrt", display="hyper")
##'---------------------------#
##'---------------------------------------------------------------------------------------------------------------------------#

# ggplot(counts, aes(x=class, y=Perc, fill=txf)) + geom_bar(stat = "identity") + facet_grid(. ~ data)

##LINE PLOT COUNT DATA -------------------------------------------------------------------------------------------------------#
hoRse_plot_counts_line <- function(df, txf_of_interest, metric="count", display="both",
                              scale=F, sort=T, filter=T) {
  require(ggplot2)
  if(!missing(txf_of_interest)) {
    df <- df[grep(paste(txf_of_interest, collapse="|"), df$txf),] 
  }
  
  if(filter == T) {
    df <- df[!is.na(df$RRTlog) & !is.infinite(df$RRTlog),]
  }
  
  if(sort==T) {
    df$txf <- as.character(df$txf)
    if(metric == "rrt") {
      df     <- df[!is.na(df$RRT) & !is.infinite(df$RRT),]
      df     <- df[order(df$data, df$RRT),]
      df$txf <- factor(df$txf, levels=df$txf)
    } else if (metric == 'rrtlog') {
      df     <- df[!is.na(df$RRTlog) & !is.infinite(df$RRTlog),]
      df     <- df[order(df$data, df$RRTlog),]
      df$txf <- factor(df$txf, levels=df$txf)
    }
    rownames(df) <- c(1:nrow(df))
  }
  
  if(display == "hyper") {
    df <- df[df$class == "hyper",] 
  } else if(display == "hypo") {
    df <- df[df$class == "hypo",]
  }
  
  if(metric == 'count') {
    g <- ggplot(df, aes(x=factor(txf), y=count, group=class, colour=class) )  
  } else if(metric == 'perc') {
    g <- ggplot(df, aes(x=factor(txf), y=Perc, group=class, colour=class) )
    if(scale == T) { g <- g + scale_y_continuous(limits=c(0, 100)) }
  } else if(metric == 'binding') {
    g <- ggplot(df, aes(x=factor(txf), y=Not_Bound, group=class, colour=class) )  
  } else if(metric == 'freq') {
    g <- ggplot(df, aes(x=factor(txf), y=Freq, group=class, colour=class) )  
  } else if(metric == 'array_perc') {
    g <- ggplot(df, aes(x=factor(txf), y=Array_Perc, group=class, colour=class) )  
  } else if(metric == 'rrt') {
    g <- ggplot(df, aes(x=factor(txf), y=RRT, group=data, colour=data) )  
  } else if(metric == 'rrtlog') {
    g <- ggplot(df, aes(x=(factor(txf)), y=RRTlog, group=data, colour=data) )  +
      facet_grid(.~data, space="free", scales="free_y") +
      geom_hline(yintercept=0)
  } else { break } 
  
  g <- g + geom_line() + geom_point() +
           theme_bw() 
#            theme(axis.text.y=element_text(size=7)) +             
    
#     facet_grid(. ~ data, space="free", scales="free_y") +
#            labs(y=metric, x="Transcription Factor")
  
  if(length(unique(df$txf)) < 10) {
    g <- g + theme(axis.text.x = element_text(angle = 90, hjust = 1, size=10))
  } else {
    g <- g + theme(axis.text.x = element_text(angle = 90, hjust = 1, size=4))
  }
    
  return(g)
}
##'TEST FUNCTION--------------#
# hoRse_plot_counts_line(counts, c("ZNF274", "RUNX3", "IRF4", "HMGN3", "EBF1"), metric="perc")
# hoRse_plot_counts_line(counts, metric="perc")
# 
# 
# hoRse_plot_counts(counts, c("ZEB1", "EZH2"), metric="perc")
# hoRse_plot_counts_line(counts, c("ZEB1", "EZH2"), metric="perc")
# 
# hoRse_plot_counts_line(df=counts, metric="perc", sort=F, scale=T, filter=T, txf_of_interest=c("ZEB1", "EZH2"))
# 
# hoRse_plot_counts_line(df=counts, metric="rrtlog", sort=T, scale=T)
# 
# hoRse_plot_counts_line(df=counts, txf_of_interest=c("ZEB1", "EZH2", "CTBP2", "SUZ12"), metric="rrtlog", sort=T)
##'---------------------------#
##'---------------------------------------------------------------------------------------------------------------------------#

##SUBET COUNT DATA FOR POTENTIALLY INTERESTING TXF ---------------------------------------------------------------------------#
hoRse_rrtlog_extract <- function(counts, data, N=5) {
  if(!missing(data)) {
    counts_tmp <- counts[counts$data == data,] #& counts$class == "hyper",]
  }
  
  counts_tmp <- counts
  
  max_temp <- sort(counts_tmp$RRTlog, decreasing=T)
  max_txf <- counts_tmp[grep(paste(max_temp[max_temp != Inf][1:(N*2)], collapse="|"), counts_tmp$RRTlog),]
  max_txf <- max_txf[1:(N*2),]
  
  min_temp <- sort(unique(counts_tmp$RRTlog), decreasing=F)
  min_txf <- counts_tmp[grep(paste(min_temp[min_temp != -Inf][1:(N*2)], collapse="|"), counts_tmp$RRTlog),]
  min_txf <- min_txf[1:(N*2),]
  
  return(counts_tmp[counts_tmp == as.vector(max_txf$txf) || counts_tmp == as.vector(min_txf$txf),])
}
##'---------------------------------------------------------------------------------------------------------------------------#