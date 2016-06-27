#-----------------------------------------------------------------------------------------------------------------------------#
#  Author      : Andrew & Ruddy                                                                                               |
#  Language    : R Statistical Programming Language                                                                           |
#  Optimised   :                                                                                                              |
#  Data Owner  : Ruddy                                                                                                        |
#  Description : hoRse Package Test                                                                                           |
#-----------------------------------------------------------------------------------------------------------------------------#

##STEP 1 SETUP DATABASES -----------------------------------------------------------------------------------------------------#
hoRse.env <- hoRse_setup()
##Or 
#replace the paths to the Rdata object
load("~/Documents/Repos/Methylation_Analysis_Project/Build_1/pre_load.Rdata")
pkgTest("rtracklayer")
pkgTest("FDb.InfiniumMethylation.hg19")
pkgTest("GenomicFeatures")
pkgTest("GenomicRanges")
pkgTest("lumi")
pkgTest("ggplot2")
pkgTest("annotate")
pkgTest("biomaRt")
pkgTest("plyr")
##'---------------------------------------------------------------------------------------------------------------------------#

##STEP 2 GET SOME TRANSCRIPTION FACTORS --------------------------------------------------------------------------------------#
#replace the paths to the Rdata objects
# load("~/Documents/Repos/Methylation_Analysis_Project/Build_1/Adipose_tissue_hypercpgs.Rdata")
# load("~/Documents/Repos/Methylation_Analysis_Project/Build_1/Adipose_tissue_hypocpgs.Rdata")
# load("~/Documents/Repos/Methylation_Analysis_Project/Build_1/OB_array2_cpglist_hyper.Rdata")
# load("~/Documents/Repos/Methylation_Analysis_Project/Build_1/OB_array2_cpglist_hypo.Rdata")
# load("~/Documents/Repos/Methylation_Analysis_Project/Build_1/hPSC_tissue_hypercpgs.Rdata")
# load("~/Documents/Repos/Methylation_Analysis_Project/Build_1/hPSC_tissue_hypocpgs.Rdata")
#Or
load("~/Documents/Repos/Methylation_Analysis_Project/Build_1/Test_Samples_Bundle.Rdata")
##'---------------------------------------------------------------------------------------------------------------------------#

##STEP 3 SET UP TRANSCRIPTION FACTORS FOR PROCESSING -------------------------------------------------------------------------#
#3 TISSUE TYPES
names_in <- c("Adipose", "hPSC", "OB", "Thymus")
hypo_in  <- list(Adipose_tissue_hypocpgs, hPSC_tissue_hypocpgs, OB_array2_cpglist_hypo, Thymus_tissue_hypocpgs)
hyper_in <- list(Adipose_tissue_hypercpgs, hPSC_tissue_hypercpgs, OB_array2_cpglist_hyper, Thymus_tissue_hypercpgs)
counts   <- hoRse_counts(names_in, hypo_in, hyper_in, hoRse.env)
#OR
#1 TISSUE TYPE
names_in <- c("hPSC")
hypo_in  <- list(hPSC_tissue_hypocpgs)
hyper_in <- list(hPSC_tissue_hypercpgs)
counts   <- hoRse_counts(names_in, hypo_in, hyper_in, hoRse.env)
#OR
#1 TISSUE TYPE
names_in <- c("OB")
hypo_in  <- list(OB_array2_cpglist_hypo)
hyper_in <- list(OB_array2_cpglist_hyper)
counts   <- hoRse_counts(names_in, hypo_in, hyper_in, hoRse.env)
#OR
#1 TISSUE TYPE
names_in <- c("Thymus")
hypo_in  <- list(Thymus_tissue_hypocpgs)
hyper_in <- list(Thymus_tissue_hypercpgs)
counts   <- hoRse_counts(names_in, hypo_in, hyper_in, hoRse.env)
##'---------------------------------------------------------------------------------------------------------------------------#

#FOR A SINGLE TISSUE TYPE
#Look at the RRTlog metric
hoRse_plot_counts_line(counts, metric="rrtlog")
#Extract the top 4 highest and lowest on the log scale that 
#are potentially interesting
counts_2 <- hoRse_rrtlog_extract(counts, N=4)
#Plot the Percentage Metric
# hoRse_plot_counts_line(counts_2, metric="rrtlog")
hoRse_plot_counts(counts_2, metric="perc")

##FOR MULTIPLE TISSUE TYPES
#Plot intividual rrt log values
hoRse_plot_counts_line(counts, metric="rrtlog", data="Adipose")
hoRse_plot_counts_line(counts, metric="rrtlog", data="hPSC")
hoRse_plot_counts_line(counts, metric="rrtlog", data="OB")
hoRse_plot_counts_line(counts, metric="rrtlog", data="Thymus")

#extract interesting TXFs and plot by tissue type
counts_2 <- hoRse_rrtlog_extract(counts, "Adipose", 4)
hoRse_plot_counts(counts_2, metric="perc")

counts_2 <- hoRse_rrtlog_extract(counts, "hPSC", 4)
hoRse_plot_counts(counts_2, metric="perc")

counts_2 <- hoRse_rrtlog_extract(counts, "OB", 4)
hoRse_plot_counts(counts_2, metric="perc")

counts_2 <- hoRse_rrtlog_extract(counts, "Thymus", 4)
hoRse_plot_counts(counts_2, metric="perc")

#extract interesting TXFs and plot against all tissue types
counts_2 <- hoRse_rrtlog_extract(counts, N=4)
hoRse_plot_counts(counts_2, metric="perc")

##'-------------------------------------------------------------------
# genenames     <- c("ACAN", "SOX9", "COL2A1", "LCN2")
# illuminanames <- c("ILMN_2326509", "ILMN_2147517", "ILMN_1712678", "ILMN_2224103")
chrIn         <- c("15:88803443:88875354:1")
ensembl <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")

mapping <- ID2TXF(x=genenames, type="gene_name", mart=ensembl, hoRse.environment=hoRse.env)
mapping[mapping$Gene == "LCN2",]

chrIn         <- c("15:88803443:88875354:1")
mapping <- ID2TXF(x=chrIn, type="chromosomal", mart=ensembl, hoRse.environment=hoRse.env)



hoRse_plot_counts_line(counts[counts$txf %in% mapping$txf,], metric="rrtlog")
hoRse_plot_counts(counts[counts$txf %in% mapping$txf,], metric="perc")
foo <- hoRse_rrtlog_extract(counts[counts$txf %in% mapping$txf,], N=4)

hoRse_plot_counts_line(foo, metric="rrtlog")
hoRse_plot_counts(foo, metric="perc")

counts
##'-------------------------------------------------------------------

