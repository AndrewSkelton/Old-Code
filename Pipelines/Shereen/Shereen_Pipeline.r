#-------------------------------------------------------------------------------------------#
#  Author      : Andrew J Skelton                                                           |
#  Language    : R Statistical Programming Language                                         |
#  Optimised   : Sjogren's Syndome                                                          |
#  Data Owner  : Newcastle University - Fai                                                 |
#  Description : A Microarray analysis pipeline for downstream analysis of a                |
#                CSV file                                                                   |
#-------------------------------------------------------------------------------------------#

source("http://bioconductor.org/biocLite.R")
biocLite()
biocLite(c("lumi", "gplots", "ggplot2", "limma", "ExiMiR"))
library(stringr)
library(sva)
library(lumi)
library(gplots)
library(ggplot2)
library(annotate)
library(lumiHumanAll.db)
library(limma)
library(pathview)
library(GOstats)
library(RamiGO)

##LOAD_DATA_&_PREPARE----------------------------------------------------------------------------------------#
setwd("/Users/andrew/Documents/Bioinformatics/Customers/Shereen/Fai_Samples")
setwd("~/Documents/Bioinformatics/Customers/Shereen/")
list.files()
filename <- "Fai_Sample_Probe_Profile.txt"

raw_data <- lumiR(filename)
sampleNames(raw_data)
raw_data <- raw_data[, -c(2, 6)]
raw_data <- raw_data[, -c(15, 52, 69, 125)]

array_names <- c("NCL-123-0_CTRL_1","CAM-006-0_CTRL_1","DER-006-1_SS_1",  "BAT-023-0_CTRL_1","NCL-084-1_SS_1",
                 "NCL-117-0_CTRL_1","BIR-051-1_SS_1",  "BAS-017-0_CTRL_1","NCL-052-1_SS_1",  "GAT-027-0_CTRL_1",
                 "NCL-113-0_CTRL_1","NCL-011-1_SS_1",  "NCL-055-1_SS_1",  "NCL-091-0_CTRL_1","NCL-024-1_SS_1",
                 "SUN-009-1_SS_1",  "NCL-130-0_CTRL_1","BIR-029-1_SS_1",  "NCL-083-1_SS_1",  "NCL-136-0_CTRL_1",
                 "NCL-060-1_SS_1",  "FIF-009-1_L_4",   "BIR-073-1_PM_4",  "NCL-054-1_SS_4",  "LEE-062-1_PM_4",
                 "UCH-20-1_L_4",    "SUR-001-1_L_4",   "NOT-019-1_PM_4",  "BIR-037-1_SS_4",  "POR-009-0_CTRL_4",
                 "BLH-006-1_C_4",   "LEE-048-1_SS_4",  "UCH-017-1_L_4",   "SWI-110-0_CTRL_4","NOT-015-1_PM_4",
                 "NCL-165-0_CTRL_4","NCL-031-0_CTRL_4","UCH-016-1_PM_4",  "BIR-039-1_C_4",   "NCL-007-1_PM_4",
                 "SWI-084-0_CTRL_4","BLH-021-1_PM_4",  "SWI-039-1_C_4",   "WIN-007-1_L_4",   "LEE-034-1_C_4",
                 "NCL-191-0_CTRL_2","NCL-103-1_C_2",   "TOR-011-0_CTRL_2","NCL-164-1_C_2",   "BAS-011-1_PM_2",
                 "BAS-020-0_CTRL_2","BIR-042-1_SS_2",  "GAT-005-1_C_2",   "NCL-078-1_PM_2",  "NCL-089-1_PM_2",
                 "NCL-081-1_PM_2",  "BIR-082-1_C_4",   "BLH-017-1_SS_4",  "BIR-024-1_SS_4",  "BIR-036-1_PM_4",
                 "TOR-007-0_CTRL_4","BIR-019-1_SS_4",  "BLH-019-1_SS_4",  "SWI-078-1_PM_4",  "TOR-010-0_CTRL_4",
                 "BIR-033-1_SS_4",  "LEE-058-0_CTRL_4","NCL-005-1_L_4",   "SWI-008-1_C_4",   "SWI-081-1_PM_4",
                 "DER-002-1_SS_4",  "BIR-043-1_SS_4",  "BIR-027-1_SS_4",  "GLA-029-1_C_4",   "NCL-097-0_CTRL_4",
                 "FIF-019-1_C_4",   "BIR-055-1_SS_4",  "FIF-018-1_L_4",   "WIN-023-1_PM_4",  "NOT-036-1_L_3",
                 "UCH-002-1_L_3",   "BIR-041-1_SS_3",  "NCL-044-1_SS_3",  "NCL-018-1_L_3",   "SWI-005-1_SS_3",
                 "GAT-002-1_L_3",   "NCL-053-1_SS_3",  "SWI-049-1_SS_3",  "GLA-015-1_L_3",   "IPS-004-1_L_3",
                 "SWI-031-1_SS_3",  "NCL-025-1_L_3",   "NCL-067-1_SS_3",  "FIF-027-1_L_3",   "NCL-100-0_CTRL_3",
                 "NOT-012-1_L_3",   "BIR-018-1_C_3",   "NCL-030-1_SS_3",  "LEE-018-1_SS_3",  "LEE-002-1_SS_3",
                 "BIR-025-1_SS_3",  "DER-004-1_SS_3",  "SOU-007-1_PM_3",  "SWI-012-1_SS_3",  "GLA-026-1_C_3",
                 "BIR-053-1_SS_3",  "GAT-010-1_PM_3",  "WIN-011-0_CTRL_3","BIR-034-1_SS_3",  "SWI-038-1_SS_3",
                 "SOU-013-1_PM_3",  "NCL-029-1_SS_3",  "BLH-003-1_SS_3",  "NCL-035-1_SS_3",  "NCL-038-0_CTRL_3",
                 "SWI-006-1_C_3",   "NCL-010-1_SS_3",  "NCL-061-1_SS_3",  "BIR-032-1_SS_3",  "FIF-014-1_C_3",
                 "CAM-008-0_CTRL_3","NCL-028-1_SS_3",  "LEE-052-1_SS_3",  "SWI-073-1_C_3",   "BIR-060-0_CTRL_3",
                 "BIR-011-1_SS_3",  "BIR-021-1_SS_3",  "NCL-056-1_SS_3",  "NCL-009-1_SS_3",  "TOR-012-0_CTRL_3",
                 "SUN-011-1_SS_3",  "NCL-101-0_CTRL_3","NOT-017-1_C_3",   "NCL-096-1_PM_3",  "BIR-022-1_SS_3",
                 "SWI-072-1_PM_3",  "NCL-033-1_SS_3",  "NCL-039-1_SS_3",  "SWI-085-0_CTRL_3","NCL-014-1_SS_3",
                 "GLA-019-1_SS_3",  "BIR-006-0_CTRL_3","DER-013-1_SS_3",  "LEE-016-1_SS_3",  "DER-019-1_C_3",
                 "NCL-036-1_SS_3",  "BIR-030-1_SS_3",  "BIR-031-1_SS_3",  "BIR-013-1_SS_3",  "WIN-016-1_C_3",
                 "LEE-044-0_CTRL_3","BIR-001-1_SS_3",  "NCL-034-1_SS_3",  "LEE-041-1_SS_3",  "LEE-012-1_SS_3",
                 "BLH-018-1_SS_3",  "BIR-008-1_PM_3",  "LEE-015-1_SS_3",  "BIR-017-1_SS_3",  "SWI-060-1_C_3",
                 "BIR-046-1_SS_3",  "HAR-004-1_PM_3")

sampleNames(raw_data) <- array_names
##'----------------------------------------------------------------------------------------------------------#

##'----------------------------------------------------------------------------------------------------------#
##'RIN SCORE EXCLUSION
##'
# Arrays Less than 7
Seven <- c("LEE-060-0","BIR-033-1","NCL-097-0","NCL-007-1","TOR-006-1","BIR-039-1","LEE-034-1", 
           "LEE-062-1","TOR-007-0","NOT-036-1","FIF-027-1","BIR-011-1","SWI-006-1","BIR-030-1", 
           "IPS-002-1","DER-019-1","FIF-014-1","WIN-016-1","LEE-016-1","LEE-012-1","GLA-019-1", 
           "BIR-041-1","SWI-031-1","NCL-053-1")

matches <- unique (grep(paste(Seven,collapse="|"), 
                        colnames(raw_data)))

# Arrays Less than 7 with lymphoma
SevenL <- c("LEE-060-0", "BIR-033-1", "NCL-097-0", "NCL-007-1", "TOR-006-1", "BIR-039-1", "LEE-034-1", 
            "LEE-062-1", "TOR-007-0", "BIR-011-1", "SWI-006-1", "BIR-030-1", 
            "IPS-002-1", "DER-019-1", "FIF-014-1", "WIN-016-1", "LEE-016-1", "LEE-012-1", "GLA-019-1", 
            "BIR-041-1", "SWI-031-1", "NCL-053-1")

matches <- unique (grep(paste(SevenL,collapse="|"), 
                        colnames(raw_data)))

# Arrays Less than 5
#Five <- c("NCL-097-0", "BIR-011-1", "SWI-006-1", "BIR-041-1", "SWI-031-1", "NCL-053-1")

#matches <- unique (grep(paste(Five,collapse="|"), 
#                       colnames(raw_data)))

#Eliminate Arrays
length(colnames(raw_data))
raw_data_L <- raw_data[,-matches]
length(colnames(raw_data))

##'----------------------------------------------------------------------------------------------------------#


##'SUBSET BATCHES OR NOT??-----------------------------------------------------------------------------------#
############################
##'TRUE  - USE BATCHES 2-4 #
##'FALSE - USE BATCHES 1-4 #
      subset <- TRUE       #
############################

if(subset == TRUE)
{
    raw_data_in <- raw_data[,grep('_2|_3|_4', sampleNames(raw_data))]
}else
{
    raw_data_in <- raw_data
}
##'----------------------------------------------------------------------------------------------------------#

##'NORMALISATION---------------------------------------------------------------------------------------------#
vst_data         <- lumiT(raw_data_in)
rsn_data         <- lumiN(vst_data, method = "rsn")
normalised_data  <- lumiQ(rsn_data)
##'QUANTILE_NORMALISATION - COMMENT OUT RSN NORMALISATION LINE IF YOU'RE USING QUANTILE
#quantile_data       <- lumiN(vst_data, method="quantile")
#normalised_data     <- lumiQ(quantile_data)
##'----------------------------------------------------------------------------------------------------------#

##'DETECTION_CALL_THRESHOLD_REMOVAL--------------------------------------------------------------------------#
exprs_data                   <- exprs(normalised_data)
present_count                <- detectionCall(normalised_data)
normalised_data              <- exprs_data[present_count > 0, ]
##'----------------------------------------------------------------------------------------------------------#

##'BATCH CORRECTION (PHASE)----------------------------------------------------------------------------------#
batches <- colnames(normalised_data)
if(subset == TRUE)
{
    batches[grep('_2', batches)] <- 1
    batches[grep('_3', batches)] <- 2
    batches[grep('_4', batches)] <- 3
}else
{
    batches[grep('_1', batches)] <- 1
    batches[grep('_2', batches)] <- 2
    batches[grep('_3', batches)] <- 3
    batches[grep('_4', batches)] <- 4
}

classification   <- sapply(strsplit(colnames(normalised_data), split='_'), function(x) x[[2]])
pheno            <- data.frame(sample = c(1:length(colnames(normalised_data))),
                               outcome = classification, batch = batches)
rownames(pheno)  <- colnames(normalised_data)
batch            <- pheno$batch
mod              <- model.matrix(~as.factor(outcome), data = pheno)
batchCorrected_data = ComBat(dat = normalised_data, batch = batch, mod = mod, numCovs = NULL,
                                      par.prior = TRUE, prior.plots = TRUE)
##'----------------------------------------------------------------------------------------------------------#

##PROBE_ANNOTATION-------------------------------------------------------------------------------------------#
probe_list       <- rownames(batchCorrected_data)
nuIDs            <- probeID2nuID(probe_list)[, "nuID"]
symbol           <- getSYMBOL(nuIDs, "lumiHumanAll.db")
name             <- unlist(lookUp(nuIDs, "lumiHumanAll.db", "GENENAME"))
anno_df          <- data.frame(ID = nuIDs, probe_list, symbol, name)
##'----------------------------------------------------------------------------------------------------------#

##EXPERIMENTAL_DESIGN----------------------------------------------------------------------------------------#
treatments  <- c("CTRL", "SS", "L", "C", "PM")
array_names <- colnames(normalised_data)
array_names[grep('_CTRL_', array_names)] <- "CTRL"
array_names[grep('_SS_',   array_names)] <- "SS"
array_names[grep('_L_',    array_names)] <- "L"
array_names[grep('_C_',    array_names)] <- "C"
array_names[grep('_PM_',   array_names)] <- "PM"

batchCorrected_data <- data.matrix(batchCorrected_data)
design              <- model.matrix(~0 + factor(array_names, levels = treatments))
colnames(design)    <- treatments
num_parameters      <- ncol(design)
fit                 <- lmFit(batchCorrected_data, design)

cont_mat         <- makeContrasts(CTRL-SS, CTRL-L, CTRL-C, CTRL-PM,
                                  SS-L,    SS-C,   SS-PM,
                                  L-C,     L-PM,
                                  C-PM,
                                  levels=treatments)
fit2             <- contrasts.fit(fit, contrasts=cont_mat)
fit2             <- eBayes(fit2)
fit2$genes = anno_df
##'----------------------------------------------------------------------------------------------------------#

##GENE_LISTS-------------------------------------------------------------------------------------------------#
comparisons <- c("CTRL - SS", "CTRL - L", "CTRL - C", "CTRL - PM",
                 "SS - L",    "SS - C",   "SS - PM",  "L - C",
                 "L - PM",    "C - PM")
p_cut_off   <- 0.05
fold_change <- 1.2
i           <- 1

#pdf("volcano_plots.pdf")#
#for(i in 1:length(comparisons))#
#{#
topTable(fit2, coef = comparisons[i], p.value = p_cut_off, lfc = log2(fold_change), adjust.method="BH")
gene_list_unfiltered <- topTable(fit2, coef = comparisons[i], number = nrow(anno_df), adjust.method="BH")
gene_list            <- topTable(fit2, coef = comparisons[i], p.value = p_cut_off, lfc = log2(fold_change),
                                 number = nrow(anno_df), adjust.method = "BH")

comparison_split <- strsplit(comparisons[i], split=" ")
controls      <- batchCorrected_data[as.character(gene_list$probe_list),
                                     grep(comparison_split[[1]][1], array_names)]
conditions    <- batchCorrected_data[as.character(gene_list $probe_list),
                                     grep(comparison_split[[1]][3], array_names)]

gene_list$ControlMean   <- apply(controls,   1, mean)
gene_list$ConditionMean <- apply(conditions, 1, mean)

write.csv(gene_list, paste("gene_list_", comparison_split[[1]][1], "-", comparison_split[[1]][3], ".csv",
                           sep="",collapse=""))
##'----------------------------------------------------------------------------------------------------------#

##'VOLCANO_PLOT----------------------------------------------------------------------------------------------#
gene_list_unfiltered$threshold = as.factor(abs(gene_list_unfiltered$logFC) > log2(fold_change) &
                                             gene_list_unfiltered$adj.P.Val < p_cut_off)
g = ggplot(data=gene_list_unfiltered, aes(x=gene_list_unfiltered$logFC,
                                          y=-log10(gene_list_unfiltered$P.Value),
                                          colour=gene_list_unfiltered$threshold)) +
  geom_point(alpha=0.4, size=1.75) +
  theme(legend.position = "none") +
  #xlim(c(-4, 4)) + ylim(c(0, 4.5)) +
  xlab("log2 fold change") + ylab("-log10 p-value") +
  ggtitle(paste("Volcano Plot of Fold Change vs P-Value", comparisons[i], sep=" "))
print(g)
#}#
#dev.off()#
##'----------------------------------------------------------------------------------------------------------#

##'PATHWAY_ANALYSIS----------------------------------------------------------------------------------------------#
##'TO BE REPLACED BY REACTOME.
sig_probes      <- as.character(gene_list$probe_list)
entrez          <- unique(unlist(lookUp(nuIDs[sig_probes],
                                        "lumiHumanAll.db", "ENTREZID")))
entrez          <- as.character(entrez[!is.na(entrez)])
entrez_universe <- unique(unlist(lookUp(nuIDs,"lumiHumanAll.db", "ENTREZID")))
entrez_universe <- as.character(entrez_universe[!is.na(entrez_universe)])

kegg_params <- new("KEGGHyperGParams",
                   geneIds=entrez,
                   universeGeneIds=entrez_universe,
                   annotation="lumiHumanAll.db",
                   pvalueCutoff= 0.01,
                   testDirection="over")
kegg_hyperg_result <- hyperGTest(kegg_params)
print(kegg_hyperg_result)

sig_logFCs <- gene_list$logFC
fc_entrez <- unlist(lookUp(nuIDs[sig_probes],
                           "lumiHumanAll.db", "ENTREZID"))
names(sig_logFCs) <- fc_entrez
sig_logFCs <- sig_logFCs[!is.na(entrez)]

summary(kegg_hyperg_result)

for(j in 1:length(sigCategories(kegg_hyperg_result)))
{
  pv <- pathview(gene.data=sig_logFCs,
                 pathway.id=sigCategories(kegg_hyperg_result)[j],
                 species="hsa",
                 limit=list(gene=c(round(min(sig_logFCs)),
                                   round(max(sig_logFCs))), cpd=1))
}
##'----------------------------------------------------------------------------------------------------------#

##GO ANNOTATION----------------------------------------------------------------------------------------------#
sig_probes      <- as.character(gene_list$probe_list)
entrez          <- unique(unlist(lookUp(nuIDs[sig_probes],
                                        "lumiHumanAll.db", "ENTREZID")))
entrez          <- as.character(entrez[!is.na(entrez)])
entrez_universe <- unique(unlist(lookUp(nuIDs,"lumiHumanAll.db", "ENTREZID")))
entrez_universe <- as.character(entrez_universe[!is.na(entrez_universe)])

params <- new("GOHyperGParams", geneIds=entrez,universeGeneIds=entrez_universe,annotation="lumiHumanAll.db",
              ontology="BP", pvalueCutoff= 0.01, conditional=FALSE, testDirection="over")
hyperg_result <- hyperGTest(params)
print(hyperg_result)

pval_go     <- pvalues(hyperg_result)
go_fdr      <- p.adjust(pvalues(hyperg_result), method="fdr")
sig_go_id   <- names(go_fdr[go_fdr < 0.01])
if(length(sig_go_id) > 0)
{
  sig_go_term <- getGOTerm(sig_go_id)#[["BP"]]
  GO_Table <- data.frame(GO.ID = sig_go_id, Adj.P.Val = go_fdr[go_fdr < 0.01], GO.Term = unlist(sig_go_term))
  length(GO_Table$GO.ID)
  write.table(GO_Table, paste(comparisons[i], "_P0.01_GO.txt", sep=""), sep="\t")
  amigo_tree = getAmigoTree(sig_go_id, pvalues=go_fdr[go_fdr < 0.01], webserver="http://amigo1.geneontology.org/cgi-bin/amigo/visualize",
                          filename=paste(comparisons[i], "_P0.001_upregulated_RamiGo.png", sep=""))
}
##'----------------------------------------------------------------------------------------------------------#
