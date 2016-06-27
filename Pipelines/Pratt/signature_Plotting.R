data_in          <- batchCorrected_data
data_in_RA       <- batchCorrected_data[, (pData(lumi.Q)$Pool_ID == "RA")]
data_in_NRA      <- batchCorrected_data[, (pData(lumi.Q)$Pool_ID == "NRA")]


RA_SD            <- apply(data_in_RA, 1,sd)
NRA_SD           <- apply(data_in_NRA,1,sd)
plot(density(RA_SD))
plot(density(NRA_SD))
abline(v=0.2, col=2)

RA_SD_filtered   <- RA_SD[RA_SD < 0.2]
NRA_SD_filtered  <- NRA_SD[NRA_SD < 0.2]

SD_Filtered      <- intersect(names(RA_SD_filtered), names(NRA_SD_filtered))
SD_Filtered_set  <- batchCorrected_data[intersect(SD_Filtered, rownames(batchCorrected_data)),]

plot(density(apply(SD_Filtered_set,1,sd)))

RA_Mean         <- apply(SD_Filtered_set[, grep('_RA_', colnames(SD_Filtered_set))], 1, mean)
NRA_Mean        <- apply(SD_Filtered_set[, grep('_NRA_', colnames(SD_Filtered_set))], 1, mean)

Filtered_Diff   <- NRA_Mean - RA_Mean
max(Filtered_Diff)
min(Filtered_Diff)
Test_Set <- names(sort(Filtered_Diff))[1:10]
Test_Set <- c(Test_Set, names(sort(Filtered_Diff, decreasing=T))[1:10])

Test_Set_In <- data_in[Test_Set, ]

##############################All In#############################################################################
RA_Mean         <- apply(data_in_RA, 1, mean)
NRA_Mean        <- apply(data_in_NRA, 1, mean)

Filtered_Diff   <- NRA_Mean - RA_Mean
max(Filtered_Diff)
min(Filtered_Diff)
Test_Set <- names(sort(Filtered_Diff))[1:10]
Test_Set <- c(Test_Set, names(sort(Filtered_Diff, decreasing=T))[1:10])

Test_Set_In <- data_in[Test_Set, ]
##############################All In#############################################################################


##############################Dual Plot#############################################################################
Test_Set_In <- data_in[rownames(gene_list_filtered)[1:10],]
pheno       <- pData(lumi.Q)
fooRA       <- batchCorrected_data[selection_IDs, (pData(lumi.Q)$Pool_ID == "RA")]
fooNRA      <- batchCorrected_data[selection_IDs, (pData(lumi.Q)$Pool_ID == "NRA")]
fooUA       <- batchCorrected_data[selection_IDs, (pData(lumi.Q)$Pool_ID == "UA")]
axis_val    <- selection_unique$symbol
axis_val[13:14] <- c('GPRIN3', 'LOC731186')
for(i in 1:nrow(fooRA)) {
  png(file=paste("exp_plot_", rownames(fooRA)[i], ".png", sep=""),width=12.53,height=6.98,units="in",res=600) 
  par( mfrow = c(1, 3))
  
  win_max  <- max(c(as.vector(fooRA[i,]), as.vector(fooNRA[i,]), as.vector(fooUA[i,])))
  win_min  <- min(c(as.vector(fooRA[i,]), as.vector(fooNRA[i,]), as.vector(fooUA[i,])))
  
  plot(x=c(1:ncol(fooRA)), y=fooRA[i,], xlab="Patient", ylab="log2(Expression)", col=2, #main=paste("RA Samples - ", rownames(fooRA)[i], sep=""),
       ylim=c(win_min, win_max))
  lines(x=c(1:length(colnames(fooRA))), y=fooRA[i,], col=2)
  abline(h=mean(fooRA[i,]), col=1)
  text(x=c(1:ncol(fooRA)), y=fooRA[i,], labels=pheno[pheno$Pool_ID == "RA",]$Sex, cex= 0.7, pos=3)
  
  plot(x=c(1:ncol(fooNRA)), y=fooNRA[i,], main=paste("NRA Samples - ", rownames(fooNRA)[i], " (", axis_val[i], ")",  sep=""), xlab="Patient", ylab="log2(Expression)", col=3,
       ylim=c(win_min, win_max))
  lines(x=c(1:length(colnames(fooNRA))), y=fooNRA[i,], col=3)
  abline(h=mean(fooNRA[i,]), col=1)
  text(x=c(1:ncol(fooNRA)), y=fooNRA[i,], labels=pheno[pheno$Pool_ID == "NRA",]$Sex, cex= 0.7, pos=3)
  
  plot(x=c(1:ncol(fooUA)), y=fooUA[i,], xlab="Patient", ylab="log2(Expression)", col=4, #main=paste("NRA Samples - ", rownames(fooUA)[i], sep=""),
       ylim=c(win_min, win_max))
  lines(x=c(1:length(colnames(fooUA))), y=fooUA[i,], col=4)
  abline(h=mean(fooUA[i,]), col=1)
  text(x=c(1:ncol(fooUA)), y=fooUA[i,], labels=pheno[pheno$Pool_ID == "UA",]$Sex, cex= 0.7, pos=3)
  
  dev.off()
#   cat ("Press [enter] to continue")
#   line <- readline()
}
##################################################################################################################


##############################Single Plot#############################################################################
# Test_Set_In <- data_in[rownames(gene_list_filtered)[1:10],]
pheno       <- pData(lumi.Q)
fooRA       <- batchCorrected_data[selection_IDs, (pData(lumi.Q)$Pool_ID == "RA")]
fooNRA      <- batchCorrected_data[selection_IDs, (pData(lumi.Q)$Pool_ID == "NRA")]
co_var      <- grep("Sex", colnames(pheno))
for(i in 1:nrow(fooRA)) {
  #   png(file=paste("exp_plot_", rownames(fooRA)[i], ".png", sep=""),width=12.53,height=6.98,units="in",res=600) 
  
  win_max  <- max(c(as.vector(fooRA[i,]), as.vector(fooNRA[i,])))
  win_min  <- min(c(as.vector(fooRA[i,]), as.vector(fooNRA[i,])))
  
  if(ncol(fooRA) > ncol(fooNRA)) {x_lim <- ncol(fooRA)} else {x_lim <- ncol(fooNRA)}
  
  plot(NULL, main=paste("RA vs NRA Samples Biosignature - ", rownames(fooRA)[i], sep=""), xlab="Patient", ylab="log2(Expression)",
       ylim=c(win_min, win_max), xlim=c(1, x_lim), xaxt="n")
  
  points(x=c(1:ncol(fooRA)), y=fooRA[i,], col=2)
  lines(x=c(1:length(colnames(fooRA))), y=fooRA[i,], col=2, lwd=2.5)
  abline(h=mean(fooRA[i,]), col=2)
  text(x=c(1:ncol(fooRA)), y=fooRA[i,], labels=pheno[pheno$Pool_ID == "RA",][,co_var], cex= 0.7, pos=3)
  
  points(x=c(1:ncol(fooNRA)), y=fooNRA[i,], col=3)
  lines(x=c(1:length(colnames(fooNRA))), y=fooNRA[i,], col=3, lwd=2.5)
  abline(h=mean(fooNRA[i,]), col=3)
  text(x=c(1:ncol(fooNRA)), y=fooNRA[i,], labels=pheno[pheno$Pool_ID == "NRA",][,co_var], cex= 0.7, pos=3, offset=.5)
  
  legend('topright', c("RA","Non-RA", paste("RA  - Mean:", mean(fooRA[i,]), "SD:", sd(fooRA[i,]), sep=" "),
                       paste("NRA - Mean:", mean(fooNRA[i,]), "SD:", sd(fooNRA[i,]), sep=" "), 
                       paste("Text Overlay:", colnames(pheno)[co_var], sep=" ")), 
         lty=c(1,1,NA,NA,NA), lwd=c(2.5,2.5,NA,NA,NA), col=c(2,3,1,1,1), cex=.6)
  
  #   dev.off()
  cat ("Press [enter] to continue")
  line <- readline()
}
##################################################################################################################

pheno       <- pData(lumi.Q)
co_var      <- grep("Sex", colnames(pheno))

plotPCA(normalised_data, groups=as.numeric(as.factor(pData(lumi.Q)[,co_var])), 
        groupnames=levels(as.factor(pData(lumi.Q)[,co_var])), pcs=c(1,2), 
        plot3d=F, main="Grouped By Treatment - Post Batch Correction", addtext=pData(lumi.Q)$Sample_Name)

plotPCA(normalised_data, groups=as.numeric(as.factor(pData(lumi.Q)[,co_var])), 
        groupnames=levels(as.factor(pData(lumi.Q)[,co_var])), pcs=c(1,2), 
        plot3d=F, main="Grouped By Treatment - Post Batch Correction")#, addtext=pData(lumi.Q)$Sex)



##'GENE SIGNATURE PLOT---------------------------------------------------------------------------------------#

png("Twelve_Gene_Suspected_Signature.png",width=12.53,height=6.98,units="in",res=600)
# foo      <- batchCorrected_data[selection_IDs, (pData(lumi.Q)$Pool_ID == "NRA")]
foo      <- normalised_data[selection_IDs,]# (pData(lumi.Q)$Pool_ID == "UA")]
plot(NULL, xlim=c(1, length(rownames(foo))), ylim=c(6, 14), ylab="log2(Expression)", xlab="Genes", main="Expression Plot of Suspected 12-Gene RA Signature",
     xaxt="n")
legend('topright', c("RA","Non-RA"), lty=c(1,1), lwd=c(2.5,2.5), col=c(2,3), cex=.6)
for(i in 1:length(colnames(foo))) {
  points(x=c(1:length(rownames(foo))), y=foo[,i], col=3)
  lines(x=c(1:length(rownames(foo))), y=foo[,i], col=3)
}
foo       <- batchCorrected_data[selection_IDs, (pData(lumi.Q)$Pool_ID == "RA")]
for(i in 1:length(colnames(foo))) {
  points(x=c(1:length(rownames(foo))), y=foo[,i], col=2)
  lines(x=c(1:length(rownames(foo))), y=foo[,i], col=2)
}
axis_val <- selection_unique$symbol
axis_val[12:13] <- c('GPRIN3', 'LOC731186')
axis(1, at=1:nrow(foo),labels=axis_val, las=2, cex.axis=.7)
dev.off()
##'----------------------------------------------------------------------------------------------------------#


##'GENE SIGNATURE PLOT UA------------------------------------------------------------------------------------#

# png("Twelve_Gene_Suspected_Signature.png",width=12.53,height=6.98,units="in",res=600)
# foo      <- batchCorrected_data[selection_IDs, (pData(lumi.Q)$Pool_ID == "NRA")]
foo      <- batchCorrected_data[selection_IDs, (pData(lumi.Q)$Pool_ID == "RA")]
plot(NULL, xlim=c(1, length(rownames(foo))), ylim=c(6, 14), ylab="log2(Expression)", xlab="Genes", main="Expression Plot of Suspected 12-Gene RA Signature",
     xaxt="n")
legend('topright', c("RA","Non-RA"), lty=c(1,1), lwd=c(2.5,2.5), col=c(2,3), cex=.6)
for(i in 1:length(colnames(foo))) {
  points(x=c(1:length(rownames(foo))), y=foo[,i], col=3)
  lines(x=c(1:length(rownames(foo))), y=foo[,i], col=3)
}
foo       <- batchCorrected_data[selection_IDs, (pData(lumi.Q)$Pool_ID == "RA")]
for(i in 1:length(colnames(foo))) {
  points(x=c(1:length(rownames(foo))), y=foo[,i], col=2)
  lines(x=c(1:length(rownames(foo))), y=foo[,i], col=2)
}
axis_val <- selection_unique$symbol
axis_val[13:14] <- c('GPRIN3', 'LOC731186')
# axis_val[12] <- c('LOC731186')
axis(1, at=1:nrow(foo),labels=axis_val, las=2, cex.axis=.7)

legend('topright', c("RA","Non-RA","UA"), lty=c(1,1,1), lwd=c(2.5,2.5,2.5), col=c(3,4,5), cex=.6)
# dev.off()
##'----------------------------------------------------------------------------------------------------------#

##'GENE SIGNATURE DISTRIBUTION PREDICTION--------------------------------------------------------------------#

fooRA       <- batchCorrected_data[gene_list_filtered$probe_list[1:5], (pData(lumi.Q)$Pool_ID == "RA")]
fooNRA      <- batchCorrected_data[gene_list_filtered$probe_list[1:5], (pData(lumi.Q)$Pool_ID == "NRA")]
fooUA       <- batchCorrected_data[gene_list_filtered$probe_list[1:5], (pData(lumi.Q)$Pool_ID == "UA")]
DE_IDs      <- gene_list_filtered$probe_list[1:5]





##'----------------------------------------------------------------------------------------------------------#