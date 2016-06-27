library(RColorBrewer)
palette <- brewer.pal(10, "Paired")

split         <- strsplit(colnames(lumi.Q),'_')
unique_arrays <- unique(grep('9477|9479|9481', unlist(split), value=T))

foo <- strsplit(colnames(lumi.Q), "_")
sample_ID  <- c()
Treatment   <- c()
Array_ID    <- c()
Array_Index <- c()
Array_Batch <- c()
Metho_State <- c()
for(i in 1:length(foo)) {
  sample_ID  <- c(sample_ID, unlist(foo[[i]])[1])
  Treatment   <- c(Treatment, unlist(foo[[i]])[2])
  Array_ID    <- c(Array_ID, unlist(foo[[i]])[3])
  Array_Index <- c(Array_Index, unlist(foo[[i]])[4])
  Array_Batch <- c(Array_Batch, unlist(foo[[i]])[5])
  Metho_State <- c(Metho_State, unlist(foo[[i]])[6])
}



png("density_Plot_Without879|948.png",width=12.53,height=6.98,units="in",res=600) 
for(i in 1:length(unique_arrays)) {
  range_in <- grep(unique_arrays[i], colnames(lumi.Q))
  plotCDF(lumi.Q[,range_in], reverse=TRUE, main=paste("Array ", unique_arrays[i], " (", i, ")", sep=""), col=brewer.pal(10, "Set3"))
  density(lumi.Q[,range_in])
}
dev.off()

range_in <- grep(unique_arrays[13], colnames(lumi.Q))
density(lumi.Q[,range_in])
plot(lumi.Q[,range_in], what='cv')
plotCDF(lumi.Q[,range_in], reverse=TRUE, main=paste("Array ", unique_arrays[i], " (", i, ")", sep=""), col=brewer.pal(10, "Set3"))



png("PCA_ArrayID_Grouping_Normalised.png",width=12.53,height=6.98,units="in",res=600) 
# svg("PCA_ArrayID_Grouping.svg") exprs(lumi.Q)
plotPCA(log2(normalised_data), groups=as.numeric(as.factor(Array_ID)), 
        groupnames=levels(as.factor(Array_ID)), pcs=c(1,2), 
        plot3d=F, main="Grouped By Array ID", squarepca=F)#, addtext=sample_ID)#xlim=c(-16,16), ylim=c(-14,12))
dev.off()

png("PCA_ScanBatch_Grouping.png",width=12.53,height=6.98,units="in",res=600) 
plotPCA(log2(exprs(lumi.Q)), groups=as.numeric(as.factor(Array_Batch)), 
        groupnames=levels(as.factor(Array_Batch)), addtext=sample_ID, pcs=c(1,2), 
        plot3d=F, main="Grouped By Array Batch")
dev.off()

png("PCA_ArrayIndex_Grouping_Normalised.png",width=12.53,height=6.98,units="in",res=600) 
plotPCA(log2(exprs(lumi.Q)), groups=as.numeric(as.factor(Array_Index)), 
        groupnames=levels(as.factor(Array_Index)), addtext=sample_ID, pcs=c(1,2), 
        plot3d=F, main="Grouped By Array Index")
dev.off()

pdf("PCA_Treatment_Grouping.pdf", width=14, height=6)
plotPCA(log2(exprs(lumi.Q)), groups=as.numeric(as.factor(Treatment)), 
        groupnames=levels(as.factor(Treatment)), pcs=c(1,2), 
        plot3d=F, main="Grouped By Treatment")
dev.off()

data_in   <- grep("ND", Metho_State, invert=T, value=F)
data_in_T <- grep("ND", Metho_State, invert=T, value=T)
pdf("PCA_Methotrexate_Grouping.pdf", width=14, height=6)
plotPCA(log2(exprs(lumi.Q[, data_in])), groups=as.numeric(as.factor(data_in_T)), 
        groupnames=levels(as.factor(data_in_T)), addtext=sample_ID[data_in], pcs=c(1,2), 
        plot3d=F, main="Grouped By Methotrexate Treatment")
dev.off()


plot(lumi.Q, what="sampleRelation", method="mds")
trans <- plotVST(vst_data)
matplot(log2(trans$untransformed), trans$transformed)





range_in <- grep(unique_arrays[12], colnames(lumi.Q))


boxplot(lumi.Q[,range_in])

plotPCA(exprs(lumi.Q[,range_in]), 
        addtext=colnames(lumi.Q[,range_in]), pcs=c(1,2,3), 
        plot3d=TRUE)

plotPCA(exprs(lumi.Q[,range_in]), 
        addtext=colnames(lumi.Q[,range_in]), pcs=c(1,2), 
        plot3d=F, legend=F)



plotPCA(exprs(lumi.Q), 
        addtext=colnames(lumi.Q), pcs=c(1,2), 
        plot3d=F, legend=F)


for(i in 1:length(unique_arrays)) {
  range_in <- grep(unique_arrays[i], colnames(lumi.Q))
  
  cat(unique_arrays[i])
  cat('\n')
  cat(Treatment[range_in])
  
  cat ("\nPress [enter] to continue")
  line <- readline()
}

