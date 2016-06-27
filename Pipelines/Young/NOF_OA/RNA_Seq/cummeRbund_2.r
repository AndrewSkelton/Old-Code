myGenes<-getGenes(cuff)
myGenes

gene.rep.matrix<-repFpkmMatrix(genes(cuff))
gene.matrix.annotation<-annotation(genes(cuff))
ele <- c(4,5)
gene.rep.matrix <- cbind(gene.rep.matrix, gene.matrix.annotation[,ele])
#log FC && Q val
head(gene.rep.matrix)
   
isoform.rep.matrix<-repFpkmMatrix(isoforms(cuff))
isoform.matrix.annotation<-annotation(isoforms(cuff))
#log FC && Q val
ele <- c(2,4,6,8,9)
isoform.rep.matrix <- cbind(isoform.rep.matrix, isoform.matrix.annotation[,ele])
head(isoform.rep.matrix)

foo <- gene.rep.matrix[gene.rep.matrix[,1:6] > 1]

foo <- c()
for(i in 1:length(rownames(gene.rep.matrix)))
{
  z <- gene.rep.matrix[i,1:6] > 1
  if(length(z[z==TRUE]) > 5)
  {
    y <- gene.rep.matrix[i,7:16] > 1
    if(length(y[y==TRUE]) > 8)
    {
      foo <- rbind(foo, gene.rep.matrix[i,])
    }
  }
}

foo_1 <- c()
for(i in 1:length(rownames(isoform.rep.matrix)))
{
  z <- isoform.rep.matrix[i,1:6] > 1
  if(length(z[z==TRUE]) > 5)
  {
    y <- isoform.rep.matrix[i,7:16] > 1
    if(length(y[y==TRUE]) > 8)
    {
      foo_1 <- rbind(foo_1, isoform.rep.matrix[i,])
    }
  }
}

write.csv(foo, "FPKM_Cutoff_Genes.csv")
write.csv(foo_1, "FPKM_Cutoff_Isoforms.csv")

diff_exp <- -log10(foo)
diff_exp[diff_exp == Inf] <- 0 
heatmap.2(cor(diff_exp), col="redblue", density.info="none",
trace="none", cexRow=0.8, cexCol=0.8,
main="NOF Vs OA Differentially Expressed Genes")

d <- dist(cor(as.matrix(diff_exp)))   # find distance matrix 
hc <- hclust(d)                # apply hirarchical clustering 
plot(hc) 

heatmap.2(cor(-log10(foo[,1:16])), col=my_palette, density.info="none",
          trace="none", cexRow=0.8, cexCol=0.8,
          main="NOF Vs OA Differentially Expressed Genes")



gene.diff<-diffData(genes(cuff))
require(ggplot2)

for(i in 1:length(rownames(gene.diff)))
{
  gene.diff$t[i] <- "1"
  if( (abs(gene.diff$log2_fold_change[i] > 2)) || ((gene.diff$log2_fold_change[i] < -2)) )
  {
    gene.diff$threshold[i] <- "4"
  }
  if(gene.diff$q_value[i] < 0.05)
  {
    gene.diff$threshold[i] <- "3"
    #print(1)
  }
  if(abs(gene.diff$log2_fold_change[i]) > 2 & gene.diff$q_value[i] < 0.05)
  {
    gene.diff$threshold[i] <- "2"
  }
}


gene.diff$threshold <- "1"
gene.diff[(gene.diff$q_value < 0.05),]$threshold <- "2"
gene.diff[(gene.diff$log2_fold_change > 2) | (gene.diff$log2_fold_change < -2),]$threshold <- "3"
gene.diff[((gene.diff$log2_fold_change > 2) | (gene.diff$log2_fold_change < -2)) & (gene.diff$q_value < 0.05),]$threshold <- "4"
cb_palette <- c("#000000", "#70DB70", "#FF3333", "#CC33FF")

g = ggplot(data=gene.diff, aes(x=(log2_fold_change), y=-log10(q_value), colour=threshold)) +
  geom_point(alpha=0.1, size=1.75) +
  theme(legend.position = "none") +
  xlab("log2 fold change") + ylab("-log10 p-value") +
  labs(title="Volcano Plot showing significantly expressed genes") +
  theme(panel.background = element_blank()) +
  scale_color_manual(values=cb_palette) +
  geom_point(size = 1.2) #+
  #geom_text(aes(x=gene.diff$log2_fold_change, y=-log10(gene.diff$q_value),
  #              label=gene.diff$gene_id, size=1.2), colour="black")
g