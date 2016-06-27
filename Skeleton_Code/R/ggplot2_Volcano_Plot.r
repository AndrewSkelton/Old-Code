#=======================================================#
#	ggplot2 Volcano Plot Template						            #
#  														                         #
#	Description : Adds 'threshold' slot to input 		    #
#	dataframe with strings (for fixed scale). Uses	   	#
#	dynamic colour palette that can be customised.	   	#
#														                           #
#=======================================================#

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
