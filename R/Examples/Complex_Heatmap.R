#-------------------------------------------------------------------------------------------#
#  Author      : Andrew J Skelton                                                           |
#  Language    : R Statistical Programming Language                                         |
#  Study       : Complex Heatmap Example Code                                               |
#  Data Owner  : Utility Code                                                               |
#  Description : Custom Plots                                                               |
#-------------------------------------------------------------------------------------------#



##'Download and install
##'-----------------------------------------------------------------------------------------#
library(devtools)
install_github("jokergoo/ComplexHeatmap")

library(ComplexHeatmap)
library(circlize)
##'-----------------------------------------------------------------------------------------#



##'Reload the package if things get buggy
##'-----------------------------------------------------------------------------------------#
detach("package:ComplexHeatmap", unload=TRUE)
library(ComplexHeatmap)
library(circlize)
##'-----------------------------------------------------------------------------------------#


##'Matrix of Data in
##'-----------------------------------------------------------------------------------------#
signature_exp <- exprs(raw_data_det)[rownames(raw_data_det) %in% signature,
                                     grep("RA", pData(raw_data_det)$Update_Collapse_2)]
##'-----------------------------------------------------------------------------------------#



##'Heatmap.2 Flavouring
##'   - Dendrogram
##'   - Scale by Row
##'-----------------------------------------------------------------------------------------#
Rowv          <- colMeans(signature_exp)
d             <- dist(t(signature_exp))
h             <- hclust(d)
ddr           <- as.dendrogram(h)
ddr           <- reorder(ddr, Rowv)

rm            <- rowMeans(signature_exp)
signature_exp <- sweep(signature_exp, 1, rm)
sx            <- apply(signature_exp, 1, sd)
signature_exp <- sweep(signature_exp, 1, sx, "/")
##'-----------------------------------------------------------------------------------------#



##'Annotation and pin strip vector 
##'-----------------------------------------------------------------------------------------#
gene_names                    <- as.vector(anno_df[match(rownames(signature_exp),
                                                         anno_df$ID),]$symbol)
gene_names[is.na(gene_names)] <- "No Annotation"
rownames(signature_exp)       <- gene_names
foo                           <- as.vector(pData(raw_data_det)$Update_Collapse_2)
foo                           <- foo[grep("RA", foo)]
##'-----------------------------------------------------------------------------------------#



##'Complex Heatmap - Annotation
##'-----------------------------------------------------------------------------------------#
df  <- data.frame(Sample     = foo)
ha  <- HeatmapAnnotation(df  = df, 
                         col = list(Sample = c("RA"  =  "black", 
                                               "NRA" =  "grey")))
##'-----------------------------------------------------------------------------------------#



##'Complex Heatmap - Main Function Call
##'-----------------------------------------------------------------------------------------#
ht2 = Heatmap(signature_exp, 
              name                     = "Row Z-Score", 
              column_title             = "", 
              bottom_annotation        = ha,
              col                      = colorRamp2(c(-2, 0, 2), 
                                                    c("blue", "white", "red")),
              show_row_hclust          = F,
              cluster_rows             = F,
              show_column_hclust       = T,
              show_column_names        = F,
              cluster_columns          = ddr,
              column_hclust_height     = unit(1.5, "cm"),
              row_names_gp             = gpar(fontsize = 22)
)
##'-----------------------------------------------------------------------------------------#



##'Complex Heatmap - Draw
##'-----------------------------------------------------------------------------------------#
# png("Heatmap_Focus_V2.png", 
#     width=3096, 
#     height=2096, 
#     units="px", 
#     res=300)
draw(ht2, 
     show_heatmap_legend    = T, 
     show_annotation_legend = T,
     legend_grid_width      = unit(8, "mm"), 
     legend_title_gp        = gpar(fontsize = 18)
     )
# dev.off()
##'-----------------------------------------------------------------------------------------#





