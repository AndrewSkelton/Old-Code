plot(analysis_ready_data[, c(102,117,118,122,123,135,138,139,140,141,142,144,159)], 
     cex=0.5, what="MAplot", main=paste(sampleNames(analysis_ready_data)[31:36], collapse=", "))

outliers_intensity <- sampleNames(analysis_ready_data)[c(102,117,118,122,123,135,138,139,140,141,142,144,159)]
paste(outliers_intensity, collapse=", ")