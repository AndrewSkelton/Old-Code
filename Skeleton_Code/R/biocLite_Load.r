#=======================================================#
#	biocLite Complete Load								#
#  														#
#	Description : Load all bioconductor Libraries		#
#														#
#=======================================================#

source("http://bioconductor.org/biocLite.R")
biocLite()

biocLite("limma")
biocLite("sva")
biocLite("lumi")
biocLite("annotate")
biocLite("lumiHumanAll.db")
biocLite("arrayQualityMetrics")
biocLite("GOstats")
biocLite("RamiGO")
biocLite("pathview")

library(limma)
library(sva)
library(lumi)
library(annotate)
library(lumiHumanAll.db)
library(arrayQualityMetrics)
library(GOstats)
library(RamiGO)
library(pathview)

install.packages("gplots")
install.packages("ggplot2")
library(gplots)
library(ggplot2)