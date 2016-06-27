#-------------------------------------------------------------------------------------------#
#  Author      : Andrew J Skelton                                                           |
#  Language    : R Statistical Programming Language                                         |
#  Study       : Methylation of Cartilage Study - Hip / Knee OA                             |
#  Data Owner  : Newcastle University - Prof. John Loughlin, Dr. Louise Reynard             |
#  Description : Using research developed by Steve Horvath at UCLA - Attempt to             |
#                differentiate between methylation age and true age.                        |
#-------------------------------------------------------------------------------------------#



##'Set the Working Directory, load essential packages, and import files
##'-----------------------------------------------------------------------------------------#
setwd("~/Loughlin/project_master_branch/methylation/")

source("http://bioconductor.org/biocLite.R")
biocLite(c("RPMM"))

library(WGCNA)
library(sqldf)
library(impute)
library(BMIQ)
library(dynamicTreeCut)


##'-----------------------------------------------------------------------------------------#



##'Horvath Functions and Prep
##'-----------------------------------------------------------------------------------------#
source("Horvath/NORMALIZATION.R")
probeAnnotation21kdatMethUsed <- read.csv("Horvath/probeAnnotation21kdatMethUsed.csv")
probeAnnotation27k            <- read.csv("Horvath/datMiniAnnotation.csv")
datClock                      <- read.csv("Horvath/AdditionalFile3.csv")

trafo      <- function(x, adult.age=20) { 
  x <- (x+1)/(1+adult.age)
  y <- ifelse(x<=1, log(x),x-1)
  y 
}

anti.trafo <- function(x, adult.age=20) { 
  ifelse(x<0, 
         (1+adult.age)*exp(x)-1, 
         (1+adult.age)*x+adult.age) 
}

asnumeric1 <- function(x) {
  as.numeric(as.character(x))
}

dat0              <- m2beta(exprs(lumi.norm[,c(1:5)]))
dat0              <- cbind(rownames(dat0),dat0)
colnames(dat0)[1] <- "CpG"
##'-----------------------------------------------------------------------------------------#


##'Horvath - Stage 1
##'-----------------------------------------------------------------------------------------#
nSamples     <- dim(dat0)[[2]]-1
nProbes      <- dim(dat0)[[1]]
dat0[,1]     <- gsub(x=dat0 [,1],pattern="\"",replacement="")

file.remove("LogFile.txt")
file.create("LogFile.txt")
DoNotProceed <- F

match1       <- na.omit(match(probeAnnotation21kdatMethUsed$Name , dat0[,1]))
dat1         <- dat0[match1,]
dat1         <- as.matrix(dat1)
cpgs         <- dat1[,1]
dat1         <- apply(as.matrix(dat1[,-1]), 2, asnumeric1)
rownames(dat1) <- cpgs

set.seed(1)
# Do you want to normalize the data (recommended)?
normalizeData <- T
# source("Horvath/StepwiseAnalysis.txt")

##'-----------------------------------------------------------------------------------------#


##'Horvath - Stage 2
##'-----------------------------------------------------------------------------------------#
meanMethBySample      <- as.numeric(apply(as.matrix(dat1), 2, mean, na.rm=T))
minMethBySample       <- as.numeric(apply(as.matrix(dat1), 2, min,  na.rm=T))
maxMethBySample       <- as.numeric(apply(as.matrix(dat1), 2, max,  na.rm=T))

datMethUsed           <- t(dat1)
# datMethUsed           <- apply(datMethUsed, 2, asnumeric1)
# colnames(datMethUsed) <- as.character(dat1[,1])
noMissingPerSample    <- apply(as.matrix(is.na(datMethUsed)), 1, sum)


gs_in <- probeAnnotation21kdatMethUsed[match(rownames(lumi.norm)[match1], 
                                             probeAnnotation21kdatMethUsed$Name),]$goldstandard2
if (normalizeData){
  datMethUsedNormalized <- BMIQcalibration(datM=datMethUsed, 
                                           goldstandard.beta=gs_in, 
                                           plots=FALSE)
}
##'-----------------------------------------------------------------------------------------#



##'Horvath - Stage 3
##'-----------------------------------------------------------------------------------------#
# datClock2        <- datClock[-(as.vector(na.omit(match(datClock$CpGmarker, rownames(annotation_snpIll))))),]

datClock        <- datClock[!(datClock$CpGmarker %in% rownames(annotation_snpIll)),]

selectCpGsClock <- is.element(dimnames(datMethUsedNormalized)[[2]], 
                              as.character(datClock$CpGmarker[-1]))
datMethClock0   <- data.frame(datMethUsedNormalized[,selectCpGsClock])
datMethClock    <- data.frame(datMethClock0[as.character(datClock$CpGmarker[-1])])
# dim(datMethClock)
predictedAge    <- as.numeric(anti.trafo(datClock$CoefficientTraining[1] + 
                                         as.matrix(datMethClock) %*% 
                                         as.numeric(datClock$CoefficientTraining[-1])))

datMethUsedNormalized2 <- data.frame(rbind(datMethUsedNormalized,
                                           datMethUsedNormalized))

datMethClock0          <- data.frame(datMethUsedNormalized2[,selectCpGsClock])
datMethClock           <- data.frame(datMethClock0[as.character(datClock$CpGmarker[-1])])
# dim(datMethClock)
predictedAge           <- as.numeric(anti.trafo(datClock$CoefficientTraining[1] + 
                                                as.matrix(datMethClock) %*% 
                                                as.numeric(datClock$CoefficientTraining[-1])))
predictedAge           <- predictedAge[1]

Comment                      <- ifelse(predictedAge <0, 
                                       "Negative DNAm age.", 
                                       ifelse(predictedAge > 100, 
                                              "Old DNAm age.", 
                                              rep("",length(predictedAge))))
Comment[is.na(predictedAge)] <- "Age prediction was not possible. "
##'-----------------------------------------------------------------------------------------#



datout <- data.frame(SampleID=colnames(dat1), 
                     DNAmAge=predictedAge, 
                     Comment, 
                     noMissingPerSample, 
                     meanMethBySample, 
                     minMethBySample,
                     maxMethBySample)


cbind(datout, pData(lumi.norm)[match(datout$SampleID, rownames(pData(lumi.norm))),]$AGE)


















