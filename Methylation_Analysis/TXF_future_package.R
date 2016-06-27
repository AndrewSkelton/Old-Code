#-----------------------------------------------------------------------------------------------------------------------------#
#  Author      : Ruddy                                                                                                        |
#  Language    : R Statistical Programming Language                                                                           |
#  Optimised   :                                                                                                              |
#  Data Owner  : Ruddy                                                                                                        |
#  Description : A pipeline to study the impact of the epigenetic modifications in the TXF binding sites                      |
#-----------------------------------------------------------------------------------------------------------------------------#


####Preparing the data required for the pipeline#-----------------------------------------------------------------------------------------------------------------------------#
#libraries
library("rtracklayer", lib.loc="/Library/Frameworks/R.framework/Versions/3.1/Resources/library")
library("FDb.InfiniumMethylation.hg19", lib.loc="/Library/Frameworks/R.framework/Versions/3.1/Resources/library")
library("GenomicFeatures", lib.loc="/Library/Frameworks/R.framework/Versions/3.1/Resources/library")
library("GenomicRanges", lib.loc="/Library/Frameworks/R.framework/Versions/3.1/Resources/library")
library("lumi", lib.loc="/Library/Frameworks/R.framework/Versions/3.1/Resources/library")
##Dowloading to the work space the info required to make the ovelaping ( CpG genomic ranges and TXFs genomic ranges). 
#I have optimized everything for humans but it should not be difficult to do the same for other spcecies

#1) first the info of the cpgs of illumina array
InfiniumMethylation <- features(FDb.InfiniumMethylation.hg19)
## we’d prefer if R would stop us from comparing across assemblies:
met <- metadata(FDb.InfiniumMethylation.hg19) ## need to fetch genome
genome(InfiniumMethylation) <- met[which(met[,"name"]=="Genome"),"value"]
InfiniumMethylation <- sort(InfiniumMethylation)
show(InfiniumMethylation)

#2) the stuff for the TXFs

ucscGenomes()[ , "db"]
UCSCFeatureDbTableSchema(genome="hg19",
                         track="Txn Factor ChIP",
                         tablename="wgEncodeRegTfbsClusteredV3")
fdb<-makeFeatureDbFromUCSC(genome="hg19", ##### the download process could take a lot of time
                           track="Txn Factor ChIP",
                           tablename="wgEncodeRegTfbsClusteredV3")
TXF_Ruddy_2014<-features(fdb)
TXF_Ruddy_2014<-sort(TXF_Ruddy_2014)
show(TXF_Ruddy_2014)

#3) distribution of the bindings in the array

Array<-subsetByOverlaps(TXF_Ruddy_2014,InfiniumMethylation)
table(names(Array))
All_161_TXFs_old<-as.data.frame(table(names(Array)))
All_161_TXFs_old

fov<-findOverlaps(TXF_Ruddy_2014,InfiniumMethylation) 
fov
queryHits(fov)
cov<-(countOverlaps(TXF_Ruddy_2014,InfiniumMethylation))#I'm going to solve the problem with this
cov
# I have compared cov and fob and is the same the advantage is that cov have the mames of the TXFs
cov_matrix<-as.matrix(cov)
cov_matrix
split_cov_matrix<-split(cov_matrix, rownames(cov_matrix))
split_cov_matrix[1]
Number_binding_TXF<-as.matrix(lapply(split_cov_matrix,sum))
Number_binding_TXF# Yesss this is the solution
All_161_TXFs<-Number_binding_TXF
All_161_TXFs<-as.data.frame(All_161_TXFs)
All_161_TXFsv2<-All_161_TXFs_old
All_161_TXFsv2[,2]<-All_161_TXFs[1]
All_161_TXFs<-All_161_TXFsv2
write.table(as.matrix(All_161_TXFsv2), file = "~/Documents/Newcastle Data/common files/All_161_TXFs.txt", append = FALSE, quote = TRUE, sep = " ",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE,
            col.names = TRUE, qmethod = c("escape", "double"))

#4) Upload of the data to be analyzed______________________________________________________________________________________________________________#
# we will need a list of CpGs hypomethylated
CH_5_solo_tissue_hypocpgs#this is the name of the list
# we will need a list of CpGs hypermethylated
CH_5_solo_tissue_hypercpgs
#As an example we can use the dataframe that I created with some of our data
#example of cpgs for testing the function
CH_5_solo_tissue_hypocpgs<-c("cg14193430", "cg24116380", "cg27645858", "cg09599062")
CH_5_solo_tissue_hypercpgs<-c("cg14193430", "cg24116380", "cg27645858", "cg09599062")
#5) Core of the pipeline__________________________________________________________________________________________________________________________#

###calculation of the transcription factors associated with these cpgs


#optimized function....see mine_Speed comparison
TFnames4cpg_v2<-function (x){
  cov_matrix<-as.matrix( countOverlaps(TXF_Ruddy_2014,InfiniumMethylation[ x ]))
  split_cov_matrix<-split(cov_matrix, rownames(cov_matrix))
  return(as.matrix(lapply(split_cov_matrix,sum)))
}



#calculation
Frecuency_CH_5_solo_tissue_hypocpgs_TXF<-as.data.frame(TFnames4cpg_v2(CH_5_solo_tissue_hypocpgs))
Frecuency_CH_5_solo_tissue_hypocpgs_TXF

Frecuency_CH_5_solo_tissue_hypercpgs_TXF<-as.data.frame(TFnames4cpg_v2(CH_5_solo_tissue_hypercpgs))
Frecuency_CH_5_solo_tissue_hypercpgs_TXF


#6) storage and process of the results__________________________________________________________________________________________________________________________#

length(CH_5_solo_tissue_hypocpgs)#
length(CH_5_solo_tissue_hypercpgs)#
CH_5_solo_tissue_TXF_hypohyper<-merge(Frecuency_CH_5_solo_tissue_hypocpgs_TXF, Frecuency_CH_5_solo_tissue_hypercpgs_TXF, by.x="row.names",by.y="row.names",all.x=TRUE,all.y=TRUE)
CH_5_solo_tissue_TXF_hypohyper

#I use the work previously done with the Array distribution of the enrichment of the TXFs in the Array

Allfatesandtissues_TXF_hypohyper_ceros<-merge(as.data.frame(All_161_TXFs),CH_5_solo_tissue_TXF_hypohyper, by.x="Var1",by.y="Row.names",all.x=TRUE,all.y=TRUE)
Allfatesandtissues_TXF_hypohyper_ceros
#Now we calculate the percentages
Allfatesandtissues_TXF_hypohyper_ceros[,5]<-(as.numeric(Allfatesandtissues_TXF_hypohyper_ceros[,3])*100)/length(CH_5_solo_tissue_hypocpgs)
Allfatesandtissues_TXF_hypohyper_ceros[,6]<-(as.numeric(Allfatesandtissues_TXF_hypohyper_ceros[,4])*100)/length(CH_5_solo_tissue_hypercpgs)
Allfatesandtissues_TXF_hypohyper_ceros[,7]<-(as.numeric(Allfatesandtissues_TXF_hypohyper_ceros[,2])*100)/487173#this is the total number of CpGs in the array
Allfatesandtissues_TXF_hypohyper_ceros
#now I add the columns for the not bindig cpgs that we will use for the statistical calculations
Allfatesandtissues_TXF_hypohyper_ceros[,8]<-length(CH_5_solo_tissue_hypocpgs)-(as.numeric(Allfatesandtissues_TXF_hypohyper_ceros[,3]))
Allfatesandtissues_TXF_hypohyper_ceros[,9]<-length(CH_5_solo_tissue_hypercpgs)-(as.numeric(Allfatesandtissues_TXF_hypohyper_ceros[,4]))
Allfatesandtissues_TXF_hypohyper_ceros[,10]<-487173-(as.numeric(Allfatesandtissues_TXF_hypohyper_ceros[,2]))#this is the total number of CpGs in the array
Allfatesandtissues_TXF_hypohyper_ceros
#now we label the columns
colnames(Allfatesandtissues_TXF_hypohyper_ceros)<-c("TXF_Name", "Array_binding","CH_5_solo_tissue_hypocpgs","CH_5_solo_tissue_hypercpgs",
                                                    "%_CH_5_solo_tissue_hypocpgs","%_CH_5_solo_tissue_hypercpgs", "%_Array",
                                                    "not_bind_CH_5_solo_tissue_hypocpgs","not_bind_CH_5_solo_tissue_hypercpgs","not_bind_Array"
)

write.table(Allfatesandtissues_TXF_hypohyper_ceros, file = "~/Documents/Newcastle Data/common files/Frequency CH_OB_HIPOA/Allfatesandtissues_TXF_hypohyper_ceros_CH_5_solo.txt", append = FALSE, quote = TRUE, sep = " ",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE,
            col.names = TRUE, qmethod = c("escape", "double"))

#7) statisticts________________________________________________________________________________________________________________________________________________#
#this function is for a given TXF in a constant list of CPGs in this case is the column 3 
fisher_stats_v2<-function(x){
  numberofTXF<-match(x,Allfatesandtissues_TXF_hypohyper_ceros[,1])
  Convictions <-matrix(c(as.numeric(Allfatesandtissues_TXF_hypohyper_ceros[numberofTXF,3]),
                         as.numeric(Allfatesandtissues_TXF_hypohyper_ceros[numberofTXF,8]), 
                         as.numeric(Allfatesandtissues_TXF_hypohyper_ceros[numberofTXF,2]), 
                         as.numeric(Allfatesandtissues_TXF_hypohyper_ceros[numberofTXF,10])),
                       nrow = 2,
                       dimnames =
                         list(c("TXF-Bind", "TXF-Notbind"),
                              c("Observed", "array")))
  Convictions_stats<-fisher.test(Convictions)
  return( Convictions_stats$p.value)
}

fisher_stats_v2("STAT3")

# this function is for determine a certain p-value for a given TXF and a given set of CpGs associated with a process
fisher_B_v2<-function(x,y){
  numberofTXF<-match(x,Allfatesandtissues_TXF_hypohyper_ceros[,1])
  numberoCpGlist<-match(y,colnames(Allfatesandtissues_TXF_hypohyper_ceros))
  Convictions <-matrix(c(  as.numeric(Allfatesandtissues_TXF_hypohyper_ceros[numberofTXF, numberoCpGlist]),
                           as.numeric(Allfatesandtissues_TXF_hypohyper_ceros[numberofTXF, numberoCpGlist+5]), 
                           as.numeric(Allfatesandtissues_TXF_hypohyper_ceros[numberofTXF,2]), 
                           as.numeric(Allfatesandtissues_TXF_hypohyper_ceros[numberofTXF,10])),                       
                       nrow = 2,
                       dimnames =
                         list(c("TXF-Bind", "TXF-Notbind"),
                              c("Observed", "array")))
  Convictions_stats<-fisher.test(Convictions)
  return( Convictions_stats$p.value)
}

fisher_B_v2("STAT3","CH_5_solo_tissue_hypocpgs" )# this function needs two arguments

############################################################################################################ For all the TXFs for a give list

fisher4all_lists_v2<-function(Loo){lapply(Allfatesandtissues_TXF_hypohyper_ceros[,1],fisher_B_v2,y=Loo)}
fisher4all_lists_v2("CH_5_solo_tissue_hypocpgs")#all the p-values of the fisher test for a certain list of CPGs

#for all the TXFs and for all the list of Cpgs at the same time
Tabla_fisher_test<-sapply(colnames(Allfatesandtissues_TXF_hypohyper_ceros)[3:4],fisher4all_lists_v2)
rownames(Tabla_fisher_test)<-Allfatesandtissues_TXF_hypohyper_ceros[,1]
Tabla_fisher_test
write.table(Tabla_fisher_test, file = "~/Documents/Newcastle Data/common files/Fisher_of_all_the_sets.txt", append = FALSE, quote = TRUE, sep = " ",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE,
            col.names = TRUE, qmethod = c("escape", "double"))