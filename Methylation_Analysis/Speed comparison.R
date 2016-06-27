#-----------------------------------------------------------------------------------------------------------------------------#
#  Author      : Ruddy                                                                                                        |
#  Language    : R Statistical Programming Language                                                                           |
#  Optimised   :                                                                                                              |
#  Data Owner  : Ruddy                                                                                                        |
#  Description : this script is intended to improve the speed of the funtions TFnames4cpg                                     |
#-----------------------------------------------------------------------------------------------------------------------------#

CH_5_solo_tissue_hypocpgs<-c("cg14193430", "cg24116380", "cg27645858", "cg09599062")
#Initially I had this funtion
TFnames4cpg<-function (x){
  B<-InfiniumMethylation[ which(names(InfiniumMethylation)==x) ]
  A<-TXF_Ruddy
  D<-as.data.frame(subsetByOverlaps(A,B))
  df<-as.data.frame (rownames(D))
  colnames(df)<-c("TXF")
  df<-t(df)
  
  return(df) 
}

#then I changed the sentence : B<-InfiniumMethylation[ which(names(InfiniumMethylation)==x) ] to B<-InfiniumMethylation[ x ] and work the same

TFnames4cpg<-function (x){
  B<-InfiniumMethylation[ x ]
  A<-TXF_Ruddy_2014
  D<-as.data.frame(subsetByOverlaps(A,B))
  df<-as.data.frame (rownames(D))
  colnames(df)<-c("TXF")
  df<-t(df)
  
  return(df) 
}

# in order to increase the speed I have also deleted the calls to as.data.frame function. To do this I used the function names of the resulting genomic range
#also I put all the funtion in the same line

TFnames4cpg_prueba<-function (x){names(subsetByOverlaps(TXF_Ruddy_2014,InfiniumMethylation[ x ]))}

#Test1
TFnames4cpg("cg24116380")
TFnames4cpg_prueba("cg24116380")
#test2
lapply(CH_5_solo_tissue_hypocpgs,TFnames4cpg)
lapply(CH_5_solo_tissue_hypocpgs,TFnames4cpg_prueba)
#test3
unlist(lapply(CH_5_solo_tissue_hypocpgs,TFnames4cpg))
unlist(lapply(CH_5_solo_tissue_hypocpgs,TFnames4cpg_prueba))
#test4
table(unlist(lapply(CH_5_solo_tissue_hypocpgs,TFnames4cpg)))
table(unlist(lapply(CH_5_solo_tissue_hypocpgs,TFnames4cpg_prueba)))

# I don't appreciate any improvement....I should use any funtion to estimate the time used by the system

#However after apply the new modifications in the function it works amazingly FASTSSSSSSS

TFnames4cpg_v2<-function (x){
  cov_matrix<-as.matrix( countOverlaps(TXF_Ruddy_2014,InfiniumMethylation[ x ]))
  split_cov_matrix<-split(cov_matrix, rownames(cov_matrix))
  return(as.matrix(lapply(split_cov_matrix,sum)))
}

TFnames4cpg_v2("cg24116380")
TFnames4cpg_v2(OB_array2_cpglist_hypo[,1])
TFnames4cpg_v2(OB_array2_cpglist_hyper[,1])

#this means that the function works with only one value or with a vector, which is quite amazing

