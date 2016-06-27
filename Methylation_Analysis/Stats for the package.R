# #-----------------------------------------------------------------------------------------------------------------------------#
#  Author      : Ruddy                                                                                                        |
#  Language    : R Statistical Programming Language                                                                           |
#  Optimised   :                                                                                                              |
#  Data Owner  : Ruddy                                                                                                        |
#  Description : Statistical calculations for the TXF pipeline (fisher test)                                                  |
#-----------------------------------------------------------------------------------------------------------------------------#

#from the output of the Allfatesandtissues_TXF_hypohyper_ceros data.frame I have created another dataframe for the calculation of the fisher test.  I just modified the 
#columns with the percentajes and I have added the values of the number of CPgs not binding to the TXFs

Allfatesandtissues_TXF_hypohyper_ceros_fixed <- read.table("~/Documents/Newcastle Data/common files/Allfatesandtissues_TXF_hypohyper_ceros_fixed.txt", header=TRUE, quote="\"")
View(Allfatesandtissues_TXF_hypohyper_ceros_fixed)# this data.frame was generated after the aplication of the main code of the package for the different lists of cpgs
Allfates_statistics_fixed <- read.table("~/Documents/Newcastle Data/common files/Allfates_statistics_fixed.txt", header=TRUE, quote="\"")
View(Allfates_statistics_fixed)# the data.frame with the not bindings
# the info about the list of cpgs and their column number in the other data.frames is in the next data.frame
df_TXF_proportions_stats <- read.table("~/Documents/Newcastle Data/common files/df_TXF_proportions_stat.txt", header=TRUE, quote="\"")
df_TXF_proportions_stats
#note that column v4 is the number of the columns for the not binding 
#note that column v5 is the number of the columnd for the binding

############################################### this function is for a given TXF in a constant list of CPGs in this case is the column 11 hypo_OB
fisher_stats<-function(x){
  numberofTXF<-match(x,Allfates_statistics_fixed[,1])
  Convictions <-matrix(c(Allfates_statistics_fixed[numberofTXF,11],Allfates_statistics_fixed[numberofTXF,13] , Allfates_statistics_fixed[numberofTXF,15], Allfates_statistics_fixed[numberofTXF,2]),
                       nrow = 2,
                       dimnames =
                         list(c("TXF-Bind", "TXF-Notbind"),
                              c("Observed", "array")))
  Convictions_stats<-fisher.test(Convictions)
  return( Convictions_stats$p.value)
}

fisher_stats("ZEB1")
# to validate the output I compare the columns and check that the output of the fisher_stats("ZEB1") is correct
TXF2plot<-"ZEB1"#viejo sin normalizar
numberofTXF<-match(TXF2plot,Allfatesandtissues_TXF_hypohyper_ceros[,1])
numberofTXF
barplot(as.numeric(Allfatesandtissues_TXF_hypohyper_ceros[numberofTXF,c(2,13,14,18,19,22,23,26,27,34,35,38,39,9,10, 8,7,46,47,42,43, 50,51)]),
        names.arg=c("Array","OB","OB","Adipose","Adipose","Muscle","Muscle","hPSC","hPSC","Pancreas","Pancreas","Thymus", "Thymus","OA", "OA", "CH","CH","CH5","CH5","hipClust","hipClust", "kneeClust", "kneeClust"),
        col=c("red", "white","black","white","black","white","black","white","black","white","black","white","black","white","blue","white","green","white","purple","white","chocolate1","white","bisque"), 
        main=TXF2plot, ylab="Percentage %",las=2 )


# this function is for determine a certain p-value for a given TXF and a given set of CpGs associated with a process
fisher_B<-function(x,y){
  numberofTXF<-match(x,Allfates_statistics_fixed[,1])
  numberoCpGlist<-match(y,df_TXF_proportions_stats[,1])
  Convictions <-matrix(c(Allfates_statistics_fixed[numberofTXF,df_TXF_proportions_stats[numberoCpGlist,5]],Allfates_statistics_fixed[numberofTXF,df_TXF_proportions_stats[numberoCpGlist,4]] , Allfates_statistics_fixed[numberofTXF,15], Allfates_statistics_fixed[numberofTXF,2]),
                       nrow = 2,
                       dimnames =
                         list(c("TXF-Bind", "TXF-Notbind"),
                              c("Observed", "array")))
  Convictions_stats<-fisher.test(Convictions)
  return( Convictions_stats$p.value)
}

fisher_B("ZEB1","Thym_hypo" )# this function needs two arguments
############################################################################################################ For all the TXFs for a give list

fisher4all_lists<-function(Loo){lapply(Allfates_statistics_fixed[,1],fisher_B,y=Loo)}
fisher4all_lists("OB_Hypo")#all the p-values of the fisher test for a certain list of CPGs

#for all the TXFs and for all the list of Cpgs at the same time
Tabla_fisher_test<-sapply(df_TXF_proportions_stats[,1],fisher4all_lists)
rownames(Tabla_fisher_test)<-Allfates_statistics_fixed[,1]
colnames(Tabla_fisher_test)<-df_TXF_proportions_stats[,1]
Tabla_fisher_test

write.table(Tabla_fisher_test, file = "~/Documents/Newcastle Data/common files/Fisher_of_all_the_sets.txt", append = FALSE, quote = TRUE, sep = " ",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE,
            col.names = TRUE, qmethod = c("escape", "double"))
