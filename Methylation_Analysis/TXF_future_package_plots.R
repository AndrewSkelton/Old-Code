#-----------------------------------------------------------------------------------------------------------------------------#
#  Author      : Ruddy                                                                                                        |
#  Language    : R Statistical Programming Language                                                                           |
#  Optimised   :                                                                                                              |
#  Data Owner  : Ruddy                                                                                                        |
#  Description : Plots designed to use associated with the package of the TXF_pipeline                                        |
#-----------------------------------------------------------------------------------------------------------------------------#


#This will plot for all the datasets in the outputdataframe
TXF2plot<-"ZEB1"#without further normalization of the number of cpgs outside of the TXFs
numberofTXF<-match(TXF2plot,Allfatesandtissues_TXF_hypohyper_ceros[,1])
numberofTXF
barplot(as.numeric(Allfatesandtissues_TXF_hypohyper_ceros[numberofTXF,c(7,5,6)]),
        names.arg=c("Array","OA_cpglist_hypo","OA_cpglist_hyper"),
        col=c("red", "white","black"), 
        main=TXF2plot, ylab="Percentage %",las=2 )


