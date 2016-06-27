#-------------------------------------------------------------------------------------------#
#  Author      : Andrew J Skelton                                                           |
#  Language    : R Statistical Programming Language                                         |
#  Study       : Utility Code                                                               |
#  Data Owner  : NA                                                                         |
#  Description : Create a .xlsx file, add sheets, data frames, and images                   |
#-------------------------------------------------------------------------------------------#



##'Set the Working Directory, load essential packages, and import files
##'-----------------------------------------------------------------------------------------#
library(xlsx)
##'-----------------------------------------------------------------------------------------#


##'Create a WorkBook data structure
##'-----------------------------------------------------------------------------------------#
wb               <- createWorkbook()
##'-----------------------------------------------------------------------------------------#


##'Create a Sheet and add it to the workbook - Then add a dataframe to the sheet
##'-----------------------------------------------------------------------------------------#
sheet            <- createSheet(wb, sheetName="RA_Expression")
addDataFrame(signature_exp[,grep("^RA$", pData(raw_data_det)$Update_Collapse_2)],
             sheet, startRow=1, startColumn=1)
##'-----------------------------------------------------------------------------------------#

##'Add the rest of the sheets
##'-----------------------------------------------------------------------------------------#
sheet            <- createSheet(wb, sheetName="RA_Pheno")
addDataFrame(pData(raw_data_det)[grep("^RA$", pData(raw_data_det)$Update_Collapse_2),],
             sheet, startRow=1, startColumn=1)

sheet            <- createSheet(wb, sheetName="NRA_Expression")
addDataFrame(signature_exp[,grep("^NRA$", pData(raw_data_det)$Update_Collapse_2)],
             sheet, startRow=1, startColumn=1)
sheet            <- createSheet(wb, sheetName="NRA_Pheno")
addDataFrame(pData(raw_data_det)[grep("^NRA$", pData(raw_data_det)$Update_Collapse_2),],
             sheet, startRow=1, startColumn=1)

sheet            <- createSheet(wb, sheetName="UA_Expression")
addDataFrame(signature_exp[,grep("^UA$", pData(raw_data_det)$Update_Collapse_2)],
             sheet, startRow=1, startColumn=1)
sheet            <- createSheet(wb, sheetName="UA_Pheno")
addDataFrame(pData(raw_data_det)[grep("^UA$", pData(raw_data_det)$Update_Collapse_2),],
             sheet, startRow=1, startColumn=1)
##'-----------------------------------------------------------------------------------------#


##'Create a Sheet and add an image to it
##'-----------------------------------------------------------------------------------------#
sheet            <- createSheet(wb, sheetName="Heatmap_WithUA")
addPicture("heatmap_all.png", sheet, scale = 1, startRow = 1, startColumn = 1)
##'-----------------------------------------------------------------------------------------#

##'Add the rest of the images...
##'-----------------------------------------------------------------------------------------#
sheet            <- createSheet(wb, sheetName="Heatmap_WithoutUA")
addPicture("heatmap_RA_NRA.png", sheet, scale = 1, startRow = 1, startColumn = 1)

sheet            <- createSheet(wb, sheetName="Probe_Profile")
addPicture("probe_profile.png", sheet, scale = 1, startRow = 1, startColumn = 1)

sheet            <- createSheet(wb, sheetName="Dendrogram_WithUA")
addPicture("dendrogram.png", sheet, scale = 1, startRow = 1, startColumn = 1)

sheet            <- createSheet(wb, sheetName="Dendrogram_WithoutUA")
addPicture("dendrogram_UA.png", sheet, scale = 1, startRow = 1, startColumn = 1)
##'-----------------------------------------------------------------------------------------#


##'Save the workbook and delete temporary images
##'-----------------------------------------------------------------------------------------#
saveWorkbook(wb, file="Array_BC_ComboV2_ComBat_Subset.xlsx")

file.remove(c("heatmap_all.png", "heatmap_RA_NRA.png", "probe_profile.png",
              "dendrogram_UA.png", "dendrogram.png"))
##'-----------------------------------------------------------------------------------------#
