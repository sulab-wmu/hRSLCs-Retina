##CytoTRACE
##organoids
library(Seurat)
library(CytoTRACE)

setwd("/share2/pub/chenchg/chenchg/SingleCell/LiHui/regeneration/CytoTRACE/")
##organoids singlecell data 
data<-readRDS("/share2/pub/chenchg/chenchg/SingleCell/LiHui/regeneration/Final_Reg_CMZ_0807.rds")

results <- CytoTRACE(as.matrix(data@assays$RNA@counts), enableFast=F)
Idents(data)<-data$final_annotated_0807_2
celltype = as.character(data$final_annotated_0807_2)
names(celltype) = row.names(data@meta.data)

##save the result
saveRDS(results,file="CytoTRACE.rds")
