##run CytoTRACE
library(Seurat)
library(CytoTRACE)

setwd("/share2/pub/chenchg/chenchg/SingleCell/LiHui/final/CytoTRACE/")
data<-readRDS("/share2/pub/chenchg/chenchg/SingleCell/LiHui/final/final_annotation.rds")

results <- CytoTRACE(as.matrix(data@assays$RNA@counts), enableFast=F)
Idents(data)<-data$wsnn_res_0.5_cell_type
celltype = as.character(data$wsnn_res_0.5_cell_type)
names(celltype) = row.names(data@meta.data)

##save the result
saveRDS(results,file="CytoTRACE.rds")
