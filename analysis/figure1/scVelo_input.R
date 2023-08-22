##scVelo input
library(Seurat)
library(Signac)

##input the single_cell data form fetal_celltypes_annotation.R 
single_data<-readRDS("C:/Users/Lenovo/Desktop/final_version/annotation/final_annotation.rds")

##In the subsequent analysis, we define Stem1 as hRSLCs, Stem2 as RPE_progenitors
DefaultAssay(single_data)<-"RNA"
Idents(single_data)<-single_data$wsnn_res_0.5_cell_type
DimPlot(single_data,group.by = "wsnn_res_0.5_cell_type",reduction = "umap.biMod",label = T,label.size = 3)

##change the Seurat object to h5ad
##devtools::install_github("https://github.com/jimhester/scater.git")
##devtools::install_github('satijalab/seurat-data')
##devtools::install_github("https://github.com/mojaveazure/seurat-disk.git")
library(scater)
library(Seurat)
library(SeuratData)
library(patchwork)
library(SeuratDisk)
library(dplyr)

slot(single_data$SCT@SCTModel.list[[1]], 'median_umi') = median(single_data$SCT@SCTModel.list[[1]]@cell.attributes$umi)
SaveH5Seurat(single_data, filename = "data.h5Seurat")
Convert("data.h5Seurat", dest = "h5ad")

##save the metadata
metadata<-single_data@meta.data
write.csv(metadata,file = "metadata.csv",quote = F,col.names = T,row.names = T)




