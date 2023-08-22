##annotation
library(Seurat)
library(SeuratObject)
library(openxlsx)
setwd("C:/Users/Lenovo/Desktop/final_version/annotation")

##Input the single-cell data
single_data<-readRDS("C:/Users/Lenovo/Desktop/final_version/RDS_merged_reclustered-final.rds")

single_data <- FindClusters(object = single_data, resolution = 0.5, graph.name="wsnn")

##31cluster
DimPlot(single_data,group.by = "wsnn_res.0.5",reduction = "umap.biMod",label = T)

##input the marker gene for each cell type
marker<-read.xlsx("C:/Users/Lenovo/Desktop/github_code/figure1/analysis/fetal Marker Genes 整理20230802 调整GJA1顺序.xlsx")
DotPlot(single_data,features=unique(marker$GENES),group.by="wsnn_res.0.5",assay="SCT")+
  scale_colour_gradient2(low="#3A71AA",mid="white",high="#B22028",midpoint=0)+
  theme_bw()+
  coord_flip()+
  theme (axis.text.x = element_text (angle = 35, hjust = 1))


Idents(single_data)<-single_data$wsnn_res.0.5
wsnn_res.0.5_cell_type<-c("0"="PCs","1"="RPCs","2"="RPCs","3"="PCs","4"="ACs","5"="Stem1","6"="BCs","7"="RPCs","8"="PC_precursors",
                          "9"="Stem2","10"="BCs","11"="PCs","12"="ACs","13"="ACs","14"="ACs","15"="HCs","16"="RGCs","17"="ACs",
                          "18"="BCs","19"="RPCs","20"="BCs","21"="PC_precursors","22"="ACs","23"="Stem2","24"="RPE",
                          "25"="MCs","26"="HCs","27"="ACs","28"="RGCs","29"="BCs","30"="RGCs") 
single_data <- RenameIdents(single_data, wsnn_res.0.5_cell_type)
DimPlot(single_data, reduction = 'umap.biMod', label = TRUE, pt.size = 0.5) + NoLegend()

##add another column
single_data$wsnn_res_0.5_cell_type<-Idents(single_data)
DimPlot(single_data,group.by = "wsnn_res_0.5_cell_type",reduction = "umap.biMod",label = T)

##subset the cluster5 which is asscociated stemness
Idents(single_data)<-single_data$wsnn_res.0.5
data_subset<-subset(single_data,idents = "5")


##combine two kinds of omics to FindMultiModalNeighbors
data_subset <- FindMultiModalNeighbors(
  object = data_subset,
  reduction.list = list("pca", "lsi"), 
  dims.list = list(1:50, 2:50),
  modality.weight.name = "RNA.weight",
  verbose = TRUE
)
data_subset <- FindClusters(object = data_subset, resolution = 0.4, graph.name="wsnn")

##build a joint UMAP visualization
data_subset <- RunUMAP(
  object = data_subset,
  nn.name = "weighted.nn",
  #assay = "RNA",
  verbose = TRUE,
  reduction.name="umap.biMod",
  reduction.key = "UMAP_biMod_"
  #dims = 1:50
)

##marker dotplot
DimPlot(data_subset,group.by = "wsnn_res.0.4",reduction = "umap.biMod",label = T)+NoLegend()
DotPlot(data_subset,features=unique(marker$GENES),group.by="wsnn_res.0.4",assay="SCT")+
  scale_colour_gradient2(low="#3A71AA",mid="white",high="#B22028",midpoint=0)+
  theme_bw()+
  coord_flip()+
  theme (axis.text.x = element_text (angle = 35, hjust = 1))

##RenameIdents the cell
Idents(data_subset)<-data_subset$wsnn_res.0.4
single_data$wsnn_res_0.5_cell_type<-as.character(single_data$wsnn_res_0.5_cell_type)
single_data$wsnn_res_0.5_cell_type[WhichCells(data_subset,idents = c("0","2"))] <- "RPCs"
single_data$wsnn_res_0.5_cell_type[WhichCells(data_subset,idents = c("3"))] <- "PCs"
single_data$wsnn_res_0.5_cell_type[WhichCells(data_subset,idents = c("1"))] <- "Stem1"
table(single_data$wsnn_res_0.5_cell_type)
DimPlot(single_data,group.by = "wsnn_res_0.5_cell_type",reduction = "umap.biMod",label = T)+NoLegend()


##RenameIdents cluster6
DimPlot(single_data,group.by = "wsnn_res.0.5",reduction = "umap.biMod",label = T)+NoLegend()
Idents(single_data)<-single_data$wsnn_res.0.5
data_subset<-subset(single_data,idents = "6")

data_subset <- FindMultiModalNeighbors(
  object = data_subset,
  reduction.list = list("pca", "lsi"), 
  dims.list = list(1:50, 2:50),
  modality.weight.name = "RNA.weight",
  verbose = TRUE
)

data_subset <- FindClusters(object = data_subset, resolution = 4, graph.name="wsnn")

data_subset <- RunUMAP(
  object = data_subset,
  nn.name = "weighted.nn",
  #assay = "RNA",
  verbose = TRUE,
  reduction.name="umap.biMod",
  reduction.key = "UMAP_biMod_"
  #dims = 1:50
)

DimPlot(data_subset,group.by = "wsnn_res.4",reduction = "umap.biMod",label = T)+NoLegend()
DotPlot(data_subset,features=unique(marker$GENES),group.by="wsnn_res.4",assay="SCT")+
  scale_colour_gradient2(low="#3A71AA",mid="white",high="#B22028",midpoint=0)+
  theme_bw()+
  coord_flip()+
  theme (axis.text.x = element_text (angle = 35, hjust = 1))

Idents(data_subset)<-data_subset$wsnn_res.4
single_data$wsnn_res_0.5_cell_type[WhichCells(data_subset,idents = "6")] <- "PC_precursors"

table(single_data$wsnn_res_0.5_cell_type)
DimPlot(single_data,group.by = "wsnn_res_0.5_cell_type",reduction = "umap.biMod",label = T)+NoLegend()


##adjust the ordering of cell types
single_data$wsnn_res_0.5_cell_type <- factor(single_data$wsnn_res_0.5_cell_type,levels = c("RPE","Stem2","Stem1","RPCs","PC_precursors",
                                                                                           "PCs","RGCs","ACs","HCs","BCs","MCs"))

DotPlot(single_data,features=rev(unique(marker$GENES)),group.by="wsnn_res_0.5_cell_type",assay="SCT")+
  scale_colour_gradient2(low="#3A71AA",mid="white",high="#B22028",midpoint=0.5)+
  theme_bw()+
  coord_flip()+
  theme (axis.text.x = element_text (angle = 35, hjust = 1))

##save the single_cell data
saveRDS(single_data,file = "final_annotation.rds")








