##Figure3I
library(Seurat)

##import the singlecell data
single_data<-readRDS("C:/Users/Lenovo/Desktop/fetal_cmz_60_90_harmnoy_30pcs.rds")
DimPlot(single_data, reduction = "umap", group.by = "final_annotated_0807_2", label = TRUE, repel = TRUE)
single_data$final_annotated_0807_2<-factor(single_data$final_annotated_0807_2,levels = c("RPE","RPE_progenitors","hRSLCs","RPCs","PC_precursors",
                                                                                         "PCs","RGCs","ACs","HCs","BCs"))

Idents(single_data)<-single_data$tissue
single_data_fetal<-subset(single_data,idents = "Fetal")
single_data_Organoid<-subset(single_data,idents = "Organoid")

##fetal
colors<-c("#71C89C","#67ADB7","#36600E", "#6A8473", "#C0BFDF","#E77A77","#7B6148","#ECBA84","#CA8C74","#A13B46","#9569AB")
DimPlot(single_data_fetal, reduction = "umap", group.by = "final_annotated_0807_2", cols = colors)

##Organoid
DimPlot(single_data_Organoid, reduction = "umap", group.by = "final_annotated_0807_2", cols = colors)
