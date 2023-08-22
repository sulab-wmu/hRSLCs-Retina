##Figure B-E plot+Supplemental_Figure S2A
library(Seurat)
library(Signac)
library(ggplot2)
library(openxlsx)

##import the scRNA data from fetal_celltypes_annotation.R
single_data<-readRDS("C:/Users/Lenovo/Desktop/final_version/annotation/final_annotation.rds")

##Figure 1B
DimPlot(single_data,group.by = "orig.ident",reduction = "umap.biMod")

##Figure 1C
colors<-c("#71C89C","#67ADB7","#36600E", "#6A8473", "#C0BFDF","#E77A77","#7B6148","#ECBA84","#CA8C74","#A13B46","#9569AB")
DimPlot(single_data,group.by = "wsnn_res_0.5_cell_type",reduction = "umap.biMod",cols = colors)

##input the markers for each cell type
marker<-read.xlsx("C:/Users/Lenovo/Desktop/github_code/analysis/figure1/fetal Marker Genes.xlsx")

##Figure 1D
DotPlot(single_data,features=rev(unique(marker$GENES)),group.by="wsnn_res_0.5_cell_type",assay="SCT")+
  scale_colour_gradient2(low="#3A71AA",mid="white",high="#B22028",midpoint=0.5)+
  theme_bw()+
  coord_flip()+
  theme (axis.text.x = element_text (angle = 35, hjust = 1))

##figure 1E
##density_plot_for individual genes
#BiocManager::install("Nebulosa")change color of plot_density
#install.packages("ks")
#devtools::install_github(repo = "samuel-marsh/scCustomize")
library(Nebulosa)
library(BiocFileCache)
library(paletteer)
library(scCustomize)

##Stem1
Plot_Density_Custom(seurat_object =single_data, features = "MECOM",reduction = "umap.biMod",
                    custom_palette = c("#B0CFE4","#FACABC","#E77A77","#DC0000FF"))
Plot_Density_Custom(seurat_object =single_data, features = "CPAMD8",reduction = "umap.biMod",
                    custom_palette = c("#B0CFE4","#FACABC","#E77A77","#DC0000FF"))

#Stem2
Plot_Density_Custom(seurat_object =single_data, features = "MITF",reduction = "umap.biMod",
                    custom_palette = c("#B0CFE4","#FACABC","#E77A77","#DC0000FF"))
Plot_Density_Custom(seurat_object =single_data, features = "GJA1",reduction = "umap.biMod",
                    custom_palette = c("#B0CFE4","#FACABC","#E77A77","#DC0000FF"))


##Supplemental_Figure S2A
##RPE
Plot_Density_Custom(seurat_object =single_data, features = "BEST1",reduction = "umap.biMod",
                    custom_palette = c("#B0CFE4","#FACABC","#E77A77","#DC0000FF"))
Plot_Density_Custom(seurat_object =single_data, features = "CD96",reduction = "umap.biMod",
                    custom_palette = c("#B0CFE4","#FACABC","#E77A77","#DC0000FF"))

##RPC
Plot_Density_Custom(seurat_object =single_data, features = "CCND1",reduction = "umap.biMod",
                    custom_palette = c("#B0CFE4","#FACABC","#E77A77","#DC0000FF"))
Plot_Density_Custom(seurat_object =single_data, features = "SPP1",reduction = "umap.biMod",
                    custom_palette = c("#B0CFE4","#FACABC","#E77A77","#DC0000FF"))

##PC_precursors
Plot_Density_Custom(seurat_object =single_data, features = "ATOH7",reduction = "umap.biMod",
                    custom_palette = c("#B0CFE4","#FACABC","#E77A77","#DC0000FF"))
Plot_Density_Custom(seurat_object =single_data, features = "NFIB",reduction = "umap.biMod",
                    custom_palette = c("#B0CFE4","#FACABC","#E77A77","#DC0000FF"))

##PCs
Plot_Density_Custom(seurat_object =single_data, features = "RCVRN",reduction = "umap.biMod",
                    custom_palette = c("#B0CFE4","#FACABC","#E77A77","#DC0000FF"))
Plot_Density_Custom(seurat_object =single_data, features = "RP1",reduction = "umap.biMod",
                    custom_palette = c("#B0CFE4","#FACABC","#E77A77","#DC0000FF"))

##RGCs
Plot_Density_Custom(seurat_object =single_data, features = "TUBB2A",reduction = "umap.biMod",
                    custom_palette = c("#B0CFE4","#FACABC","#E77A77","#DC0000FF"))
Plot_Density_Custom(seurat_object =single_data, features = "EBF1",reduction = "umap.biMod",
                    custom_palette = c("#B0CFE4","#FACABC","#E77A77","#DC0000FF"))

##ACs
Plot_Density_Custom(seurat_object =single_data, features = "MEIS2",reduction = "umap.biMod",
                    custom_palette = c("#B0CFE4","#FACABC","#E77A77","#DC0000FF"))
Plot_Density_Custom(seurat_object =single_data, features = "MIR181A1HG",reduction = "umap.biMod",
                    custom_palette = c("#B0CFE4","#FACABC","#E77A77","#DC0000FF"))

##HCs
Plot_Density_Custom(seurat_object =single_data, features = "PROX1",reduction = "umap.biMod",
                    custom_palette = c("#B0CFE4","#FACABC","#E77A77","#DC0000FF"))
Plot_Density_Custom(seurat_object =single_data, features = "TMOD1",reduction = "umap.biMod",
                    custom_palette = c("#B0CFE4","#FACABC","#E77A77","#DC0000FF"))

##BCs
Plot_Density_Custom(seurat_object =single_data, features = "CA10",reduction = "umap.biMod",
                    custom_palette = c("#B0CFE4","#FACABC","#E77A77","#DC0000FF"))
Plot_Density_Custom(seurat_object =single_data, features = "VSX1",reduction = "umap.biMod",
                    custom_palette = c("#B0CFE4","#FACABC","#E77A77","#DC0000FF"))

##MCs
Plot_Density_Custom(seurat_object =single_data, features = "CLU",reduction = "umap.biMod",
                    custom_palette = c("#B0CFE4","#FACABC","#E77A77","#DC0000FF"))
Plot_Density_Custom(seurat_object =single_data, features = "GPX3",reduction = "umap.biMod",
                    custom_palette = c("#B0CFE4","#FACABC","#E77A77","#DC0000FF"))
































