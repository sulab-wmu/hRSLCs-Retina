##Supplemental_Figure S3C1-C4
library(BiocManager)
library(Seurat)
library(monocle)
library(tidyverse)
library(patchwork)
library(sctransform)
library(S4Vectors)
library(Signac)
library(dplyr)
library(ggplot2)

##import the result from fetal_monocle2.R
cds<-readRDS("C:/Users/Lenovo/Desktop/final_version/monocle2/cds.rds")

##import the singlecell data from fetal_celltypes_annotation.R
data_fetal_scRNA<-readRDS("C:/Users/Lenovo/Desktop/final_version/annotation/final_annotation.rds")

##plot by celltypes(C1)
pdf("monocle cell_type3.pdf",height = 7,width=7)
plot_cell_trajectory(cds,color_by="wsnn_res_0.5_cell_type", size=1,show_backbone=TRUE)+
  scale_colour_manual(values =c("#ECBA84","#A13B46","#CA8C74","#36600E","#9569AB","#C0BFDF","#E77A77","#7B6148","#6A8473","#71C89C","#67ADB7"))
dev.off()

##facet_wrap(C3,C4)
pdf("monocle cell_type.faceted3.pdf",height = 7,width=7)
plot_cell_trajectory(cds, cell_size = 0.5, color_by = "wsnn_res_0.5_cell_type") + facet_wrap("~wsnn_res_0.5_cell_type", nrow = 4 )+
  scale_colour_manual(values = c("#ECBA84","#A13B46","#CA8C74","#36600E","#9569AB","#C0BFDF","#E77A77","#7B6148","#6A8473","#71C89C","#67ADB7"))
dev.off()


##choose hRSLCs as root
GM_state <- function(cds){
  if (length(unique(pData(cds)$State)) > 1){
    T0_counts <- table(pData(cds)$State, pData(cds)$wsnn_res_0.5_cell_type)[,"hRSLCs"]
    return(as.numeric(names(T0_counts)[which
                                       (T0_counts == max(T0_counts))]))
  } else {
    return (1)
  }
}
cds1 <- orderCells(cds, root_state = GM_state(cds))

##plot by Pseudotime(C2)
pdf("monocle pseudotime3.pdf",height = 7,width = 7)
plot_cell_trajectory(cds1,color_by="Pseudotime", size=1,show_backbone=TRUE)
dev.off()






