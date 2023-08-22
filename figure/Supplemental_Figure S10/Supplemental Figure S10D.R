##Supplemental Figure S10D
library(BiocManager)
library(Seurat)
library(monocle)
library(tidyverse)
library(patchwork)
library(sctransform)
library(S4Vectors)
library(Signac)
library(dplyr)

##organoids singlecell data
cds<-readRDS("organoids_monocle2_cds.rds")

##Supplemental Figure S10D1
pdf("monocle cell_type2.pdf",height = 7,width=7)
plot_cell_trajectory(cds,color_by="final_annotated_0807_2", size=1,show_backbone=TRUE)+
  scale_colour_manual(values = c("#ECBA84", "#A13B46", "#CA8C74", "#36600E","#C0BFDF","#E77A77","#7E6148FF","#6A8473","#85C89C","#67ADB7"))
dev.off()

##Supplemental Figure S10D3
pdf("monocle cell_type.faceted2.pdf",height = 7,width=7)
plot_cell_trajectory(cds, cell_size = 0.5, color_by = "final_annotated_0807_2") + facet_wrap("~final_annotated_0807_2", nrow = 4 )+
  scale_colour_manual(values = c("#ECBA84", "#A13B46", "#CA8C74", "#36600E","#C0BFDF","#E77A77","#7E6148FF","#6A8473","#85C89C","#67ADB7"))
dev.off()

##choose hRSLCs as root
GM_state <- function(cds){
  if (length(unique(pData(cds)$State)) > 1){
    T0_counts <- table(pData(cds)$State, pData(cds)$final_annotated_0807_2)[,"hRSLCs"]
    return(as.numeric(names(T0_counts)[which
                                       (T0_counts == max(T0_counts))]))
  } else {
    return (1)
  }
}
cds1 <- orderCells(cds, root_state = GM_state(cds))

##Supplemental Figure S10D2
pdf("monocle pseudotime2.pdf",height = 7,width = 7)
plot_cell_trajectory(cds1,color_by="Pseudotime", size=1,show_backbone=TRUE)
dev.off()