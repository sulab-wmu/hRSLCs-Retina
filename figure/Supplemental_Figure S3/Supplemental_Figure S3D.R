##Supplemental_Figure S3D
library(cicero)
library(monocle3)
library(Seurat)
library(tidyverse)
library(patchwork)
library(SingleCellExperiment)
set.seed(1234)

setwd("C:/Users/Lenovo/Desktop/final_version/monocle3_scATAC/")

##import the result from fetal_monocle3_scATAC.R
cds<-readRDS("scATAC_cds.rds")

##import the singlecell data from fetal_celltypes_annotation.R
data_fetal_scRNA<-readRDS("C:/Users/Lenovo/Desktop/final_version/annotation/final_annotation.rds")

##choose hRSLCs as root
##D1
library(ggplot2)
cds <- order_cells(cds)
p5<-plot_cells(cds,color_cells_by="pseudotime",label_groups_by_cluster=FALSE,
               label_leaves=FALSE,label_branch_points=FALSE)+
  scale_color_gradient(low="#B0CFE4", high = "#B22028")
ggsave("monocle3_trajectory_pseudotime.pdf",p5,width = 6,height = 5)

##pseudotime
pseudotime_atac<-pseudotime(cds, reduction_method ="UMAP")

data_fetal_scRNA$pseudotime_atac <- pseudotime_atac

head(data_fetal_scRNA$pseudotime_atac)

pseudotime_atac <- as.numeric(data_fetal_scRNA$pseudotime_atac)
celltypes <- as.character(data_fetal_scRNA$wsnn_res_0.5_cell_type)
data_pseudotime_atac <- as.data.frame(cbind(celltypes,pseudotime_atac))
names(data_pseudotime_atac) <- c("celltypes_0.05","pseudotime_atac")

### subset of four cell types in Lineage 1
data_pseudotime_subset_atac <- data_pseudotime_atac[which(data_pseudotime_atac$celltypes %in% c("hRSLCs","PC_precursors","PCs","RPCs")),]
head(data_pseudotime_subset_atac)

### subset of two cell types in Lineage 2
data_pseudotime_subset2_atac <- data_pseudotime_atac[which(data_pseudotime_atac$celltypes %in% c("RPE","RPE_progenitors")),]
head(data_pseudotime_subset2_atac)


data_pseudotime_subset_atac$celltypes_0.05<-factor(data_pseudotime_subset_atac$celltypes_0.05,levels = c("hRSLCs","RPCs","PC_precursors","PCs"))

library(ggplot2)
##D2
p1 = ggplot(data_pseudotime_subset_atac, aes(x=as.numeric(pseudotime_atac),fill=celltypes_0.05)) + 
          theme_classic()+
          geom_density()+
          scale_fill_manual(values = c("#af2157","#BD8098","#C0BFDF","#67ADB7"))+
          theme(axis.title = element_text(size=18),
          axis.text = element_text(size=14),
          axis.text.x = element_text(angle = 35,
                                     hjust=1,
                                     size = 10,
                                     color = "black"))
p1

##D3                                 
data_pseudotime_subset2_atac$celltypes_0.05<-factor(data_pseudotime_subset2_atac$celltypes_0.05,levels = c("RPE","RPE_progenitors"))
p2 = ggplot(data_pseudotime_subset2_atac, aes(x=as.numeric(pseudotime_atac),fill=celltypes_0.05)) + 
           theme_classic()+
           geom_density()+
           scale_fill_manual(values = c("#67ADB7","#BD8098"))+
           scale_y_sqrt(limits = c(0,100),breaks = c(0,3,5,10,20,30,50,100))+
           theme(axis.title = element_text(size=18),
           axis.text = element_text(size=14),
           axis.text.x = element_text(angle = 35,
                                      hjust=1,
                                      size = 10,
                                      color = "black"))
p2 +xlim(0,15)















