##monocle3
##organoids
library(BiocManager)
library(Seurat)
library(monocle3)
library(SingleCellExperiment)
library(tidyverse)
library(patchwork)
library(ggplot2)
set.seed(123)    

setwd("C:/Users/Lenovo/Desktop/regeneration/monocle3")
##organoids singlecell data
single_data<-readRDS("C:/Users/Lenovo/Desktop/regeneration/Final_Reg_CMZ_0807.rds")

##create cds object
data <- GetAssayData(single_data,assay='RNA',slot='counts')
cell_metadata <- single_data@meta.data
gene_annotation <- data.frame(gene_short_name=rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,cell_metadata=cell_metadata,gene_metadata=gene_annotation)

##preprocess_cds
cds <- preprocess_cds(cds, num_dim = 50)
cds <- reduce_dimension(cds,preprocess_method='PCA')
colData(cds)

##umap
p1=plot_cells(cds, reduction_method="UMAP",color_cells_by = "final_annotated_0807_2") + ggtitle('monocle3.umap')

##raw umap
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(single_data,reduction='umap')
cds@int_colData$reducedDims$UMAP <- int.embed[rownames(cds.embed),]

p2=plot_cells(cds,reduction_method='UMAP',color_cells_by = "final_annotated_0807_2")+ggtitle('raw.umap')
p2
p=p1|p2
ggsave("combine_monocle3_raw_umap.pdf",p,width = 9.5)

###cluster_cells
cds <- cluster_cells(cds)

##learn_graph
cds <- learn_graph(cds)

p3 = plot_cells(cds, color_cells_by = "final_annotated_0807_2", 
                label_groups_by_cluster=FALSE,label_leaves=FALSE, label_branch_points=FALSE)
ggsave("Monocle3_Trajectory_withlabel.pdf", p3,width = 6,height = 5)


p4 <-plot_cells(cds,color_cells_by="final_annotated_0807_2",label_groups_by_cluster=FALSE,
                label_leaves=TRUE,label_branch_points=FALSE)
ggsave("Monocle3_trajectory_leaves_withoutlabel.pdf",p4,width = 6,height = 5)

##save the result
saveRDS(cds,file = "organoids_monocle3_cds.rds")
