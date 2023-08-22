##monocle3  Pseudotime trajectory analysis
#BiocManager::install("DelayedMatrixStats")
#BiocManager::install("SingleCellExperiment")
#BiocManager::install(c( 'Matrix.utils'))
#devtools::install_github('cole-trapnell-lab/leidenbase')
#devtools::install_github('cole-trapnell-lab/monocle3')
library(BiocManager)
library(Signac)
library(BiocGenerics)
library(monocle3)
library(SingleCellExperiment)
library(tidyverse)
library(patchwork)
library(ggplot2)
set.seed(123)    

##input the single_cell data from fetal_celltypes_annotation.R
single_data<-readRDS("C:/Users/Lenovo/Desktop/final_version/annotation/final_annotation.rds")

##create CDS object
data <- GetAssayData(single_data,assay='RNA',slot='counts')
cell_metadata <- single_data@meta.data
gene_annotation <- data.frame(gene_short_name=rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,cell_metadata=cell_metadata,gene_metadata=gene_annotation)

##data preprocess
cds <- preprocess_cds(cds, num_dim = 50)
cds <- reduce_dimension(cds,preprocess_method='PCA')
colData(cds)

##umap
p1=plot_cells(cds, reduction_method="UMAP",color_cells_by = "wsnn_res_0.5_cell_type") + ggtitle('monocle3.umap')

##input the raw umap to the cds
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(single_data,reduction='umap.biMod')
cds@int_colData$reducedDims$UMAP <- int.embed[rownames(cds.embed),]

p2=plot_cells(cds,reduction_method='UMAP',color_cells_by = "wsnn_res_0.5_cell_type")+ggtitle('raw.umap')
p2
p=p1|p2
#ggsave("combine_monocle3_raw_umap.pdf",p,width = 9.5)

##cluster_cells
cds <- cluster_cells(cds)

##learn_graph
cds <- learn_graph(cds)

p3 = plot_cells(cds, color_cells_by = "wsnn_res_0.5_cell_type", 
                label_groups_by_cluster=FALSE,label_leaves=FALSE, label_branch_points=FALSE)
#ggsave("Monocle3_Trajectory_withlabel.pdf", p3,width = 6,height = 5)


p4 <-plot_cells(cds,color_cells_by="wsnn_res_0.5_cell_type",label_groups_by_cluster=FALSE,
                label_leaves=TRUE,label_branch_points=FALSE)
#ggsave("Monocle3_trajectory_leaves_withoutlabel.pdf",p4,width = 6,height = 5)


##choose the root
cds <- order_cells(cds)
p5<-plot_cells(cds,color_cells_by="pseudotime",label_groups_by_cluster=FALSE,
               label_leaves=FALSE,label_branch_points=FALSE)+
  scale_color_gradient(low="#B0CFE4", high = "#B22028")
#ggsave("monocle3_trajectory_pseudotime.pdf",p5,width = 6,height = 20)

##save the result
saveRDS(cds, file = "scRNA_cds.rds")
