##monocle3 scATAC
#install.packages("devtools")
#devtools::install_github("cole-trapnell-lab/cicero-release", ref = "monocle3")
library(cicero)
library(monocle3)
library(Seurat)
library(tidyverse)
library(patchwork)
library(SingleCellExperiment)
set.seed(1234)

setwd("/share2/pub/chenchg/chenchg/SingleCell/LiHui/final/monocle3_scATAC/")

##import the singlecell from fetal_celltypes_annotation.R
single_data<-readRDS("/share2/pub/chenchg/chenchg/SingleCell/LiHui/final/final_annotation.rds")

##create cds object
DefaultAssay(single_data)<-"peaks"
data <- GetAssayData(single_data, assay = 'peaks', slot = 'counts')
cell_metadata <- single_data@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,cell_metadata = cell_metadata,gene_metadata = gene_annotation)


cds <- monocle3::detect_genes(cds)
cds <- cds[Matrix::rowSums(exprs(cds)) != 0,]


cds <- detect_genes(cds)
cds <- estimate_size_factors(cds)

##LSI
cds <- preprocess_cds(cds, method = "LSI")
##UMAP
cds <- reduce_dimension(cds, reduction_method = 'UMAP', preprocess_method = "LSI")

cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(single_data, reduction = "umap.biMod")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed
plot_cells(cds, reduction_method="UMAP", color_cells_by="wsnn_res_0.5_cell_type")

##cluster_cells
cds <- cluster_cells(cds)
##learn_graph
cds <- learn_graph(cds)

p3 = plot_cells(cds, color_cells_by = "wsnn_res_0.5_cell_type", 
                label_groups_by_cluster=FALSE,label_leaves=FALSE, label_branch_points=FALSE)
ggsave("Monocle3_Trajectory_withlabel.pdf", p3,width = 6,height = 5)

p4 <-plot_cells(cds,color_cells_by="wsnn_res_0.5_cell_type",label_groups_by_cluster=FALSE,
                label_leaves=TRUE,label_branch_points=FALSE)
ggsave("Monocle3_trajectory_leaves_withoutlabel.pdf",p4,width = 6,height = 5)

##save the result
saveRDS(cds,file = "scATAC_cds.rds")