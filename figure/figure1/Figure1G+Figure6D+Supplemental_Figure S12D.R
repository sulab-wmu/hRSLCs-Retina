##fetal scRNA monocle3 plot 
library(BiocManager)
library(Signac)
library(BiocGenerics)
library(monocle3)
library(SingleCellExperiment)
library(tidyverse)
library(patchwork)
library(ggplot2)

##input cds from fetal_scRNA_monocle3.R
cds<-readRDS("scRNA_cds.rds")

##choose the root
cds <- order_cells(cds)
p5<-plot_cells(cds,color_cells_by="pseudotime",label_groups_by_cluster=FALSE,
               label_leaves=FALSE,label_branch_points=FALSE)+
  scale_color_gradient(low="#B0CFE4", high = "#B22028")

##save the monocle3_trajectory_pseudotime plot
ggsave("monocle3_trajectory_pseudotime.pdf",p5,width = 6,height = 20)

##show how key genes change over pseudotime
Track_genes_sig<-c("MECOM","MEIS1","TBX20","FOXP1","COL9A1","RELN","CPAMD8")
p6<-plot_genes_in_pseudotime(cds[Track_genes_sig,], color_cells_by="pseudotime",min_expr=0.01, ncol = 2,cell_size = 0)+
  scale_color_gradient(low="#B0CFE4", high = "#B22028")
ggsave("fetal_monocle3_Genes_pseudotime.pdf", plot = p6, width = 6, height =20 )