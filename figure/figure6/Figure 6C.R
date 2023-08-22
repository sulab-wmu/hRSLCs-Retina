##Figure 6C
##Pando plot
library(Pando)
library(Seurat)
library(BSgenome.Hsapiens.UCSC.hg38)
library(Signac)
library(tidyverse)
library(patchwork)

##import the result from Pando.R
single_data<-readRDS("C:/Users/Lenovo/Desktop/final_version/Pando/Pando.rds")

##umap 
single_data <- get_network_graph(single_data, graph_name='umap_graph',umap_method = "corr",random_seed = 123)


pdf("network_with_label.pdf",height = 15,width = 15)
plot_network_graph(single_data, graph='umap_graph',
                   layout = "umap",
                   edge_color="grey",
                   node_color = "#6A8372",
                   node_size = c(5,5.1),
                   edge_width = 0.01,
                   color_nodes = T,
                   label_nodes = T
)
dev.off()


pdf("network_without_label.pdf",height = 15,width = 15)
plot_network_graph(single_data, graph='umap_graph',
                   layout = "umap",
                   edge_color="grey",
                   node_color = "#6A8372",
                   node_size = c(5,5.1),
                   edge_width = 0.01,
                   color_nodes = T,
                   label_nodes = F
)
dev.off()

single_data <- get_network_graph(
  single_data, 
  graph_name = 'full_graph', 
  umap_method = 'none'
)

##visualization key TF
##MECOM
single_data <- get_tf_network(single_data, tf="MECOM", graph='full_graph', keep_all_edges = T)
plot_tf_network(
  single_data, 
  tf = "MECOM",
  edge_color = "#E77A77",
  circular = T,
  edge_width = 1,
  label_nodes = "all"
)

##TBX20
single_data <- get_tf_network(single_data, tf="TBX20", graph='full_graph', keep_all_edges = T)
plot_tf_network(
  single_data, 
  tf = "TBX20",
  edge_color = "#E77A77",
  circular = T,
  edge_width = 1,
  label_nodes = "all"
)



