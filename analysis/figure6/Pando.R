##Pando code
##devtools::install_github('quadbio/Pando')
library(Pando)
library(Seurat)
library(BSgenome.Hsapiens.UCSC.hg38)
library(Signac)

##import singlecell data from fetal_celltypes_annotation.R
single_data<-readRDS("C:/Users/Lenovo/Desktop/final_version/annotation/final_annotation.rds")
DimPlot(single_data,group.by = "wsnn_res_0.5_cell_type",reduction = "umap.biMod",label = T,label.size = 3)

##import motifs
data(motifs)

##FindVariableFeatures
single_data <- Seurat::FindVariableFeatures(single_data, assay='RNA')

##initiate_grn
single_data <- initiate_grn(single_data)

##find_motifs
single_data <- find_motifs(single_data,
                           pfm = motifs,
                           genome = BSgenome.Hsapiens.UCSC.hg38)

##infer_grn
single_data <- infer_grn(single_data,
                         peak_to_gene_method ="Signac",
                         method = 'glm')

##save the result
save(single_data,file = "Pando.rds")


