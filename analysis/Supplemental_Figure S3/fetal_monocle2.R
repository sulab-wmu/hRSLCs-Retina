##monocle2 
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

##import singlecell data from fetal_celltypes_annotation.R
data<-readRDS("/share2/pub/chenchg/chenchg/SingleCell/LiHui/final/final_annotation.rds")

expr_matrix<-as(as.matrix(data@assays$RNA@counts),'sparseMatrix')
p_data<-data@meta.data
p_data$wsnn_res_0.5_cell_type = as.character(p_data$wsnn_res_0.5_cell_type)
f_data <- data.frame(gene_short_name = row.names(data@assays$RNA@counts),row.names = row.names(data@assays$RNA@counts))

##create cds object
pd <- new('AnnotatedDataFrame', data = p_data)
fd <- new('AnnotatedDataFrame', data = f_data)
cds <- newCellDataSet(expr_matrix,
                      phenoData = pd,
                      featureData = fd,
                      lowerDetectionLimit = 0.5,
                      expressionFamily = negbinomial.size())


##estimateSizeFactors and estimateDispersions
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

##input the DEG
deg.cluster <- FindAllMarkers(data)
diff.genes <- subset(deg.cluster,p_val_adj<0.05)$gene
diff.genes <- unique(diff.genes)
diff <- differentialGeneTest(cds[diff.genes,],fullModelFormulaStr="~wsnn_res_0.5_cell_type",cores=16)

deg <- subset(diff, qval < 0.01) 
deg <- deg[order(deg$qval,decreasing=F),]
ordergene <- rownames(deg)

##setOrderingFilter
cds <- setOrderingFilter(cds, ordergene)
##reduceDimension
cds <- reduceDimension(cds, max_components = 2, method = 'DDRTree')
##orderCells
cds <- orderCells(cds)

##save the result
saveRDS(cds,file = "cds.rds")

