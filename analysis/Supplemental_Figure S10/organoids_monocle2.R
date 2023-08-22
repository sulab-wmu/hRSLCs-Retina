##organoids monocle2
library(BiocManager)
library(Seurat)
library(monocle)
library(tidyverse)
library(patchwork)
library(sctransform)
library(S4Vectors)
library(Signac)
library(dplyr)

setwd("/share2/pub/chenchg/chenchg/SingleCell/LiHui/regeneration/monocle2/")
##organoids singlecell data
data<-readRDS("/share2/pub/chenchg/chenchg/SingleCell/LiHui/regeneration/Final_Reg_CMZ_0807.rds")

expr_matrix<-as(as.matrix(data@assays$RNA@counts),'sparseMatrix')
p_data<-data@meta.data
p_data$final_annotated_0807_2 = as.character(p_data$final_annotated_0807_2)
f_data <- data.frame(gene_short_name = row.names(data@assays$RNA@counts),row.names = row.names(data@assays$RNA@counts))

##create cds object
pd <- new('AnnotatedDataFrame', data = p_data)
fd <- new('AnnotatedDataFrame', data = f_data)
cds <- newCellDataSet(expr_matrix,
                      phenoData = pd,
                      featureData = fd,
                      lowerDetectionLimit = 0.5,
                      expressionFamily = negbinomial.size())


##estimateSizeFactors,estimateDispersions
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

##input the DEG
disp_table <- dispersionTable(cds)
disp.genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
diff <- differentialGeneTest(cds[disp.genes,],fullModelFormulaStr="~final_annotated_0807_2",cores=16)

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
saveRDS(cds,file = "organoids_monocle2_cds.rds")