##enrichment analysis
library(Seurat)
library(Seurat)
library(tidyverse)
library(clusterProfiler)
library(pathview)
library(enrichplot)
library(msigdbr)

##import the singlecell data from fetal_celltypes_annotation.R
single_data<-readRDS("C:/Users/Lenovo/Desktop/final_version/annotation/final_annotation.rds")
Idents(single_data)<-single_data$wsnn_res_0.5_cell_type
DimPlot(single_data,group.by ="wsnn_res_0.5_cell_type",reduction = "umap.biMod" )

##hRSLCs
hRSLCs_markers<-FindMarkers(single_data,assay = "RNA",ident.1 = "hRSLCs",only.pos = T,min.pct = 0.25, logfc.threshold = 0.25)
hRSLCs_markers<-hRSLCs_markers[hRSLCs_markers$p_val_adj<0.05,]

##GO BP
GO_BP <- enrichGO(gene = row.names(hRSLCs_markers),
                  OrgDb = 'org.Hs.eg.db',
                  keyType = 'SYMBOL',
                  ont = "BP", 
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.01,
                  qvalueCutoff = 0.05)
GO_BP_result<-GO_BP@result
GO_BP_result<-GO_BP_result[GO_BP_result$p.adjust<0.05,]
write.csv(GO_BP_result,file = "C:/Users/Lenovo/Desktop/github_code/analysis/hRSLCs_GO_BP_result.csv",quote = F,row.names = F)


##RPE_progenitors
RPE_progenitors_markers<-FindMarkers(single_data,assay = "RNA",ident.1 = "RPE_progenitors",only.pos = T,min.pct = 0.25, logfc.threshold = 0.25)
RPE_progenitors_markers<-RPE_progenitors_markers[RPE_progenitors_markers$p_val_adj<0.05,]

GO_BP <- enrichGO(gene = row.names(RPE_progenitors_markers),
                  OrgDb = 'org.Hs.eg.db',
                  keyType = 'SYMBOL',
                  ont = "BP", 
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.01,
                  qvalueCutoff = 0.05)
GO_BP_result<-GO_BP@result
GO_BP_result<-GO_BP_result[GO_BP_result$p.adjust<0.05,]
write.csv(GO_BP_result,file = "C:/Users/Lenovo/Desktop/github_code/analysis/RPE_progenitors_GO_BP_result.csv",quote = F,row.names = F)
