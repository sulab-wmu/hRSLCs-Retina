####--------hRSLCs-------------------------------------------------------Figure 5H -----
library(Seurat)
library(ggplot2)
##organoids singlecell data  
data_regen_scRNA<-readRDS("C:/Users/Lenovo/Desktop/regeneration/Final_Reg_CMZ_0807.rds")

Idents(data_regen_scRNA) <- data_regen_scRNA$final_annotated_0807_2
subset_regen_hRSLCs <- subset(data_regen_scRNA,idents = c("hRSLCs"))

Idents(subset_regen_hRSLCs) <- subset_regen_hRSLCs$treatment
table(Idents(subset_regen_hRSLCs))

DefaultAssay(subset_regen_hRSLCs)<-"RNA"
markers_subset_regen_hRSLCs_CMZ5_vs_0 <- FindMarkers(object = subset_regen_hRSLCs, ident.1 ="CMZ-05",
                                                     ident.2 = "CMZ-0", 
                                                     min.pct = 0,
                                                     logfc.threshold=0)
print(x = head(markers_subset_regen_hRSLCs_CMZ5_vs_0))
write.csv(markers_subset_regen_hRSLCs_CMZ5_vs_0,"markers_subset_regen_hRSLCs_CMZ5_vs_02_based_on_RNA.csv")

