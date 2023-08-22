#-------------------------Fetal four main cell types-----------------------------Figure 1H
##import singlecell data from fetal_celltypes_annotation.R
data_fetal_scRNA<-readRDS("C:/Users/Lenovo/Desktop/final_version/annotation/final_annotation.rds")

##import CytoTRACE result from fetal_CytoTRACE.R
data_cytoTRACE<-readRDS("C:/Users/Lenovo/Desktop/final_version/CytoTRACE/CytoTRACE.rds")
data_fetal_scRNA$cytoTRACE <- data_cytoTRACE$CytoTRACE


Idents(data_fetal_scRNA)<-data_fetal_scRNA$wsnn_res_0.5_cell_type
table(Idents(data_fetal_scRNA))
subset_fetal_scRNA <- subset(data_fetal_scRNA,idents = c("hRSLCs","RPCs","PC_precursors","PCs"))
  
table(subset_fetal_scRNA$wsnn_res_0.5_cell_type)
  
#####################################################---------------Start-------Figure 1H--Left panel
#Re-PCA and UMAP for all severe
subset_fetal_scRNA = FindVariableFeatures(subset_fetal_scRNA, do.plot = F, display.progress = FALSE)
subset_fetal_scRNA = ScaleData(subset_fetal_scRNA, display.progress = FALSE)
subset_fetal_scRNA <- RunPCA(subset_fetal_scRNA, verbose = FALSE)
ElbowPlot(subset_fetal_scRNA, ndims = 10)
subset_fetal_scRNA <- FindNeighbors(object = subset_fetal_scRNA, dims = 1:10)
subset_fetal_scRNA <- FindClusters(object = subset_fetal_scRNA, resolution = 0.1)
subset_fetal_scRNA <- RunUMAP(object = subset_fetal_scRNA, dims = 1:10)
  
#plot by using UMAP
Idents(subset_fetal_scRNA)<-subset_fetal_scRNA$wsnn_res_0.5_cell_type
DimPlot(subset_fetal_scRNA, reduction = "umap")
  
##change the color
pal <-c("#af2157","#BD8098","#C0BFDF","#67ADB7")
Idents(subset_fetal_scRNA)<-subset_fetal_scRNA$wsnn_res_0.5_cell_type
DimPlot(subset_fetal_scRNA, reduction = "umap",label = F,  cols= pal, pt.size = 1, repel = T)
  
  
######################---------------------------------------------------Figure 1H--Middle panel
library(tidyverse)
library(RColorBrewer)
W2<-FeaturePlot(subset_fetal_scRNA, reduction = "umap",features = "cytoTRACE", label = FALSE, repel = TRUE) +
    scale_colour_gradient2(high = "#DC0000FF",mid = "#EDE361",low = "#69B9DA",midpoint = 0.3)
W2
  

#################four cell types in Lineage 1 for CytoTRACE boxplot
cytoTRACE <- as.numeric(data_fetal_scRNA$cytoTRACE)
celltypes <- as.character(data_fetal_scRNA$wsnn_res_0.5_cell_type)
data_cytoTRACE <- as.data.frame(cbind(celltypes,cytoTRACE))
names(data_cytoTRACE) <- c("celltypes_0.05","cytoTRACE_score")

data_cytoTRACE_subset <- data_cytoTRACE[which(data_cytoTRACE$celltypes %in% c("hRSLCs","PC_precursors","PCs","RPCs")),]
head(data_cytoTRACE_subset)
  
  
  
###########-----------------------------------------------------Boxplot------------------Figure 1H right panel
data_cytoTRACE_subset$cytoTRACE_score
data_cytoTRACE_subset$celltypes_0.05<-factor(data_cytoTRACE_subset$celltypes_0.05,levels = c("hRSLCs","RPCs","PC_precursors","PCs"))
  
library(ggplot2)
#boxplot3
p = ggplot(data_cytoTRACE_subset, aes(x=celltypes_0.05, y=as.numeric(cytoTRACE_score))) + 
    stat_boxplot(geom = "errorbar",width=0.2, size=0.5,position=position_dodge(0.6),color= "black")+
    theme_classic()+
    geom_boxplot(position = position_dodge(0.6),
                 size = 0.5,
                 width = 0.8,
                 fill = c("#af2157","#BD8098","#C0BFDF","#67ADB7"),
                 color = "black",
                 outlier.color = "black",
                 outlier.fill = "black",
                 outlier.shape = 19,
                 outlier.size = 1.5,
                 outlier.stroke = 0.5,
                 outlier.alpha = 45,
                 notch = F,
                 notchwidth = 0.5)+
    xlab("")+
    ylab("Stemness score")+
    theme(axis.title = element_text(size=18),
          axis.text = element_text(size=14),
          axis.text.x = element_text(angle = 35,
                                     hjust=1,
                                     size = 10))


