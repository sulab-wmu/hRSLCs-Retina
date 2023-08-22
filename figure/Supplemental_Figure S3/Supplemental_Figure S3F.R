#####Supplemental_Figure S3F

##import singlecell data from fetal_celltypes_annotation.R
data_fetal_scRNA<-readRDS("C:/Users/Lenovo/Desktop/final_version/annotation/final_annotation.rds")

##import the result from fetal_Stem_ID.R
load("C:/Users/Lenovo/Desktop/final_version/Stem_ID/entropy.Rdata")
data_fetal_scRNA$entropy <- entropy

Idents(data_fetal_scRNA)<-data_fetal_scRNA$wsnn_res_0.5_cell_type
table(Idents(data_fetal_scRNA))
subset_fetal_scRNA <- subset(data_fetal_scRNA,idents = c("hRSLCs","RPCs","PC_precursors","PCs"))

table(subset_fetal_scRNA$wsnn_res_0.5_cell_type)

#####################################################---------------Start---------Supplemental_Figure S3F Left panel
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

library(tidyverse)
library(RColorBrewer)
W2<-FeaturePlot(subset_fetal_scRNA, reduction = "umap",features = "entropy", label = FALSE, repel = TRUE) +
  scale_colour_gradient2(high = "#DC0000FF",mid = "#EDE361",low = "#69B9DA",midpoint = 0.6)
W2


#####################################################---------------Start---------Supplemental_Figure S3F Right panel
entropy <- as.numeric(data_fetal_scRNA$entropy)
celltypes <- as.character(data_fetal_scRNA$wsnn_res_0.5_cell_type)
data_entropy <- as.data.frame(cbind(celltypes,entropy))
names(data_entropy) <- c("celltypes_0.05","entropy")

data_entropy_subset <- data_entropy[which(data_entropy$celltypes %in% c("hRSLCs","PC_precursors","PCs","RPCs")),]
head(data_entropy_subset)

data_entropy_subset$entropy
data_entropy_subset$celltypes_0.05<-factor(data_entropy_subset$celltypes_0.05,levels = c("hRSLCs","RPCs","PC_precursors","PCs"))

library(ggplot2)
#boxplot3
p = ggplot(data_entropy_subset, aes(x=celltypes_0.05, y=as.numeric(entropy))) + 
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
p






