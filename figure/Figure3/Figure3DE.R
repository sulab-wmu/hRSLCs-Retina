library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(ggpubr)
scdata<-readRDS("H9_D60_2_stdata.rds")
pdf("Figure3D.pdf")
SpatialDimPlot(scdata)
dev.off()
scdata<-subset(scdata,subset = seurat_clusters!=6 & seurat_clusters!=7)
#############
#load marker genes
###########
read.csv("cellmarker_important_select.csv",header = T)->marker_check
scdata$seurat_clusters<-as.character(scdata$seurat_clusters)
scdata$seurat_clusters[scdata$seurat_clusters==0]<-6
scdata$New_celltype<-scdata$seurat_clusters
scdata$New_celltype[scdata$New_celltype=="6"]<-"RPC/PC"
scdata$New_celltype[scdata$New_celltype=="1"]<-"RGC/RPC"
scdata$New_celltype[scdata$New_celltype=="2"]<-"RPC/PC_precursors/MC"
scdata$New_celltype[scdata$New_celltype=="3"]<-"RPE/RPE_progenitors"
scdata$New_celltype[scdata$New_celltype=="4"]<-"hRSLCs"
scdata$New_celltype[scdata$New_celltype=="5"]<-"RPC/PC/AC_HC/BC/MC"

saveRDS(scdata,"H9_D60_2_final_stdata.rds")
metadata<-scdata@meta.data
write.csv(metadata,"H9_D60_2_final_metadata.csv",row.names = T)

cell_names<-c("RGC/RPC","RPC/PC_precursors/MC","RPE/RPE_progenitors","hRSLCs","RPC/PC/AC_HC/BC/MC","RPC/PC")
scdata$New_celltype<-factor(scdata$New_celltype,levels = rev(cell_names))
GENES<-intersect(marker_check$GENES,rownames(scdata))
pdf("Figure3E_dotplot.pdf",width = 14,height = 4)
DotPlot(scdata,features=unique(GENES),group.by="New_celltype",assay="SCT")+
      scale_colour_gradient2(low="#3A71AA",mid="white",high="#B22028",midpoint=0)+
      theme_bw()+
      theme (axis.text.x = element_text (angle = 45, hjust = 1))+theme(legend.position = "top")
dev.off()