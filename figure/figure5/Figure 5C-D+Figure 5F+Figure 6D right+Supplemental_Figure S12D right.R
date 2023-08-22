###Figure 5C-D+Figure 5F+Figure 6D right+Supplemental_Figure S12D right
library(BiocManager)
library(Seurat)
library(monocle3)
library(SingleCellExperiment)
library(tidyverse)
library(patchwork)
library(ggplot2)
set.seed(123)

##organoids singlecell data
single_data<-readRDS("C:/Users/Lenovo/Desktop/regeneration/Final_Reg_CMZ_0807.rds")

##import the result from organoids_monocle3.R
cds<-readRDS("organoids_monocle3_cds.rds")

##choose hRSLCs as root
cds <- order_cells(cds)

##Figure 5C
p5<-plot_cells(cds,color_cells_by="pseudotime",label_groups_by_cluster=FALSE,
               label_leaves=FALSE,label_branch_points=FALSE)+
  scale_color_gradient(low="#B0CFE4", high = "#B22028")
ggsave("monocle3_trajectory_pseudotime.pdf",p5,width = 6,height = 5)


##Figure 6D right+Supplemental_Figure S12D right
Track_genes_sig<-c("MECOM","MEIS1","TBX20","FOXP1","COL9A1","RELN","CPAMD8")
p<-plot_genes_in_pseudotime(cds[Track_genes_sig,], color_cells_by="pseudotime",min_expr=0.01, ncol = 2,cell_size = 0)+
  scale_color_gradient(low="#B0CFE4", high = "#B22028")
ggsave("regeneration_monocle3_Genes_pseudotime1.pdf", plot = p, width = 6, height = 20)


##pseudotime
pseudotime<-pseudotime(cds, reduction_method ="UMAP")
metadata<-single_data@meta.data

metadata$pseudotime<-pseudotime
metadata<-metadata[,c("final_annotated_0807_2","pseudotime")]

##subset four celltypes
metadata<-metadata[metadata$final_annotated_0807_2=="hRSLCs"|metadata$final_annotated_0807_2=="RPCs"|metadata$final_annotated_0807_2=="PC_precursors"|metadata$final_annotated_0807_2=="PCs",]

library(dplyr)
metadata <- filter(metadata, pseudotime!= Inf)
metadata$final_annotated_0807_2<-factor(metadata$final_annotated_0807_2,levels = c("hRSLCs","RPCs","PC_precursors","PCs"))
colnames(metadata)[1]<-"Cell_types"

##Figure 5D
library(ggplot2)
p1 <- ggplot(metadata, aes(x=pseudotime,fill=Cell_types)) + 
  theme_classic()+
  geom_density()+
  scale_fill_manual(values = c("#af2157","#BD8098","#C0BFDF","#67ADB7"))+
  scale_y_sqrt(limits = c(0,0.5),breaks = c(0,0.1,0.2,0.3))+
  theme(axis.title = element_text(size=18),
        axis.text = element_text(size=14),
        axis.text.x = element_text(angle = 35,
                                   hjust=1,
                                   size = 10,
                                   color = "black"))+
  xlim(0,30)
p1
ggsave(filename = "density_plot.pdf",p1)


##Figure 5F
##pseudotime
pseudotime<-pseudotime(cds, reduction_method ="UMAP")
metadata<-single_data@meta.data

metadata$pseudotime<-pseudotime
metadata<-metadata[,c("treatment","pseudotime")]

library(dplyr)
metadata <- filter(metadata, pseudotime!= Inf)

metadata$treatment<-factor(metadata$treatment,levels = c("CMZ-0","CMZ-05", "CMZ-10", "CMZ-20", "CMZ-40"))

#boxplot3
p = ggplot(metadata, aes(x=treatment, y=as.numeric(pseudotime))) + 
  stat_boxplot(geom = "errorbar",width=0.2, size=0.5,position=position_dodge(0.6),color= "black")+
  theme_classic()+
  geom_boxplot(position = position_dodge(0.6),
               size = 0.5,
               width = 0.8,
               fill = c("#af2157","#BD8098","#C0BFDF","#CEE2EC","#67ADB7"),
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
  ylab("Pseudotimes")+
  theme(axis.title = element_text(size=18),
        axis.text = element_text(size=14),
        axis.text.x = element_text(angle = 35,
                                   hjust=1,
                                   size = 10))
p








