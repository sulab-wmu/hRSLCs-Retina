##Figure 5A-B+Supplemental_Figure S10A-C
library(Seurat)

##organoids singlecell data
single_data<-readRDS("C:/Users/Lenovo/Desktop/regeneration/Final_Reg_CMZ_0807.rds")

table(single_data$final_annotated_0807)
table(single_data$final_annotated_0807_2)
table(single_data$treatment)

DefaultAssay(single_data)<-"RNA"
Idents(single_data)<-single_data$final_annotated_0807_2
single_data$final_annotated_0807_2<-factor(single_data$final_annotated_0807_2,levels = c("RPE","RPE_progenitors","hRSLCs","RPCs","PC_precursors",
                                                                                         "PCs","RGCs","ACs","HCs","BCs"))

##Figure 5A
pal <- c("#85C89C","#67ADB7","#36600E","#6A8473","#C0BFDF","#E77A77","#7E6148FF","#ECBA84","#CA8C74","#A13B46")
pdf("regeneration umap.pdf")
DimPlot(single_data,group.by = "final_annotated_0807_2",reduction = "umap",cols = pal,label = T,label.size = 3)
dev.off()


##Figure 5B
library(openxlsx)
marker<-read.xlsx("C:/Users/Lenovo/Desktop/github_code/analysis/figure5/Markers_organoids.xlsx")

library(ggplot2)
pdf("regeneration marker dotplot.pdf")
DotPlot(single_data,features=rev(unique(marker$GENES)),group.by="final_annotated_0807_2",assay="SCT")+
  scale_colour_gradient2(low="#3A71AA",mid="white",high="#B22028",midpoint=0)+
  theme_bw()+
  coord_flip()+
  theme (axis.text.x = element_text (angle = 35, hjust = 1))
dev.off()


##Supplemental_Figure S10A
library(Nebulosa)
library(BiocFileCache)
library(paletteer)
library(scCustomize)
##hRSLCs
Plot_Density_Custom(seurat_object =single_data, features = "MECOM",reduction = "umap",
                    custom_palette = c("#B0CFE4","#FACABC","#E77A77","#DC0000FF"))
Plot_Density_Custom(seurat_object =single_data, features = "CPAMD8",reduction = "umap",
                    custom_palette = c("#B0CFE4","#FACABC","#E77A77","#DC0000FF"))
Plot_Density_Custom(seurat_object =single_data, features = "COL9A1",reduction = "umap",
                    custom_palette = c("#B0CFE4","#FACABC","#E77A77","#DC0000FF"))


##RPCs
Plot_Density_Custom(seurat_object =single_data, features = "SPP1",reduction = "umap",
                    custom_palette = c("#B0CFE4","#FACABC","#E77A77","#DC0000FF"))
Plot_Density_Custom(seurat_object =single_data, features = "SOX2",reduction = "umap",
                    custom_palette = c("#B0CFE4","#FACABC","#E77A77","#DC0000FF"))


##PC_precursors
Plot_Density_Custom(seurat_object =single_data, features = "ATOH7",reduction = "umap",
                    custom_palette = c("#B0CFE4","#FACABC","#E77A77","#DC0000FF"))

##PCs
Plot_Density_Custom(seurat_object =single_data, features = "RCVRN",reduction = "umap",
                    custom_palette = c("#B0CFE4","#FACABC","#E77A77","#DC0000FF"))
Plot_Density_Custom(seurat_object =single_data, features = "RP1",reduction = "umap",
                    custom_palette = c("#B0CFE4","#FACABC","#E77A77","#DC0000FF"))

##RPE_progenitors
Plot_Density_Custom(seurat_object =single_data, features = "GJA1",reduction = "umap",
                    custom_palette = c("#B0CFE4","#FACABC","#E77A77","#DC0000FF"))
Plot_Density_Custom(seurat_object =single_data, features = "MITF",reduction = "umap",
                    custom_palette = c("#B0CFE4","#FACABC","#E77A77","#DC0000FF"))

##RPE
Plot_Density_Custom(seurat_object =single_data, features = "CD96",reduction = "umap",
                    custom_palette = c("#B0CFE4","#FACABC","#E77A77","#DC0000FF"))
Plot_Density_Custom(seurat_object =single_data, features = "BEST1",reduction = "umap",
                    custom_palette = c("#B0CFE4","#FACABC","#E77A77","#DC0000FF"))


##Supplemental_Figure S10B
colors<-c("#9A3B45","#DC7774","#E6B682","#C8D7AA","#6A8372")
single_data$treatment<-factor(single_data$treatment,levels=c("CMZ-0","CMZ-05","CMZ-10","CMZ-20","CMZ-40"))
DimPlot(single_data,group.by = "treatment",reduction = "umap",cols = colors,label = F,label.size = 3)

##Supplemental_Figure S10C
Idents(single_data)<-single_data$treatment
Cellratio <- prop.table(table(Idents(single_data), single_data$final_annotated_0807_2),margin = 2)
Cellratio <- as.data.frame(Cellratio)

library(ggplot2)
library(cowplot)
ggplot(Cellratio) + 
  geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity",width = 0.7,size = 0.5,colour = NA)+ 
  scale_fill_manual(values = colors)+
  theme_cowplot()+
  theme(panel.border = element_rect(fill=NA,color=NA, size=0.5, linetype="solid"),
        axis.text.x = element_text(angle = 35, vjust = 0.85,size = 9))+
  ylab("Percentage")+
  xlab("")+
  scale_y_continuous(expand = c(0,0))








