##Supplemental_Figure S1A-D
library(Seurat)
library(Signac)

##import singlecell data from fetal_celltypes_annotation.R
single_data<-readRDS("C:/Users/Lenovo/Desktop/final_version/annotation/final_annotation.rds")

##change the color
colors<-c("#71C89C","#67ADB7","#36600E", "#6A8473", "#C0BFDF","#E77A77","#7B6148","#ECBA84","#CA8C74","#A13B46","#9569AB")

##Figure S1A
##RNA only
DimPlot(single_data,group.by = "wsnn_res_0.5_cell_type1",reduction = "umap.rna",cols = colors)

##ATAC only
DimPlot(single_data,group.by = "wsnn_res_0.5_cell_type1",reduction = "umap.atac",cols = colors)

##combine RNA and ATAC
DimPlot(single_data,group.by = "wsnn_res_0.5_cell_type1",reduction = "umap.biMod",cols = colors)

##Figure S1B top
table(single_data$wsnn_res_0.5_cell_type1)
cell_number<-as.data.frame(table(single_data$wsnn_res_0.5_cell_type1))

library(ggplot2)
library(cowplot)
ggplot(cell_number,aes(Var1,Freq))+
  geom_col(aes(fill=Var1),width = 0.7)+
  scale_fill_manual(values = colors)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        plot.title = element_text (hjust = 0.5),
        axis.text.x = element_text (face="bold", size=10),
        legend.position = "none")+
  xlab("")+
  ylab("Cell counts")+
  geom_text(aes(label = Freq), vjust = -0.5)+
  scale_y_sqrt(limits = c(0,5000),breaks = c(0,50,100,200,1000,2000,4000,5000))
  


##Figure S1B bottom
Idents(single_data)<-single_data$orig.ident
Cellratio <- prop.table(table(Idents(single_data), single_data$wsnn_res_0.5_cell_type1), margin = 2)
Cellratio <- as.data.frame(Cellratio)
Cellratio$Var1<-factor(Cellratio$Var1,levels = c("1_CMZ_1","1_NR_1"))

ggplot(Cellratio) + 
  geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity",width = 0.7,size = 0.5,colour = NA)+ 
  scale_fill_manual(values = c("#E77A77","#67ADB7"))+
  theme_cowplot()+
  theme(panel.border = element_rect(fill=NA,color=NA, size=0.5, linetype="solid"),
        axis.text.x = element_text(angle = 35, vjust = 0.85,size = 9))+
  ylab("Percentage")+
  xlab("")+
  scale_y_continuous(expand = c(0,0))



##Figure S1C left
##subset cells belong to CMZ
Idents(single_data)<-single_data$orig.ident
single_data_CMZ<-subset(single_data,idents = "1_CMZ_1")
DimPlot(single_data_CMZ,group.by = "wsnn_res_0.5_cell_type1",reduction = "umap.biMod",cols = colors)

##Figure S1C right
cell_number_CMZ<-as.data.frame(table(single_data_CMZ$wsnn_res_0.5_cell_type1))
cell_number_CMZ<-cell_number_CMZ[-11,]
ggplot(cell_number_CMZ,aes(Var1,Freq))+
  geom_bar(aes(fill=Var1),width = 0.7,stat = "identity")+
  scale_fill_manual(values = colors)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        plot.title = element_text (hjust = 0.5),
        axis.text.x = element_text (face="bold", size=10),
        legend.position = "none")+
  xlab("")+
  ylab("Cell counts")+
  geom_text(aes(label = Freq),hjust=-0.3,vjust =0)+
  coord_flip()+
  scale_y_sqrt(limits = c(0,3000),breaks = c(0,50,100,200,1000,2000,3000))






##Figure S1D
##subset cells belong to NR
single_data_NR<-subset(single_data,idents = "1_NR_1")
colors1<-c("#71C89C","#67ADB7", "#6A8473", "#C0BFDF","#E77A77","#7B6148","#ECBA84","#CA8C74","#A13B46","#9569AB")
DimPlot(single_data_NR,group.by = "wsnn_res_0.5_cell_type1",reduction = "umap.biMod",cols = colors1)

##Figure S1D right
cell_number_NR<-as.data.frame(table(single_data_NR$wsnn_res_0.5_cell_type1))
cell_number_NR<-cell_number_NR[-3,]
ggplot(cell_number_NR,aes(Var1,Freq))+
  geom_bar(aes(fill=Var1),width = 0.7,stat = "identity")+
  scale_fill_manual(values = colors1)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        plot.title = element_text (hjust = 0.5),
        axis.text.x = element_text (face="bold", size=10),
        legend.position = "none")+
  xlab("")+
  ylab("Cell counts")+
  geom_text(aes(label = Freq),hjust=-0.3,vjust =0)+
  coord_flip()+
  scale_y_sqrt(limits = c(0,3000),breaks = c(0,50,100,200,1000,2000,3000))















