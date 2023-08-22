##Figure 4E
library(openxlsx)

##import the gene expression profile about regeneration
regeneration<-read.xlsx("regeneration gene.xlsx")
rownames(regeneration)<-regeneration[,1]
regeneration<-regeneration[,-1]
row_annotation<-cbind(rownames(regeneration),regeneration[,7])
regeneration<-regeneration[,-7]

regeneration<-regeneration[,-1]
regeneration<-t(scale(t(regeneration),center = T,scale = T))

##pheatmep
library(ggplot2)
library(pheatmap)

range(regeneration)
bk<-c(seq(-2,-0.1,by=0.01),seq(0,2,by=0.01))

##row_annotation
rownames(row_annotation)<-row_annotation[,1]
row_annotation<-as.data.frame(row_annotation[,-1])
colnames(row_annotation)<-"Cell_type"

##change the color
ann_colors<-list(Cell_type=c("RPC"="#ECBA84","RGC"="#67ADB7", "HC/AC"="#CA8C74","BC"="#9569AB","pan-PRC"="#C0BFDF",
                             "Cone"="#E77A77","Rod"="#6A8473","MC"="#71C89C"))

p1<-pheatmap(regeneration,cluster_cols = F,cluster_rows = F,
             show_rownames = F,show_colnames = T,
             cellwidth = 30,cellheight = 1.5,
             color = c(colorRampPalette(colors = c("#3A71AA","white"))(length(bk)/2),colorRampPalette(colors = c("white","#B22028"))(length(bk)/2)),
             breaks = bk,fontsize = 7,
             annotation_row =row_annotation,
             annotation_colors = ann_colors,
             angle_col = 90
)

ggsave(p1,file="regeneration gene heatmap.pdf")