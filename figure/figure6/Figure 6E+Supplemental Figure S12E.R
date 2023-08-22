###################################-----------------------------------------Figure 6E+Supplemental Figure S12E

##organoids singlecell data
data_regen_scRNA<-readRDS("C:/Users/Lenovo/Desktop/regeneration/Final_Reg_CMZ_0807.rds")
DefaultAssay(data_regen_scRNA)<-"SCT"

library(tidyverse)
library(RColorBrewer)
library(Seurat)
library(ggplot2)
mecom_targeted_genes <- c("COL9A1","RELN","FBN1","ROR1","PLCH1","CACNB4","CPLX4",
                          "CPAMD8","UTRN","PLD1","SDK1","COLEC12","HSD17B2")

data_regen_scRNA <- AddModuleScore(data_regen_scRNA,
                                   features = list(mecom_targeted_genes),
                                   name="mecom_targeted_genes")

#######Supplemental Figure S12E1
W2<-FeaturePlot(data_regen_scRNA, reduction = "umap",
                features = "mecom_targeted_genes1", label = FALSE, repel = TRUE) +
  scale_colour_gradient2(high = "#DC0000FF",mid = "#EDE361",low = "#69B9DA",midpoint =2)
W2


data_matrix <-data_regen_scRNA@assays$SCT@data 
#class(data_matrix)
#head(data.frame(data_matrix))
data_MECOM <- data_matrix[which(rownames(data_matrix) == "MECOM"),]
class(data_MECOM)
data_expression <- data_MECOM
celltypes <- as.character(data_regen_scRNA$final_annotated_0807_2)
data_exp <- as.data.frame(cbind(celltypes,data_expression))
time_points <- as.character(data_regen_scRNA$treatment)
data_exp$timepoints <- time_points

mecom_targeted_genes_score <- data_regen_scRNA$mecom_targeted_genes1
data_exp$mecom_targeted_genes_score  <- mecom_targeted_genes_score 
names(data_exp) <- c("celltypes_10","expression","timepoints","mecom_targeted_genes_score")

table(data_exp$timepoints)


data_exp$timepoints<-factor(data_exp$timepoints,levels = c("CMZ-0","CMZ-05","CMZ-10","CMZ-20","CMZ-40"))

#######Figure 6E
p = ggplot(data_exp, aes(x=timepoints, y=as.numeric(expression))) + 
  stat_boxplot(geom = "errorbar",width=0.2, size=0.5,position=position_dodge(0.6),color= "black")+
  theme_classic()+
  geom_boxplot(position = position_dodge(0.6),
               size = 0.5,
               width = 0.8,
               fill = c("#af2157","#BD8098","#C0BFDF","#CEE1EB","#67ADB7"),
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
  ylab("MECOM expression")+
  theme(axis.title = element_text(size=18),
        axis.text = element_text(size=14),
        axis.text.x = element_text(angle = 35,
                                   hjust=1,
                                   size = 10))
p



####Supplemental Figure S12E2
p1 = ggplot(data_exp, aes(x=timepoints, y=as.numeric(mecom_targeted_genes_score))) + 
  stat_boxplot(geom = "errorbar",width=0.2, size=0.5,position=position_dodge(0.6),color= "black")+
  theme_classic()+
  geom_boxplot(position = position_dodge(0.6),
               size = 0.5,
               width = 0.8,
               fill = c("#af2157","#BD8098","#C0BFDF","#CEE1EB","#67ADB7"),
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
  ylab("Mecom targeted genes score")+
  theme(axis.title = element_text(size=18),
        axis.text = element_text(size=14),
        axis.text.x = element_text(angle = 35,
                                   hjust=1,
                                   size = 10))
p1
