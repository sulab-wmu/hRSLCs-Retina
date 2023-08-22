##########Proportion of four cell types across five time points
#Subset of distinct time points------------------------------------------------Figure 5E

##organoids singlecell data
data_regen_scRNA<-readRDS("C:/Users/Lenovo/Desktop/regeneration/Final_Reg_CMZ_0807.rds")

Idents(data_regen_scRNA) <- data_regen_scRNA$treatment
table(Idents(data_regen_scRNA))
  
##CMZ-0
subset_regen_CMZ0 <- subset(data_regen_scRNA,idents = c("CMZ-0"))
pal <- c("#85C89C","#67ADB7","#36600E","#6A8473","#C0BFDF",
         "#E77A77","#7E6148FF","#ECBA84","#CA8C74","#A13B46")
  
Idents(subset_regen_CMZ0) <- subset_regen_CMZ0$final_annotated_0807_2
DimPlot(subset_regen_CMZ0, reduction = "umap", label = F,  cols= pal, pt.size = 0.2, repel = T)
table1<-table(subset_regen_CMZ0$final_annotated_0807_2)
table2<- as.data.frame(table1)
table2$freq2 <- table2$Freq/length(subset_regen_CMZ0$final_annotated_0807_2)
  
  
#CMZ-05
subset_regen_CMZ5 <- subset(data_regen_scRNA,idents = c("CMZ-05"))
pal <- c("#85C89C","#67ADB7","#36600E","#6A8473","#C0BFDF",
         "#E77A77","#7E6148FF","#ECBA84","#CA8C74","#A13B46")
  
Idents(subset_regen_CMZ5) <- subset_regen_CMZ5$final_annotated_0807_2
DimPlot(subset_regen_CMZ5, reduction = "umap", label = F,  cols= pal, pt.size = 0.2, repel = T)
table(subset_regen_CMZ5$final_annotated_0807_2)
table12<-table(subset_regen_CMZ5$final_annotated_0807_2)
table22<- as.data.frame(table12)
table22$freq2 <- table22$Freq/length(subset_regen_CMZ5$final_annotated_0807_2)
table22
  
#CMZ-10 
subset_regen_CMZ10 <- subset(data_regen_scRNA,idents = c("CMZ-10"))
pal <- c("#85C89C","#67ADB7","#36600E","#6A8473","#C0BFDF",
         "#E77A77","#7E6148FF","#ECBA84","#CA8C74","#A13B46")
  
Idents(subset_regen_CMZ10) <- subset_regen_CMZ10$final_annotated_0807_2
DimPlot(subset_regen_CMZ10, reduction = "umap", label = F,  cols= pal, pt.size = 0.2, repel = T)
table(subset_regen_CMZ10$final_annotated_0807_2)
table13<-table(subset_regen_CMZ10$final_annotated_0807_2)
table23<- as.data.frame(table13)
table23$freq2 <- table23$Freq/length(subset_regen_CMZ10$final_annotated_0807_2)
table23
  
    
#CMZ-20 
subset_regen_CMZ20 <- subset(data_regen_scRNA,idents = c("CMZ-20"))
pal <- c("#85C89C","#67ADB7","#36600E","#6A8473","#C0BFDF",
         "#E77A77","#7E6148FF","#ECBA84","#CA8C74","#A13B46")
  
Idents(subset_regen_CMZ20) <- subset_regen_CMZ20$final_annotated_0807_2
DimPlot(subset_regen_CMZ20, reduction = "umap", label = F,  cols= pal, pt.size = 0.2, repel = T)
table(subset_regen_CMZ20$final_annotated_0807_2)
table14<-table(subset_regen_CMZ20$final_annotated_0807_2)
table24<- as.data.frame(table14)
table24$freq2 <- table24$Freq/length(subset_regen_CMZ20$final_annotated_0807_2)
table24
  
  
#CMZ-40 
subset_regen_CMZ40 <- subset(data_regen_scRNA,idents = c("CMZ-40"))
pal <- c("#85C89C","#67ADB7","#36600E","#6A8473","#C0BFDF",
         "#E77A77","#7E6148FF","#ECBA84","#CA8C74","#A13B46")
  
Idents(subset_regen_CMZ40) <- subset_regen_CMZ40$final_annotated_0807_2
DimPlot(subset_regen_CMZ40, reduction = "umap", label = F,  cols= pal, pt.size = 0.2, repel = T)
table(subset_regen_CMZ40$final_annotated_0807_2)
table15<-table(subset_regen_CMZ40$final_annotated_0807_2)
table25<- as.data.frame(table15)
table25$freq2 <- table25$Freq/length(subset_regen_CMZ40$final_annotated_0807_2)
table25
  
###Construct a data.frame for all time points
all_prop <- table2
names(all_prop)<-c("Celltype","CMZ0_count","CMZ0_prop")
all_prop$CMZ5_count <- table22$Freq
all_prop$CMZ5_prop <- table22$freq2
all_prop$CMZ10_count <- table23$Freq
all_prop$CMZ10_prop <- table23$freq2
all_prop$CMZ20_count <- table24$Freq
all_prop$CMZ20_prop <- table24$freq2
all_prop$CMZ40_count <- table25$Freq
all_prop$CMZ40_prop <- table25$freq2
###Final_version
all_prop


###scale cell proportion---------------------------------------------------------------Figure 5E
all_prop_melt_sub2 <-all_prop[,which(colnames(all_prop) %in% c("Celltype","CMZ0_prop","CMZ5_prop","CMZ10_prop","CMZ20_prop","CMZ40_prop"))]

data <- all_prop[,1]

all_prop_melt_sub3 <- all_prop_melt_sub2[,-1]
rownames(all_prop_melt_sub3) <- data
all_prop_melt_sub4 <- all_prop_melt_sub3/all_prop_melt_sub3[,1] #Using the first time point as scale
names(all_prop_melt_sub4)<-c("CMZ0","CMZ5","CMZ10","CMZ20","CMZ40")
all_prop_melt_sub4$Celltype <- data

#####---------------
library(reshape2)
all_prop_melt_new <- melt(all_prop_melt_sub4,id="Celltype")
names(all_prop_melt_new) <- c("Celltype","timepoints","proportion")
all_prop_melt_new_sub<-all_prop_melt_new[which(all_prop_melt_new$Celltype %in% c("hRSLCs","RPCs","PC_precursors","PCs")),]

#ggplot-------------------------------
library(ggplot2)
p <- ggplot(all_prop_melt_new_sub, aes(x=timepoints,y=proportion,group = Celltype))+
  geom_line(aes(color=Celltype))+
  scale_color_manual(values = c("#af2157","#BD8098","#C0BFDF","#67ADB7"))+
  theme_classic()+
  geom_point(aes(color=Celltype),size=5)+
  ylab("Fold Change")+
  xlab("Time points")+
  theme(legend.title = element_text(size = 14),
        legend.text = element_text(size=12),
        plot.title = element_text(face="bold",size = 20),
        axis.title.x = element_text(face="bold",size = 12),
        axis.title.y = element_text(face="bold",size = 12),
        axis.text.y = element_text(size=10),
        axis.text.x = element_text(angle = 50,
                                   hjust=1,
                                   size = 10,
                                   color = "black"))
p      














  