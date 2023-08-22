##Figure 5G
library(Seurat)

##organoids singlecell data
single_data<-readRDS("C:/Users/Lenovo/Desktop/regeneration/Final_Reg_CMZ_0807.rds")

##import the result from organoids_CytoTRACE.R
CytoTRACE<-readRDS("C:/Users/Lenovo/Desktop/regeneration/CytoTRACE/CytoTRACE.rds")
CytoTRACE<-as.data.frame(CytoTRACE$CytoTRACE)
single_data$CytoTRACE<-CytoTRACE

metadata<-single_data@meta.data
metadata<-metadata[,c("treatment","CytoTRACE")]

metadata$treatment<-factor(metadata$treatment,levels = c("CMZ-0","CMZ-05", "CMZ-10", "CMZ-20", "CMZ-40"))

library(ggplot2)
p = ggplot(metadata, aes(x=treatment, y=as.numeric(CytoTRACE))) + 
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
  ylab("Stemness")+
  theme(axis.title = element_text(size=18),
        axis.text = element_text(size=14),
        axis.text.x = element_text(angle = 35,
                                   hjust=1,
                                   size = 10))
p
