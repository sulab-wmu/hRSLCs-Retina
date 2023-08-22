##Supplemental Figure S10E
library(Seurat)

##organoids singlecell data
single_data<-readRDS("C:/Users/Lenovo/Desktop/regeneration/Final_Reg_CMZ_0807.rds")

##import the result from organoids_Stem_ID.R
load("C:/Users/Lenovo/Desktop/regeneration/Stem_ID/entropy.Rdata")
single_data<-subset(single_data,idents = c("hRSLCs"))
single_data$entropy<-entropy

entropy <- as.numeric(single_data$entropy)
time_points <- as.character(single_data$treatment)
data_entropy <- as.data.frame(cbind(time_points,entropy))
colnames(data_entropy) <- c("time_points","entropy")

head(data_entropy)
data_entropy$time_points<-factor(data_entropy$time_points,levels = c("CMZ-0","CMZ-05", "CMZ-10", "CMZ-20", "CMZ-40"))

library(ggplot2)
#boxplot3
p = ggplot(data_entropy, aes(x=time_points, y=as.numeric(entropy))) + 
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
  ylab("Stemness score")+
  theme(axis.title = element_text(size=18),
        axis.text = element_text(size=14),
        axis.text.x = element_text(angle = 35,
                                   hjust=1,
                                   size = 10))
p
