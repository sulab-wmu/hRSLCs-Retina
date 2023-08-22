######------------------------------------------------------Figure 7B
#fluorescence data
data_fluo <- read.csv("Fluorescence data.csv")

data_fluo <-data_fluo[,3:7]

###########Raw plot_for cell proportion
library(reshape2)
library(ggplot2)
all_fluo_melt <- melt(data_fluo,id="ID")
names(all_fluo_melt) <- c("Celltype","timepoints","values")

#ggplot
p <- ggplot(all_fluo_melt, aes(x=timepoints,y=values,group = Celltype))+
  geom_line(aes(color=Celltype))+
  scale_color_manual(values = c("#36600E","#36600E","#36600E","#36600E","#36600E","#36600E","#36600E",
                                "#36600E","#A13B46","#A13B46","#A13B46","#A13B46",
                                "#A13B46","#A13B46","#A13B46","#A13B46","#A13B46","#A13B46","#A13B46"))+
  theme_classic()+
  geom_point(aes(color=Celltype),size=3)+
  ylab("Values")+
  xlab("Time points")+
  theme(legend.title = element_text(size = 14),
        legend.text = element_text(face="italic",size=12),
        plot.title = element_text(face="bold",size = 20),
        axis.title.x = element_text(face="bold",size = 12),
        axis.title.y = element_text(face="bold",size = 12),
        axis.text.y = element_text(size=10),
        axis.text.x = element_text(angle = 35,
                                   hjust=1,
                                   size = 10,
                                   color = "black"))
p      