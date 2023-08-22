#-------------------------------Supplemental Figure S12C
library(ggplot2)
a <- ggplot(data_node,aes(x=UMAP_1,y=UMAP_2,color=hRSLCs))+
  geom_point(size=2.5)+scale_color_gradientn(values = seq(0,1,0.2),
                                             colors = c("grey","yellow","orange","red"))

a+theme_bw()+theme(panel.border=element_blank(),panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),axis.line = element_line(color = "black"))



a <- ggplot(data_node,aes(x=UMAP_1,y=UMAP_2,color=RPCs))+
  geom_point(size=2.5)+scale_color_gradientn(values = seq(0,1,0.2),
                                             colors = c("grey","yellow","orange","red"))

a+theme_bw()+theme(panel.border=element_blank(),panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),axis.line = element_line(color = "black"))


a <- ggplot(data_node,aes(x=UMAP_1,y=UMAP_2,color=PC_precursors))+
  geom_point(size=2.5)+scale_color_gradientn(values = seq(0,1,0.2),
                                             colors = c("grey","yellow","orange","red"))

a+theme_bw()+theme(panel.border=element_blank(),panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),axis.line = element_line(color = "black"))



a <- ggplot(data_node,aes(x=UMAP_1,y=UMAP_2,color=PCs))+
  geom_point(size=2.5)+scale_color_gradientn(values = seq(0,1,0.2),
                                             colors = c("grey","yellow","orange","red"))

a+theme_bw()+theme(panel.border=element_blank(),panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),axis.line = element_line(color = "black"))