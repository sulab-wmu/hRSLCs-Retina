##############--------------------------------Heatmap---plot-------------------Figure 3J
library(ggplot2)
jacd_mat_min <- readr::read_csv("jacad_list_0816.csv")

ggplot(jacd_mat_min,aes(x=X2,y=X1,fill= value))+
  geom_tile()+
  scale_fill_gradient2(low="#67ADB7",high="#af2157")+ggtitle('min')+
  theme(axis.text.x = element_text(angle=90,
                                   hjust = 1,
                                   vjust = 1))