##Figure 4G
##load the result from bulk_PCA_analysis.R
load("PCA_result.Rdata")

##pca plot
library(ggplot2)
library(ggrepel)
library(cowplot)
pdf("PCA.pdf")
ggplot(pca_sample, aes(x = Dim.1, y = Dim.2, color = Group, shape = Group),arrow = arrow(length=unit(0.2, "cm"))) + 
  geom_point(size=3) + 
  scale_color_manual(values = c("#E77A77","#67ADB7" )) +  
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.key = element_rect(fill = 'transparent')) +   
  labs(x =  paste('PCA1(', pca_eig1, '%)'), y = paste('PCA2(', pca_eig2, '%)'))+  
  geom_text_repel(aes(label = Sample), size = 4, show.legend = FALSE, 
                  box.padding = unit(0.2, 'lines'))+
  geom_hline(aes(yintercept=0),color="grey",linetype="dashed")+  
  geom_vline(aes(xintercept=0),color="grey",linetype="dashed")+  
  geom_smooth(se = F, method = 'loess',span=1) 
dev.off()




