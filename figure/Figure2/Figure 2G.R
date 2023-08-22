###############-------------------------Heatmap plot------------------------for Figure 2G.
library(pheatmap)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
load("DEseq_0813_2.RData")
colnames(txi.rsem2_sub)
txi.rsem2_sub2<- as.data.frame(txi.rsem2_sub)
colnames(txi.rsem2_sub2) <- c("A1","A2","A3","B1","B2","B3")
rownames(samples) <-  c("A1","A2","A3","B1","B2","B3")

samples$time_points <- c("Day21","Day21","Day21","Day0","Day0","Day0")

pheatmap(txi.rsem2_sub2,scale = 'row',
         annotation_col = samples[6],
         annotation_row = de_list[3],
         cluster_rows = F,
         cluster_cols = T,
         show_colnames = F,
         fontsize = 6,
         border_color = "white")