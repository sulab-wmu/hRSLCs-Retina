################Figure 5H
markers_subset_regen_hRSLCs_CMZ5_vs_0 <- read.csv("markers_subset_regen_hRSLCs_CMZ5_vs_0.csv")

head(markers_subset_regen_hRSLCs_CMZ5_vs_0)

DEG_genes2<-markers_subset_regen_hRSLCs_CMZ5_vs_0
plot(DEG_genes2$avg_log2FC,-log10(DEG_genes2$p_val_adj),
     xlim=c(-2, 2), ylim=c(0,45),
     xlab="Log2 (fold change)", 
     ylab="-Log10(FDR)")
r10 <- DEG_genes2[which(DEG_genes2$avg_log2FC >=0.25  & DEG_genes2$p_val_adj<=0.05),] 
x <- r10$avg_log2FC
y <- -log10(r10$p_val_adj)
points(x,y,pch = 21,cex = 1.25, lwd = 1.3, col = "black", bg =  "#BD8098")

r11 <- DEG_genes2[which(DEG_genes2$avg_log2FC<=-0.25 & DEG_genes2$p_val_adj<=0.05),] 
x <- r11$avg_log2FC
y <- -log10(r11$p_val_adj)
points(x,y,pch = 21,cex = 1.25, lwd = 1.3, col = "black", bg =  "#0F4450")


abline(v=0.25,col="red",lty="longdash")
abline(v=-0.25,col="red",lty="longdash")
abline(h=-log(0.05,10),col="red",lty="longdash")