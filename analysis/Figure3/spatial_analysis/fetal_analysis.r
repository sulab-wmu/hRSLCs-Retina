library(Seurat)
setwd("/share/pub/dengcy/STanalysis/fetusdata/analysis/")
setwd("D:\\OneDrive\\Retina_spatial\\fetal_analysis")
fetal_data<-readRDS("W19-fetal_all_analysis_bin50.rds")
#setwd("/share/pub/dengcy/STanalysis/fetusdata/analysis/")
fetal_data<-RunTSNE(object = fetal_data,assay = "SCT", reduction = 'pca', dims = 1:30,check_duplicates = FALSE)

pdf("fetal_data_umap_cluster.pdf")
DimPlot(fetal_data, reduction = "umap", label = TRUE)
dev.off()
#obj2 <- FindNeighbors(object=fetal_data, reduction = reduction, dims = 1:50)
fetal_data <- FindClusters(object=fetal_data, resolution = 4)

# 创建一个简单的 SVG 图像
svg("W19-fetal_dimplot_cluster4.svg", width =12, height =9)
SpatialDimPlot(fetal_data,label = FALSE, label.size = 3,stroke=0)
dev.off()ß
saveRDS(fetal_data1,file="sub_W19-fetal_scdata.rds")
fetal_data<-readRDS("sub_W19-fetal_scdata.rds")
saveRDS(fetal_data,file="W19-fetal_scdata.rds")

############
#可视化一些特殊基因
############
fetal_data<-readRDS("W19-fetal_scdata.rds")

pdf("BAZ2B_spatial_W19-fetal_scdata.pdf")
SpatialFeaturePlot(fetal_data, features = c('BAZ2B'))
dev.off()
pdf("DDC_spatial_W19-fetal_scdata.pdf")
SpatialFeaturePlot(fetal_data, features = c('DDC'))
dev.off()
pdf("FKBP5_spatial_W19-fetal_scdata.pdf")
SpatialFeaturePlot(fetal_data, features = c('FKBP5'))
dev.off()

pdf("MECOM_spatial_W19-fetal_scdata.pdf",width=8,height=8)
SpatialFeaturePlot(fetal_data, features = c('MECOM'))
dev.off()

fetal_data<-readRDS("W19-fetal_scdata.rds")

svg("W19-fetal_dimplot_fetal_cheminggenes.svg",width=50,height=40)
SpatialFeaturePlot(fetal_data, features =c('CRYAB', 'CRYBB1','CRYBB3', 'CRYBB2' ,'CRYGS')
, stroke = 0)
dev.off()

#######
#计算差异markers
#############
de.markers <- FindAllMarkers(fetal_data)
write.csv(de.markers,file="cluster4_demarkers.csv")

#################
#可视化细胞类型
###################
anno<-read.csv("annotation_2.csv")
fetal_data$annotation<-""
anno<-unique(anno[,1:2])
for(i in anno$CLUSTER){
    fetal_data$annotation[fetal_data$seurat_clusters==i]<-anno$CELL.TYPES[anno$CLUSTER==i]
}
svg("W19-fetal_dimplot_annotation.svg",width=10,height=10)
SpatialDimPlot(fetal_data,group.by = "annotation",label = FALSE, label.size = 3,stroke=0)
dev.off()
colors<-c('#1f77b4', '#ff7f0e', '#279e68', '#bec1d4', "#36600E", "#9569AB", '#7d87b9', "#CA8C74",  '#8e063b', '#8c6d31', '#d62728', '#f7b6d2', '#dbdb8d', '#c49c94', '#c5b0d5', '#ff9896', '#98df8a', '#ffbb78', '#aec7e8', '#17becf','#b5bd61', "#E77A77", '#aa40fc', '#9edae5','#023fa5')
library(ggplot2)
library(dplyr)
Idents(fetal_data)<-fetal_data$annotation
names(colors) <- Idents(fetal_data) %>% levels()
svg("W19-fetal_dimplot_annotation.svg",width=10,height=10)
SpatialDimPlot(fetal_data,group.by = "annotation",label = FALSE, label.size = 3,stroke=0,cols =colors)
dev.off()

saveRDS(fetal_data,file="W19-fetal_scdata.rds")

fetal_data$cluster_52<-"other"
fetal_data$cluster_52[fetal_data$seurat_clusters==52]<-"52"

svg("W19-fetal_dimplot_52.svg",width=10,height=10)
SpatialDimPlot(fetal_data,group.by = "cluster_52",label = FALSE, label.size = 3,stroke=0)+scale_colour_manual(values=colors)
dev.off()
############
#sub spatial plot
#############
summary(fetal_data@meta.data$x)
fetal_data<-readRDS("W19-fetal_scdata.rds")
fetal_data1<-subset(x = fetal_data, subset = x < 5000)
fetal_data1<-subset(x = fetal_data1, subset = y < 11000)
fetal_data1<-subset(x = fetal_data1, subset = y > 3000)
#fetal_data1<-subset(x = fetal_data1, subset = x > 3200)
svg("Sub2_W19-fetal_dimplot_cluster4.svg",width=2,height=2)
SpatialDimPlot(fetal_data1,label = FALSE, label.size = 3,stroke=0)
dev.off()

saveRDS(fetal_data1,file="sub_W19-fetal_scdata.rds")
saveRDS(fetal_data,file="W19-fetal_scdata.rds")
##
fetal_data1<-subset(x = fetal_data, subset = seurat_clusters %in% 43)

#fetal_data1_1<-subset(x = fetal_data1, subset = x < 5000)
#fetal_data1_1<-subset(x = fetal_data1_1, subset = y < 11000)
#fetal_data1_1<-subset(x = fetal_data1_1, subset = y > 3000)
svg("Sub2_W19-fetal_dimplot_cluster4.svg",width=2,height=2)
SpatialDimPlot(fetal_data1_1,label = FALSE, label.size = 3,stroke=0)
dev.off()

#重新聚类
fetal_data1<-RunPCA(fetal_data1, assay = "SCT")
fetal_data1<-FindNeighbors(object=fetal_data1, reduction = 'pca', dims = 1:30)
fetal_data1<- RunUMAP(object = fetal_data1, assay ="SCT", reduction = 'pca', dims = 1:30,check_duplicates = FALSE)
fetal_data1<-RunTSNE(object = fetal_data1,assay = "SCT", reduction = 'pca', dims = 1:30,check_duplicates = FALSE)
fetal_data1 <- FindClusters(object=fetal_data1, resolution = 1)
library(ggplot2)
pdf("umap_clusterCMZ.pdf",width=5,height=5)
DimPlot(fetal_data1, reduction = "umap", label = TRUE)
dev.off()

svg("CMZ_W19-fetal_dimplot_cluster1.svg")
SpatialDimPlot(fetal_data1,label = FALSE, label.size = 3,stroke=0)
dev.off()

####将cmz区域进一步注释
#0归为LCs，1归为RPE，2归为hNRSCs
meta.data=fetal_data@meta.data
meta.data[colnames(fetal_data1)[fetal_data1$seurat_clusters==0],'annotation']<-'LCs'
meta.data[colnames(fetal_data1)[fetal_data1$seurat_clusters==1],'annotation']<-'RPE.stem-like cells'
meta.data[colnames(fetal_data1)[fetal_data1$seurat_clusters==2],'annotation']<-'hNRSCs'

table(meta.data$annotation)
fetal_data$annotation<-meta.data$annotation

colors<-c('#1f77b4', '#ff7f0e', '#279e68', '#B4B4B3', '#aec7e8', "#9569AB", '#7d87b9', "#CA8C74",  '#8e063b',  '#c5b0d5','#023fa5', '#f7b6d2', '#E9B824', '#c49c94', '#d62728', '#ff9896', '#98df8a', '#ffbb78', "#36600E", '#17becf','#b5bd61', "#E77A77", '#aa40fc', '#9edae5','#bec1d4','#67ADB7')
Idents(fetal_data)<-fetal_data$annotation
names(colors) <- Idents(fetal_data) %>% levels()
colors['RPE.stem-like cells']<-'#CD5C08'
colors['hNRSCs']<-'#F5E8B7'
colors['LCs']<-'#6A9C89'
colors['Unknown']<-'#B4B4B3'
svg("W19-fetal_dimplot_annotation.svg",width=10,height=10)
SpatialDimPlot(fetal_data,group.by = "annotation",label = FALSE, label.size = 3,stroke=0,cols =colors)
dev.off()


pdf("sub_fetal_data_MECOM.pdf")
FeaturePlot(fetal_data1_1,feature='MECOM')
dev.off()

###########
#区域伪时序分析
##############
fetal_data1_1<-readRDS("sub_cmzcluster_W19-fetal_scdata.rds")
library(slingshot, quietly = FALSE)
sce <- SingleCellExperiment(assays = List(counts = fetal_data1_1@assays$SCT@counts))

assays(sce)$norm <- fetal_data1_1@assays$SCT@data

#pca <- prcomp(t(log1p(assays(sce)$norm)), scale. = FALSE)
rd1 <- fetal_data1_1@reductions$pca[,1:2]

#plot(rd1, col = rgb(0,0,0,.5), pch=16, asp = 1)

rd2 <- as.data.frame(Embeddings(fetal_data1_1, reduction = "tsne"))
colnames(rd2) <- c("tsne_1", "tsne_2")

#plot(rd2, col = rgb(0,0,0,.5), pch=16, asp = 1)

reducedDims(sce) <- SimpleList(PCA = rd2, UMAP = rd2)
#library(mclust, quietly = TRUE)
#cl1 <- Mclust(rd1)$classification
colData(sce)$seurat_clusters<- fetal_data1_1$seurat_clusters
#library(RColorBrewer)
#plot(rd1, col = brewer.pal(9,"Set1")[cl1], pch=16, asp = 1)
#cl2 <- kmeans(rd1, centers = 4)$cluster
#colData(sce)$kmeans <- cl2

sce <- slingshot(sce, clusterLabels = 'seurat_clusters', reducedDim = 'PCA')

library(grDevices)
library(RColorBrewer)
colors <- colorRampPalette(c('#FBFACD','#DEBACE','#BA94D1','#7F669D'))(20)
plotcol <- colors[cut(sce$slingPseudotime_1, breaks=20)]
pdf("Slingshot_sub_cmz_data.pdf")
plot(reducedDims(sce)$PCA, col = plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(sce), lwd=2, col='black')
plot(reducedDims(sce)$PCA, col = brewer.pal(9,'Set1')[sce$GMM], pch=16, asp = 1)
lines(SlingshotDataSet(sce), lwd=2, type = 'lineages', col = 'black')
dev.off()
save(sce,file='D:\\OneDrive\\Retina_spatial\\fetal_analysis\\monocle_analysi\\cma_subdata_slingPseudotime.RData')
###做伪时序boxplot

gg<-data.frame(cluster=fetal_data1_1$seurat_clusters,Pseudotime=sce$slingPseudotime_1)
# 加载必要的包
library(ggplot2)
library(dplyr)
library(ggpubr)
# 绘制boxplot
pdf("cmz_Pseudotime_boxplot.pdf")
ggplot(gg, aes(x = cluster, y = Pseudotime,fill = cluster, color = cluster)) +
  geom_boxplot() +
  theme_minimal() +
  labs(x = "Cluster", y = "Pseudotime") +
  scale_y_continuous(breaks = seq(0, max(gg$Pseudotime), by = 1)) +
 stat_compare_means( method = "wilcox.test", label = "p.signif")

dev.off()


#marker genes
de.markers <- FindAllMarkers(fetal_data1_1,assay ='SCT',slot = "data",only.pos = TRUE,logfc.threshold = 0.2,return.thresh = 0.05)
#de.markers[de.markers$cluster==2,'gene'] == 'MECOM'

write.csv(de.markers,file="sub_cmz_cluster_demarkers.csv")
pdf("dotplot_sub_spatial_CMZ_data_marker.pdf",height=3,width=15)
DotPlot(fetal_data1_1, features =de.markers$gene, cols = c("#FFE4D6", "#B0578D")) + RotatedAxis()
dev.off()

saveRDS(fetal_data1_1,file="sub_cmzcluster_W19-fetal_scdata.rds")


##############
#gene marker plot
############
fetal_data<-readRDS("W19-fetal_scdata.rds")
genemarker<-read.csv('singlecell_fetal_genemarkers.csv')
gene_list<-lapply(unique(genemarker$CELL.TYPLE),function(x) genemarker$GENES[genemarker$CELL.TYPLE==x])
names(gene_list)<-unique(genemarker$CELL.TYPLE)
fetal_data<-AddModuleScore( object=fetal_data,  features=gene_list)

new_colnames <- c("orig.ident","nCount_Spatial","nFeature_Spatial", "cell","x", "y","nCount_SCT","nFeature_SCT","SCT_snn_res.0.8","seurat_clusters","SCT_snn_res.4",unique(genemarker$CELL.TYPLE))

# 将新的列名向量赋值给 meta.data 列名
colnames(fetal_data@meta.data) <- new_colnames

svg("W19-fetal_dimplot_fetal_genemarker.svg",width=50,height=40)
SpatialFeaturePlot(fetal_data, features =unique(genemarker$CELL.TYPLE), stroke = 0)
dev.off()


fetal_data1_1<-AddModuleScore( object=fetal_data1_1,  features=gene_list[c("RPR_progenitors","CMZ_RSCs")],nbin = 10)

colnames(fetal_data1_1@meta.data)[17:18]<-c("RPR_progenitors","CMZ_RSCs")
svg("sub_W19-fetal_dimplot_fetal_genemarker.svg",width=4,height=4)
SpatialFeaturePlot(fetal_data1_1, features =c("RPR_progenitors","CMZ_RSCs"), stroke = 0)
dev.off()

#################
#整体数据的伪时间分析
#用python环境
library(Seurat)
library(feather)
setwd("/share/pub/dengcy/STanalysis/fetusdata/analysis/")
fetal_data<-readRDS("W19-fetal_scdata.rds")
x<-as.data.frame(fetal_data@assays$SCT@counts)
write_feather(x, "sct_matrix.feather")
write.csv(data.frame(rownames(fetal_data@assays$SCT@counts)), file="genenames.csv")

write.csv(fetal_data@meta.data, file="metadata.csv")
intersect(rownames(fetal_data),'MECOM')
######################
#读取伪时间数据
############
setwd("/share/pub/dengcy/STanalysis/fetusdata/analysis/")
fetal_data<-readRDS("W19-fetal_scdata.rds")
obs<-read.csv('D:\\OneDrive\\Retina_spatial\\fetal_analysis\\monocle_analysi\\obs_dpt_pseudotime_louvain.csv')
fetal_data$dpt_pseudotime<-obs$dpt_pseudotime
fetal_data$dpt_pseudotime[fetal_data$dpt_pseudotime>0.2] <- 0.2
pdf("D:\\OneDrive\\Retina_spatial\\fetal_analysis\\monocle_analysi\\dpt_pseudotime_featureplot.pdf")
FeaturePlot(fetal_data, features = 'dpt_pseudotime', cols = c("#E55604","#EBE4D1" ,"#26577C"))
dev.off()

#######################
#R 环境的伪时间分析
######################
#slingshot
library(Seurat)

#setwd("/share/pub/dengcy/STanalysis/fetusdata/analysis/")
setwd("D:\\OneDrive\\Retina_spatial\\fetal_analysis")
fetal_data<-readRDS("W19-fetal_scdata.rds")

sub_fetal_data<-subset(fetal_data,subset= annotation %in% c("PC","RPC","RSC/RPE"))
#hRSLCs到RPC到PC
sub_fetal_data<-RunPCA(sub_fetal_data, assay = "SCT")
sub_fetal_data<-FindNeighbors(object=sub_fetal_data, reduction = 'pca', dims = 1:30)
sub_fetal_data<- RunUMAP(object = sub_fetal_data, assay ="SCT", reduction = 'pca', dims = 1:30,check_duplicates = FALSE)
sub_fetal_data<-RunTSNE(object = sub_fetal_data,assay = "SCT", reduction = 'pca', dims = 1:30,check_duplicates = FALSE)
sub_fetal_data <- FindClusters(object=sub_fetal_data, resolution = 1)
library(slingshot, quietly = FALSE)
sce <- SingleCellExperiment(assays = List(counts = sub_fetal_data@assays$SCT@counts))

assays(sce)$norm <- sub_fetal_data@assays$SCT@data

#pca <- prcomp(t(log1p(assays(sce)$norm)), scale. = FALSE)
rd1 <- sub_fetal_data@reductions$pca[,1:2]

plot(rd1, col = rgb(0,0,0,.5), pch=16, asp = 1)

rd2 <- as.data.frame(Embeddings(sub_fetal_data, reduction = "umap"))
colnames(rd2) <- c("UMAP_1", "UMAP_2")

plot(rd2, col = rgb(0,0,0,.5), pch=16, asp = 1)

reducedDims(sce) <- SimpleList(PCA = rd2, UMAP = rd2)
#library(mclust, quietly = TRUE)
#cl1 <- Mclust(rd1)$classification
colData(sce)$annotation <- sub_fetal_data$annotation
#library(RColorBrewer)
#plot(rd1, col = brewer.pal(9,"Set1")[cl1], pch=16, asp = 1)
#cl2 <- kmeans(rd1, centers = 4)$cluster
#colData(sce)$kmeans <- cl2

sce <- slingshot(sce, clusterLabels = 'annotation', reducedDim = 'PCA')

library(grDevices)
library(RColorBrewer)
colors <- colorRampPalette(c('#FBFACD','#DEBACE','#BA94D1','#7F669D'))(20)
plotcol <- colors[cut(sce$slingPseudotime_1, breaks=20)]
pdf("Slingshot_sub_fetal_data.pdf")
plot(reducedDims(sce)$PCA, col = plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(sce), lwd=2, col='black')
plot(reducedDims(sce)$PCA, col = brewer.pal(9,'Set1')[sce$GMM], pch=16, asp = 1)
lines(SlingshotDataSet(sce), lwd=2, type = 'lineages', col = 'black')
dev.off()
save(sce,file='D:\\OneDrive\\Retina_spatial\\fetal_analysis\\monocle_analysi\\sce_slingPseudotime.RData')

sub_fetal_data$slingPseudotime<-sce$slingPseudotime_1

library(ggplot2)
pdf("umap_sub_fetal_data_slingPseudotime.pdf",width=5,height=5)
FeaturePlot(sub_fetal_data, features = "slingPseudotime", min.cutoff = 1, max.cutoff = 3)
dev.off()

pdf("umap_sub_fetal_data_annotation.pdf",width=5,height=5)
DimPlot(sub_fetal_data, group.by='annotation')
dev.off()
sub_fetal_data <- FindClusters(object=sub_fetal_data, resolution = 0.5)
pdf("umap_sub_fetal_data_clusters.pdf",width=5,height=5)
DimPlot(sub_fetal_data)
dev.off()


de.markers <- FindAllMarkers(sub_fetal_data,assay ='SCT',slot = "data",only.pos = TRUE,logfc.threshold = 1,return.thresh = 0.05)
#de.markers[de.markers$cluster==2,'gene'] == 'MECOM'

write.csv(de.markers,file="sub_fetal_data_demarkers.csv")
pdf("dotplot_sub_fetal_data_marker.pdf",height=3,width=13)
DotPlot(sub_fetal_data, features =de.markers$gene, cols = c("#FFE4D6", "#B0578D")) + RotatedAxis()
dev.off()

saveRDS(sub_fetal_data,file="sub_fetal_data_pc_rsc.rds")
fetal_data1_1<-readRDS("sub_fetal_data.rds")

saveRDS(fetal_data,file="W19-fetal_scdata.rds")

#################################
#将cmz区域的重新伪时序分析结果整合到整体的结果中
#############################


pdf("W19-fetal_dimplot_annotation.pdf",width=17,height=15)
SpatialDimPlot(fetal_data,group.by = "annotation",label = FALSE, label.size = 3,stroke=0,cols =colors,label.box = FALSE,pt.size.factor = 1.3)
dev.off()

svg("sub_W19-fetal_spatial_slingPseudotime.svg",width=13,height=13)
SpatialFeaturePlot(sub_fetal_data, features = "slingPseudotime")
dev.off()

#gg<-sub_fetal_data@meta.data


de.markers <- FindAllMarkers(fetal_data,assay ='SCT',slot = "data",only.pos = TRUE,logfc.threshold = 0.07,return.thresh = 0.05)
write.csv(de.markers,file="fetal_data_annotation_demarkers0.1.csv")
de.markers<-de.markers[!duplicated(de.markers$gene),]
genes_list<-c()
for(i in as.vector(unique(de.markers$cluster))){
  genes<-de.markers$gene[de.markers$cluster==i]
  genes_list<-c(genes_list,genes[1:5])
}

pdf("dotplot_fetal_annotation_marker.pdf",height=8,width=24)
DotPlot(fetal_data, features =unique(genes_list), cols = c("#FFE4D6", "#B0578D")) + RotatedAxis()+ theme(legend.position = "top")
dev.off()
################################
#点图测试
##############################
library(dplyr)
library(ggpubr)
setwd("D:\\OneDrive\\Retina_spatial\\fetal_analysis")
library(Seurat)
check_gene<-read.csv("D:\\OneDrive\\Retina_spatial\\fetal_analysis\\marker20231109.csv")
celllist<-read.csv("D:\\OneDrive\\Retina_spatial\\fetal_analysis\\更换cell types顺序.csv")
fetal_data$annotation <- factor(fetal_data$annotation,levels=unique(celllist$cluster))
pdf("dotplot_fetal_check_marker2.0.pdf",height=8,width=24)
DotPlot(fetal_data,assay = 'SCT', features =unique(check_gene$gene), group.by='annotation',
  cols = c("#FFE4D6", "#B0578D")) + 
RotatedAxis()+ theme(legend.position = "top")
dev.off()


check_gene<-read.csv("fetal Marker Genes 整理20230802 调整GJA1顺序.csv")
pdf("fetal Marker Genes 整理20230802 调整GJA1顺序.pdf",height=8,width=15)
DotPlot(fetal_data,assay = 'SCT', features =unique(check_gene$GENES), group.by='annotation',
  cols = c("#FFE4D6", "#B0578D")) + 
RotatedAxis()+ theme(legend.position = "top")
dev.off()
check_gene<-read.csv("fetal 转录因子 20230809 gene list.csv")
pdf("fetal 转录因子 20230809 gene list.pdf",height=8,width=15)
DotPlot(fetal_data,assay = 'SCT', features =unique(check_gene$GENES), group.by='annotation',
  cols = c("#FFE4D6", "#B0578D")) + 
RotatedAxis()+ theme(legend.position = "top")
dev.off()
check_gene<-read.csv("Markers_organoids 整理20230808-2.csv")
pdf("Markers_organoids 整理20230808-2.pdf",height=8,width=20)
DotPlot(fetal_data,assay = 'SCT', features =unique(check_gene$GENES), group.by='annotation',
  cols = c("#FFE4D6", "#B0578D")) + 
RotatedAxis()+ theme(legend.position = "top")
dev.off()


fetal_data<-readRDS("W19-fetal_scdata.rds")

check_gene<-read.csv("点图marker20231109-3.csv")
#fetal_data$annotation<- factor(fetal_data$annotation,levels=rev(c(unique(check_gene$cluster),"Unknown")))
fetal_data<-fetal_data[,fetal_data$annotation!="Unknown"]

fetal_data$annotation<- factor(fetal_data$annotation,levels=rev(unique(check_gene$cluster)))

pdf("点图marker20231109-3.pdf",height=8,width=20)
DotPlot(fetal_data,assay = 'SCT', features =unique(check_gene$gene), group.by='annotation',
  cols = c("#FFE4D6", "#B0578D")) + 
RotatedAxis()+ theme(legend.position = "top")
dev.off()


#################
#RPE轨迹
############################
library(Seurat)

#setwd("/share/pub/dengcy/STanalysis/fetusdata/analysis/")
setwd("D:\\OneDrive\\Retina_spatial\\fetal_analysis")
fetal_data<-readRDS("W19-fetal_scdata.rds")

sub_fetal_data<-subset(fetal_data,subset= annotation %in% c("RPE","RPE/PCs","RPE.stem-like cells"))
sub_fetal_data<-RunPCA(sub_fetal_data, assay = "SCT")
sub_fetal_data<-FindNeighbors(object=sub_fetal_data, reduction = 'pca', dims = 1:30)
sub_fetal_data<- RunUMAP(object = sub_fetal_data, assay ="SCT", reduction = 'pca', dims = 1:30,check_duplicates = FALSE)
sub_fetal_data<-RunTSNE(object = sub_fetal_data,assay = "SCT", reduction = 'pca', dims = 1:30,check_duplicates = FALSE)
library(slingshot, quietly = FALSE)
sce <- SingleCellExperiment(assays = List(counts = sub_fetal_data@assays$SCT@counts))

assays(sce)$norm <- sub_fetal_data@assays$SCT@data

#pca <- prcomp(t(log1p(assays(sce)$norm)), scale. = FALSE)
rd1 <- sub_fetal_data@reductions$pca[,1:2]

plot(rd1, col = rgb(0,0,0,.5), pch=16, asp = 1)

rd2 <- as.data.frame(Embeddings(sub_fetal_data, reduction = "umap"))
colnames(rd2) <- c("UMAP_1", "UMAP_2")

plot(rd2, col = rgb(0,0,0,.5), pch=16, asp = 1)

reducedDims(sce) <- SimpleList(PCA = rd2, UMAP = rd2)
#library(mclust, quietly = TRUE)
#cl1 <- Mclust(rd1)$classification
colData(sce)$annotation <- sub_fetal_data$annotation

sce <- slingshot(sce, clusterLabels = 'annotation', reducedDim = 'PCA',start.clus="RPE")

library(grDevices)
library(RColorBrewer)
colors <- colorRampPalette(c('#FBFACD','#DEBACE','#BA94D1','#7F669D'))(20)
plotcol <- colors[cut(sce$slingPseudotime_1, breaks=20)]
pdf("Slingshot_sub_RPE_fetal_data.pdf",width=4,height=4)
plot(reducedDims(sce)$UMAP, col = plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(sce), lwd=2, col='black')
plot(reducedDims(sce)$PCA, col = brewer.pal(9,'Set1')[sce$GMM], pch=16, asp = 1)
lines(SlingshotDataSet(sce), lwd=2, type = 'lineages', col = 'black')
dev.off()
save(sce,file='D:\\OneDrive\\Retina_spatial\\fetal_analysis\\monocle_analysi\\sce_slingPseudotime.RData')

sub_fetal_data$slingPseudotime<-sce$slingPseudotime_1

library(ggplot2)
pdf("umap_sub_RPE_slingPseudotime.pdf",width=6,height=5)
FeaturePlot(sub_fetal_data, features = "slingPseudotime", min.cutoff = 0, max.cutoff = 12)+scale_colour_gradient(low = '#FBFACD',high ='#7F669D')
dev.off()

pdf("umap_sub_RPE_annotation.pdf",width=7,height=5)
DimPlot(sub_fetal_data, group.by='annotation')
dev.off()


