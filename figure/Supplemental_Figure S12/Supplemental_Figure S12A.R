#####Supplemental_Figure S12A
library(SCopeLoomR)
library(AUCell)
library(SCENIC)
library(KernSmooth)
library(RColorBrewer)
library(plotly)
library(BiocParallel)
library(grid)
library(ComplexHeatmap)
library(data.table)
library(Seurat)

setwd("C:/Users/Lenovo/Desktop/regeneration/pyscenic")

##organoids singlecell data
single_data<-readRDS("C:/Users/Lenovo/Desktop/regeneration/Final_Reg_CMZ_0807.rds")

cellinfo <- single_data@meta.data[,c('final_annotated_0807_2',"nFeature_RNA","nCount_RNA")]
colnames(cellinfo)=c('celltype','nGene' ,'nUMI')
cellTypes <-  as.data.frame(subset(cellinfo,select = 'celltype'))
selectedResolution <- "celltype"

##import the loom form organoids_pyscenic.sh
loom <- open_loom("organoid-genes_x_cells.csv.loom.out.auc_mtx.loom") 
regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")
regulons <- regulonsToGeneLists(regulons_incidMat)
class(regulons)
regulonAUC <- get_regulons_AUC(loom, column.attr.name='RegulonsAUC')
regulonAucThresholds <- get_regulon_thresholds(loom)

sub_regulonAUC <- regulonAUC
cellsPerGroup <- split(rownames(cellTypes), 
                       cellTypes[,selectedResolution])
regulonActivity_byGroup <- sapply(cellsPerGroup,
                                  function(cells) 
                                    rowMeans(getAUC(sub_regulonAUC)[,cells]))

regulonActivity_byGroup_Scaled <- t(scale(t(regulonActivity_byGroup),
                                          center = T, scale=T)) 
regulonActivity_byGroup_Scaled=na.omit(regulonActivity_byGroup_Scaled)


library(openxlsx)
TF<-read.xlsx("organoids TF list.xlsx")
TF$GENES<-paste(TF$GENES,"(+)",sep = "")
regulonActivity_byGroup_Scaled<-regulonActivity_byGroup_Scaled[match(TF$GENES,row.names(regulonActivity_byGroup_Scaled)),]
regulonActivity_byGroup_Scaled<-regulonActivity_byGroup_Scaled[,c(3:10,1:2)]
regulonActivity_byGroup_Scaled<-regulonActivity_byGroup_Scaled[c(5:20,1:4),]

library(pheatmap)
library(grid)
library(gtable)
range(regulonActivity_byGroup_Scaled)
bk<-c(seq(-2,-0.1,by=0.01),seq(0,2,by=0.01))

p<-pheatmap(t(regulonActivity_byGroup_Scaled),cluster_cols = F,cluster_rows = F,
            show_rownames = T,show_colnames = T,
            color = c(colorRampPalette(colors = c("#3A71AA","white"))(length(bk)/2),colorRampPalette(colors = c("white","#B22028"))(length(bk)/2)),
            breaks = bk,
            cellwidth = 15,cellheight = 15,
            border_color="black")
ggsave(filename = "regeneration_regulonActivity.pdf",p)