##Supplemental_Figure S3C5-C6
library(Seurat)

##import the result from fetal_scRNA_monocle3.R
cds<-readRDS("C:/Users/Lenovo/Desktop/final_version/monocle3_scRNA/scRNA_cds.rds")

##input the single_cell data from fetal_celltypes_annotation.R
data_fetal_scRNA<-readRDS("C:/Users/Lenovo/Desktop/final_version/annotation/final_annotation.rds")

##pseudotime
pseudotime<-pseudotime(cds, reduction_method ="UMAP")

data_fetal_scRNA$pseudotime <- pseudotime

head(data_fetal_scRNA$pseudotime)

Pseudotime <- as.numeric(data_fetal_scRNA$pseudotime)
celltypes <- as.character(data_fetal_scRNA$wsnn_res_0.5_cell_type)
data_pseudotime <- as.data.frame(cbind(celltypes,pseudotime))
names(data_pseudotime) <- c("celltypes_0.05","pseudotime")

### subset of four cell types in Lineage 1
data_pseudotime_subset <- data_pseudotime[which(data_pseudotime$celltypes %in% c("hRSLCs","PC_precursors","PCs","RPCs")),]
head(data_pseudotime_subset)

### subset of two cell types in Lineage 2
data_pseudotime_subset2 <- data_pseudotime[which(data_pseudotime$celltypes %in% c("RPE","RPE_progenitors")),]
head(data_pseudotime_subset2)

##C5
#data_pseudotime_subset$pseudotime
data_pseudotime_subset$celltypes_0.05<-factor(data_pseudotime_subset$celltypes_0.05,levels = c("hRSLCs","RPCs","PC_precursors","PCs"))
library(ggplot2)
#Density plot
p = ggplot(data_pseudotime_subset, aes(x=as.numeric(pseudotime),fill=celltypes_0.05)) + 
  theme_classic()+
  geom_density()+
  scale_fill_manual(values = c("#af2157","#BD8098","#C0BFDF","#67ADB7"))+
  theme(axis.title = element_text(size=18),
        axis.text = element_text(size=14),
        axis.text.x = element_text(angle = 35,
                                   hjust=1,
                                   size = 10,
                                   color = "black"))
p


###########----------------------------Lineage 2---------density plots-----------Supplemental Figure S3C----
#data_pseudotime_subset$pseudotime
data_pseudotime_subset2$celltypes_0.05<-factor(data_pseudotime_subset2$celltypes_0.05,levels = c("RPE","RPE_progenitors"))
library(ggplot2)
#Density plot
##C6
p = ggplot(data_pseudotime_subset2, aes(x=as.numeric(pseudotime),fill=celltypes_0.05)) + 
  theme_classic()+
  geom_density()+
  scale_fill_manual(values = c("#67ADB7","#BD8098"))+
  scale_y_sqrt(limits = c(0,100),breaks = c(0,3,5,10,20,30,50,100))+
  theme(axis.title = element_text(size=18),
        axis.text = element_text(size=14),
        axis.text.x = element_text(angle = 35,
                                   hjust=1,
                                   size = 10,
                                   color = "black"))
p +xlim(0,5)
