setwd("/share/pub/dengcy/liuhui_wk")
library("Seurat")
library("ggplot2")
library("patchwork")
library("dplyr")
library("hdf5r")
data_dir <- "/share2/pub/zhouyj/zhouyj/Liu/2023-5-17/result/H9_D60_2/2.Basic_analysis/2.3.h5_files/"
file_name <- "raw_feature_bc_matrix.h5"
scdata<- Load10X_Spatial(data.dir=data_dir,filename=file_name,slice="H9_D60_2")
scdata <- normalize_spacial(Spatial_data=scdata,lowqspot=0.01,mitper=25,geneExprMin=10,spot_meta="H9_D60")
scdata<-cluster_pca_umap(scdata, assay = "SCT",reduction="pca",cluster_res = 0.6)
saveRDS(scdata,file="H9_D60_2_stdata.rds")


cluster_pca_umap <- function(obj,assay=NULL, reduction,cluster_res = 0.3){
  #obj2 <- RunPCA(obj, assay = "SCT", reduction = "harmony",verbose = F)
  obj2 <- RunTSNE(object = obj,assay = assay, reduction = reduction, dims = 1:30,check_duplicates = FALSE)
  obj2 <- RunUMAP(object = obj2, assay =assay, reduction = reduction, dims = 1:30,check_duplicates = FALSE)
  obj2 <- FindNeighbors(object=obj2, assay = assay, reduction = reduction, dims = 1:30)
  obj2 <- FindClusters(object=obj2, resolution = cluster_res)
  return(obj2)
}

normalize_spacial<-function(Spatial_data,lowqspot=0.02,mitper=25,geneExprMin=15,spot_meta="RST2bei"){
    ###标准化
    ##Quantification of spots and genes were carried out by using subset function. 
    #We filtered out ~2% of low-quality spots. Spots with over 25% mitochondrial gene expression were also discarded
    feature_nums <- colSums(as.matrix(Spatial_data@assays$Spatial@counts) >0)
    mycut_feature <- as.numeric(quantile(feature_nums, lowqspot))
    count_nums <- colSums(as.matrix(Spatial_data@assays$Spatial@counts))
    mycut_count <- as.numeric(quantile(count_nums, lowqspot))
    Spatial_data_f1 <- subset(Spatial_data, subset = nFeature_Spatial>=mycut_feature | nCount_Spatial>=mycut_count)# & percent.mt < mitper )
    message("cut the count is success!")

    #Genes expressed in fewer than 15 spots were excluded.
    min.spots <- geneExprMin
    num.spots <- rowSums(as.matrix(Spatial_data@assays$Spatial@counts) > 0)
    genes.use <- names(num.spots[which(num.spots >= min.spots)])
    mykeepgene <- c(1:nrow(Spatial_data_f1))[rownames(Spatial_data_f1)%in%as.character(genes.use)]
    Spatial_data_f2 <- subset(Spatial_data_f1,features=mykeepgene)
    message("cut the spot is success!")
    # Filter out contaminated genes with specified names.
    #Genes related to hemoglobin (considerable variation from blood contents) and Y-chromosome linked genes were removed.
    filter_genes <-c("MALAT1", "SLC4A1", "KDM5D", "ANK1", "DDX3Y", "EIF2AK1", "HBQ1", "FTL", "GATA1", "KLF1", "USP9Y", "NFE2", "MT1G", "RPS4Y1", "HBZ", "GYPC", "HEMGN", "SLC25A37", "ALAS2", "EPB41", "AHSP", "GYPA", "UTY", "HBA2", "HBG2", "EIF1AY", "HBA1", "HBM", "HBE1", "HBG1", "MTRNR2L4", "HBB", "MTRNR2L5", "MTRNR2L8", "MTRNR2L10", "MTRNR2L3", "MTRNR2L1", "MTRNR2L7", "MTRNR2L12", "MTRNR2L11", "MTRNR2L13", "MTRNR2L6")
    keepgenes <- c(1:nrow(Spatial_data_f2))[!(rownames(Spatial_data_f2)%in%as.character(filter_genes))]
    Spatial_data_f3 <- subset(Spatial_data_f2,features=keepgenes)

    write.csv(Spatial_data_f3@meta.data, file = paste0(spot_meta,"-spots-metadata-clean.csv")) # Export the clean spots metadata
    message("cut the gene is success!")
    # SCT normalization
    #The clean expression matrix data were normalized using regularized negative binomial regression.
    Spatial_data_f3  <- SCTransform(Spatial_data_f3 , assay = "Spatial", return.only.var.genes = FALSE)
    DefaultAssay(Spatial_data_f3) <- "SCT"
    message("SCT normalization is success!")
    # Dimensionality reduction (PCA)
    #PCA was performed and the 10 most significant components was determined by the DimHeatmap and ElbowPlot function. 
    Spatial_data_f3 <- RunPCA(Spatial_data_f3)
    message("runpca is success!")
    message("DimHeatmap is success!")
    return(Spatial_data_f3)
}


