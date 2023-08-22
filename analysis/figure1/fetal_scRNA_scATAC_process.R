library(DropletUtils)
library(Matrix)
library(dplyr)
library(tidyr)
library(patchwork)
library(Seurat)
library(Hmisc)
library(ggplot2)
library(rtracklayer)
library(sctransform)
library(clustree)
library(harmony)
library(cowplot)
library(stringr)
library(ggrepel)
#library(EnhancedVolcano)
library(ComplexHeatmap) #conda install -c conda-forge r-cairo, cairo-devel-cos6-x86_64, cairo
library(gplots)
library(clusterProfiler)
library(tibble)


file_pre=paste0(projectFloder, "/03-mapping/")
output_dir_pre=paste0(projectFloder, "/04-Seurat/")
output_pre_Seurat_obj=paste(output_dir_pre, "01_Seurat_obj_0_0/",sep="/")
output_pre=paste(output_dir_pre, "03_Res_0_0-seuIntegrate/",sep="/")
if(!file.exists(output_pre_Seurat_obj)){
  dir.create(output_pre_Seurat_obj)
}
if(!file.exists(output_pre)){
  dir.create(output_pre)
}

projectName=c("5_NR-0_1","5_CMZ-0_1","4_CMZ-05_1","4_CMZ-10_1","4_CMZ-20_1","4_CMZ-40_1" ) #,"3_CMZ-D45_1","3_NR-D45_1")
color_sample=c("#56B4E9","#CC6633","#CC6633","#CC6633","#CC6633","#CC6633","#CC6633", "#56B4E9")

data.10x = list(); # first declare an empty list in which to hold the feature-barcode matrices
scrna.list = list(); 
pancreas_merged=c();




#### step1: data prepration
print("step1: Data prepration ##################################################################")
for(idx in 1:length(projectName)){
  print(paste0("%%%%%%%%%%%%%%%%%: ",projectName[idx]));
  
  
  
  
  data.10x[[idx]] <- Read10X(data.dir = paste(file_pre,"/",dge_file[idx],sep="")); 
  #### step1.1: create Seurat object
  scrna.list[[idx]] = CreateSeuratObject(counts = data.10x[[idx]], min.cells=0, min.features=0, project=projectName[idx] );
  #scrna.list[[idx]] = CreateSeuratObject(counts = data.10x[[idx]], min.cells=3, min.features=200, project=projectName[idx] );
  saveRDS(scrna.list[[idx]], file = paste0(output_pre_Seurat_obj,  "/", projectName[idx], ".rds"))
  
  
  
  
  print("step1.1: raw cell number ####################################################")
  print(scrna.list[[idx]])
  
  #### step1.2: mitochondrial percentage and rRNA percentage
  scrna.list[[idx]][["percent.mt"]] <- PercentageFeatureSet(scrna.list[[idx]], pattern = "^MT-")
  scrna.list[[idx]][["percent.ribo"]] <- PercentageFeatureSet(scrna.list[[idx]], pattern = "^RP[SL][[:digit:]]")
  
  
  
  
  #### QC cutoff
  scrna.list[[idx]] <- subset(scrna.list[[idx]], subset = nFeature_RNA > 200  & nCount_RNA > 500 & percent.mt < 10 & percent.ribo < 10 )
  print("step1.2: after QC cell number ####################################################")
  print(scrna.list[[idx]])
  
  
  
  ## add meta-information
  batch_treatment=as.data.frame(t(as.data.frame(strsplit(as.character(scrna.list[[idx]]@meta.data$orig.ident), "_"))))
  names(batch_treatment)=c("batch","treatment")
  scrna.list[[idx]]@meta.data=cbind(scrna.list[[idx]]@meta.data,batch_treatment[,1:2])
}



#### step2: merge batches, add meta-information
print("step2: merge samples, add meta-information ##################################################################")

pancreas_merged <- merge(scrna.list[[1]], y = scrna.list[2:length(scrna.list)], project = "merged", merge.data = TRUE)

pancreas_merged <- NormalizeData(pancreas_merged, normalization.method = "LogNormalize", scale.factor = 10000)
pancreas_merged <- CellCycleScoring(pancreas_merged,  s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes, set.ident = TRUE)
pancreas_merged = SCTransform(pancreas_merged, vars.to.regress = c('nFeature_RNA', 'nCount_RNA',"percent.mt","percent.ribo","S.Score", "G2M.Score"), verbose = FALSE )

#######################################

saveRDS(pancreas_merged, file = paste0(output_pre,"/", "RDS_merged.rds"))



### dimension reduction
pancreas_merged <- FindVariableFeatures(pancreas_merged, assay = pancreas_merged@active.assay, selection.method = "vst", nfeatures = 3000)
pancreas_merged <- RunPCA(object = pancreas_merged, assay = pancreas_merged@active.assay,  npcs = 50)
# }
############## UMAP
pancreas_merged <- RunUMAP(object = pancreas_merged, assay = pancreas_merged@active.assay, dims = 1:50)

print("step3: Hamony: remove batch effect ##################################################################")

pancreas_merged <- RunHarmony(object = pancreas_merged,
                              assay.use = pancreas_merged@active.assay,
                              reduction = "pca",
                              dims.use = 1:50,
                              group.by.vars = "batch",
                              plot_convergence = TRUE)
reduction_name="harmony"

pancreas_merged <- RunTSNE(object = pancreas_merged, assay = pancreas_merged@active.assay, reduction = reduction_name, dims = 1:50)

############## UMAP
pancreas_merged <- RunUMAP(object = pancreas_merged, assay = pancreas_merged@active.assay, reduction = reduction_name, dims = 1:50)

############## cluster cells
pancreas_merged <- FindNeighbors(object = pancreas_merged, assay = pancreas_merged@active.assay, reduction = reduction_name, dims = 1:50)
pancreas_merged <- FindClusters(object = pancreas_merged, resolution = 0.4)


## add cell_type column
pancreas_merged$cell_type_annotated=pancreas_merged$seurat_clusters

saveRDS(pancreas_merged, file = paste(output_pre,  "/",  "RDS_merged_clustered.rds",sep=""))


#### identify cluster markers
clusterID = ""
conserved_markers  <- FindConservedMarkers(pancreas_merged, ident.1 = clusterID, grouping.var = "treatment", assay = pancreas_merged@active.assay, logfc.threshold = 0.25, only.pos = TRUE, verbose = FALSE, min.cells.group=1) 

pancreas_merged = readRDS( file = paste0("/share2/pub/xiongyc/xiongyc/project/scRNA/HuiLiu/data/10x_snRNA-seq-20211116_20220107_final.rds" ))
write.csv(pancreas_merged@meta.data, "/share2/pub/xiongyc/xiongyc/project/scRNA/HuiLiu/data/10x_snRNA-seq-20211116_20220107_final.rds.meta.data.csv" )

##combine scRNA and scATAC data
library(Signac)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
set.seed(1234)

output_pre=paste(output_dir_pre, "03_Res_0_0/",sep="/")



batch_treatment = c("1_NR_1","1_CMZ_1")
projectName = c("NR","CMZ")
color_sample=c("#56B4E9","#CC6633")


#### variates
scrna.list = list() 
pancreas_merged=c()
pbmc = c()


#### get gene annotations for hg38
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotation) <- "UCSC"
#### ATAC blacklist region from ENCODE project
blacklist_hg38 <- import("/share2/pub/xiongyc/xiongyc/project/scRNA/workspace/ATAC/Blacklist-master/lists/hg38-blacklist.v2.bed.gz")


for(idx in 1:length(projectName) ){
  
  
  
  # load the RNA and ATAC data
  counts <- Read10X_h5( paste0(file_pre, projectName[idx],"/outs/filtered_feature_bc_matrix.h5"))
  fragpath <- paste0(file_pre, projectName[idx],"/outs/atac_fragments.tsv.gz")
  
  
  # create a Seurat object containing the RNA adata
  scrna.list[[idx]] <- CreateSeuratObject(
    counts = counts$`Gene Expression`,
    assay = "RNA",
    project = projectName[idx]
  )
  scrna.list[[idx]]$orig.ident = batch_treatment[idx]
  scrna.list[[idx]][["percent.mt"]] <- PercentageFeatureSet(scrna.list[[idx]], pattern = "^MT-")
  scrna.list[[idx]][["percent.ribo"]] <- PercentageFeatureSet(scrna.list[[idx]], pattern = "^RP[SL][[:digit:]]")
  
  
  # create ATAC assay and add it to the object
  scrna.list[[idx]][["ATAC"]] <- CreateChromatinAssay(
    counts = counts$Peaks,
    sep = c(":", "-"),
    fragments = fragpath,
    annotation = annotation
  )
  
  print("create obj before QC:")
  print(scrna.list[[idx]])
  
  
  #### QC
  DefaultAssay(scrna.list[[idx]]) <- "ATAC"
  
  ## compute nucleosome signal score per cell
  scrna.list[[idx]] <- NucleosomeSignal(scrna.list[[idx]])
  ## compute TSS enrichment score per cell
  scrna.list[[idx]] <- TSSEnrichment(scrna.list[[idx]])
  
  
  #### call peak
  # call peaks using MACS2
  peaks <- CallPeaks(scrna.list[[idx]], macs2.path = "/share/pub/xiongyc/program/conda/install/envs/jinli/bin/macs2")
  
  # remove peaks on nonstandard chromosomes and in genomic blacklist regions
  peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
  peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)
  
  # quantify counts in each peak
  macs2_counts <- FeatureMatrix(
    fragments = Fragments(scrna.list[[idx]]),
    features = peaks,
    cells = colnames(scrna.list[[idx]])
  )
  
  # create a new assay using the MACS2 peak set and add it to the Seurat object
  scrna.list[[idx]][["peaks"]] <- CreateChromatinAssay(
    counts = macs2_counts,
    fragments = fragpath,
    annotation = annotation
  )
  
  
  ## Counting fraction of reads in peaks
  total_fragments <- CountFragments(fragpath)
  rownames(total_fragments) = total_fragments$CB
  scrna.list[[idx]]$fragments <- total_fragments[colnames(scrna.list[[idx]]), "frequency_count"]
  scrna.list[[idx]] <- FRiP( object = scrna.list[[idx]], assay = 'peaks', total.fragments = 'fragments' )
  print("FRiP calculate end")
  
  ## calculate blacklist fraction 
  scrna.list[[idx]]$blacklist_fraction <- FractionCountsInRegion(
    object = scrna.list[[idx]], 
    assay = 'peaks',
    regions = blacklist_hg38
  )
  
  ## plot
  
  
  
  saveRDS(scrna.list[[idx]], file = paste0(output_pre_Seurat_obj,"/", batch_treatment[idx], ".rds"))
  
  
  
  
  print("obj before QC:")
  print(scrna.list[[idx]])
  
  ## filter out low quality cells
  scrna.list[[idx]] <- subset(
    x = scrna.list[[idx]],
    subset = nCount_ATAC < 20000 &
      nCount_ATAC > 2000 &
      nCount_RNA < 20000 &
      nCount_RNA > 1000 &
      nucleosome_signal < 2 &
      TSS.enrichment > 2 &
      FRiP > 0.35 &
      blacklist_fraction < 0.05 &
      percent.mt < 5 &
      percent.ribo < 5
  )
  
  print("obj after QC:")
  print(scrna.list[[idx]])
  
  ## plot
  
  
}

#### step2: merge batches, add meta-information
print("step2: merge samples, add meta-information ##################################################################")

pbmc <- merge(scrna.list[[1]], y = scrna.list[2:length(scrna.list)], project = "merged", merge.data = TRUE)

## add meta-information
batch_treatment=as.data.frame(t(as.data.frame(strsplit(as.character(pbmc@meta.data$orig.ident), "_"))))
names(batch_treatment)=c("batch","treatment")
pbmc@meta.data=cbind(pbmc@meta.data,batch_treatment[,1:2])


saveRDS(pbmc, file = paste0(output_pre,"/", "RDS_merged.rds"))





########################################## clustering 

cluster_removed = FALSE;
if(!file_test("-f", paste0(output_pre, "/", "RDS_merged_clustered.rds") ) ){
  
  #### GEX
  DefaultAssay(pbmc) <- "RNA"
  #pbmc <- SCTransform(pbmc)
  pbmc <- CellCycleScoring(pbmc,  s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes, set.ident = TRUE)
  pbmc = SCTransform(pbmc, vars.to.regress = c('nFeature_RNA', 'nCount_RNA',"percent.mt","percent.ribo","S.Score", "G2M.Score"), verbose = FALSE, return.only.var.genes= FALSE )
  
  
  pbmc <- RunPCA(pbmc, npcs = 50)
  pbmc <- RunUMAP(pbmc, dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'UMAP.RNA_')
  
  #### ATAC
  DefaultAssay(pbmc) <- "peaks"
  pbmc <- FindTopFeatures(pbmc, min.cutoff = 5)
  pbmc <- RunTFIDF(pbmc)
  pbmc <- RunSVD(pbmc)
  pbmc <- RunUMAP(pbmc, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "UMAP.ATAC_")
  
  
  
  #### Joint UMAP visualization
  # build a joint neighbor graph using both assays
  pbmc <- FindMultiModalNeighbors(
    object = pbmc,
    reduction.list = list("pca", "lsi"), 
    dims.list = list(1:50, 2:50),
    modality.weight.name = "RNA.weight",
    verbose = TRUE
  )
  pbmc <- FindClusters(object = pbmc, resolution = 0.4, graph.name="wsnn")
  
  
  # build a joint UMAP visualization
  pbmc <- RunUMAP(
    object = pbmc,
    nn.name = "weighted.nn",
    #assay = "RNA",
    verbose = TRUE,
    reduction.name="umap.biMod",
    reduction.key = "UMAP.biMod_"
    #dims = 1:50
  )
  
  pbmc <- RunTSNE(
    object = pbmc_predicted,
    nn.name = "weighted.nn",
    #assay = "RNA",
    verbose = TRUE,
    reduction.name="tsne.biMod",
    reduction.key = "TSNE.biMod_",
    dims = 1:50
  )
  
  # print(plot_grid(plot_count,plot_nCount_ATAC,plot_tss,plot_nucl, ncol=2))
  
  p1 <- DimPlot(pbmc, reduction = "umap.rna",  label = TRUE, repel = TRUE) + ggtitle("RNA")
  p2 <- DimPlot(pbmc, reduction = "umap.atac", label = TRUE, repel = TRUE) + ggtitle("ATAC")
  p3 <- DimPlot(pbmc, reduction = "umap.biMod", label = TRUE, repel = TRUE) + ggtitle("biMod")
  p4 <- DimPlot(pbmc, reduction = "umap.biMod", group.by ="treatment", label = TRUE, repel = TRUE) + ggtitle("biMod")
  
  #DimPlot(pbmc,  label = TRUE, repel = TRUE, reduction = "umap_BiMod") #+ NoLegend()
  pdf(file=paste(output_pre, "01_UMAP-clustered.pdf",sep=""), width = 14, height = 14)
  print(p1 + p2 + p3 +p4 & theme(plot.title = element_text(hjust = 0.5)))
  dev.off()
  
  #### plot by TSNE
  p3 <- DimPlot(pbmc, reduction = "tsne.biMod", label = TRUE, repel = TRUE) + ggtitle("biMod")
  p4 <- DimPlot(pbmc, reduction = "tsne.biMod", group.by ="treatment", label = TRUE, repel = TRUE) + ggtitle("biMod")
  pdf(file=paste(output_pre, "01_TSNE-clustered.pdf",sep=""), width = 14, height = 7)
  print(p3 +p4 & theme(plot.title = element_text(hjust = 0.5)))
  dev.off()
  
  
  #### plot each cluster
  dir_each_cluster = paste0(output_pre,"/color_each_cluster")
  if(!file.exists(dir_each_cluster)){
    dir.create(dir_each_cluster)
  }
  oCol=c("0"="#DCDCDC", "1"="#DCDCDC", "2"="#DCDCDC", "3"="#DCDCDC", "4"="#DCDCDC", "5"="#DCDCDC", "6"="#DCDCDC", "7"="#DCDCDC", "8"="#DCDCDC", "9"="#DCDCDC", "10"="#DCDCDC", "11"="#DCDCDC", "12"="#DCDCDC", "13"="#DCDCDC", "14"="#DCDCDC", "15"="#DCDCDC", "16"="#DCDCDC", "17"="#DCDCDC", "18"="#DCDCDC", "19"="#DCDCDC", "20"="#DCDCDC", "21"="#DCDCDC","22"="#DCDCDC", "23"="#DCDCDC","24"="#DCDCDC","25"="#DCDCDC")
  for(ct in unique(pbmc$seurat_clusters)){
    
    t_col=c(ct='#FF0000', oCol )
    names(t_col)[1]=ct
    p=DimPlot(pbmc, group.by = c("seurat_clusters"), combine = FALSE, label = TRUE, repel=TRUE, reduction = "tsne.biMod",
              cols = t_col )
    pdf(file=paste(dir_each_cluster,"/",ct,"_TSNE.pdf",sep=""))
    print(p)
    dev.off()
  }
  
  saveRDS(pbmc, file = paste(output_pre,  "/",  "RDS_merged_clustered.rds",sep=""))
  
  
}else{
  
  pbmc = readRDS( file = paste0(output_pre,  "/",  "RDS_merged_clustered.rds" ))
  
}

#### label cell type name
pbmc_predicted = readRDS( file = paste0(output_pre,  "/",  "RDS_merged_clustered-predicted_celltype.rds" ))







