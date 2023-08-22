```R
library(Seurat)
library(Signac)
library(dplyr)
library(ggplot2)
library(Matrix)
#library(paletteer)
library(patchwork)
library(scales)
library(cowplot)
```

    Attaching SeuratObject
    
    Attaching sp
    
    
    Attaching package: ‘dplyr’
    
    
    The following objects are masked from ‘package:stats’:
    
        filter, lag
    
    
    The following objects are masked from ‘package:base’:
    
        intersect, setdiff, setequal, union
    
    
    
    Attaching package: ‘cowplot’
    
    
    The following object is masked from ‘package:patchwork’:
    
        align_plots
    
    



```R
setwd('/share2/pub/zhouyj/zhouyj/Liu/analysis/Reg/drived_data')
```


```R
#Load fetal scRNA-seq data
Fetal <- readRDS('/share2/pub/chenchg/chenchg/SingleCell/LiHui/final/final_annotation.rds')
table(Fetal@meta.data$treatment)
```


    
     CMZ   NR 
    7562 8978 



```R
#Select tissue CMZ
Fetal <- Fetal %>% subset(treatment == 'CMZ')
Fetal$final_annotated_0807_2 <- Fetal$wsnn_res_0.5_cell_type
Fetal$tissue <- 'Fetal'
```


```R
#Load development organoid snRNA-seq data and select tissue CMZ D0(60+0), D90
Reg2 <- readRDS('/share2/pub/zhouyj/zhouyj/Liu/analysis/Dev/drived_data/Final_Dev.rds') %>% subset(treatment %in% c('CMZ-0','CMZ-90'))
Reg2$final_annotated_0807_2 <- Reg2$raw_annotated6
Reg2$tissue <- 'Organoid'
```


```R
#Set defaultAssay RNA
DefaultAssay(Reg2) <- 'RNA'
DefaultAssay(Fetal) <- 'RNA'
```


```R
#Preparation for merge data
ifnb.list <- list()
ifnb.list[[1]] <- Reg2
ifnb.list[[2]] <- Fetal
```


```R
ifnb.list <- lapply(X = ifnb.list, FUN = SCTransform)
features <- SelectIntegrationFeatures(object.list = ifnb.list, nfeatures = 3000)
ifnb.list <- PrepSCTIntegration(object.list = ifnb.list, anchor.features = features)
ifnb.list <- lapply(X = ifnb.list, FUN = RunPCA, features = features)
```

    Calculating cell attributes from input UMI matrix: log_umi
    
    Variance stabilizing transformation of count matrix of size 26726 by 7227
    
    Model formula is y ~ log_umi
    
    Get Negative Binomial regression parameters per gene
    
    Using 2000 genes, 5000 cells
    


      |======================================================================| 100%


    Found 84 outliers - those will be ignored in fitting/regularization step
    
    
    Second step: Get residuals using fitted parameters for 26726 genes
    


      |======================================================================| 100%


    Computing corrected count matrix for 26726 genes
    


      |======================================================================| 100%


    Calculating gene attributes
    
    Wall clock passed: Time difference of 2.04239 mins
    
    Determine variable features
    
    Place corrected count matrix in counts slot
    
    Centering data matrix
    
    Set default assay to SCT
    
    Calculating cell attributes from input UMI matrix: log_umi
    
    Variance stabilizing transformation of count matrix of size 26562 by 7562
    
    Model formula is y ~ log_umi
    
    Get Negative Binomial regression parameters per gene
    
    Using 2000 genes, 5000 cells
    


      |======================================================================| 100%


    Found 74 outliers - those will be ignored in fitting/regularization step
    
    
    Second step: Get residuals using fitted parameters for 26562 genes
    


      |======================================================================| 100%


    Computing corrected count matrix for 26562 genes
    


      |======================================================================| 100%


    Calculating gene attributes
    
    Wall clock passed: Time difference of 1.973898 mins
    
    Determine variable features
    
    Place corrected count matrix in counts slot
    
    Centering data matrix
    
    Set default assay to SCT
    
    PC_ 1 
    Positive:  DCBLD2, CDH6, GEM, TENM3, HES1, EFNA5, AKAP12, LRP1B, NECTIN3-AS1, SQSTM1 
    	   HMGA2, RASAL2, PLXDC2, MT-RNR2, SLC5A3, TOX, ARPP21, PCDH11X, CNTN5, SERPINE1 
    	   CEMIP2, GLIS3, RP11-71N10.1, COL25A1, TRIO, SFRP2, RP11-152K4.2, WEE1, HS6ST3, EGFR 
    Negative:  FSTL5, EYS, NRG1, DPP10, CADM2, CADPS, ROBO2, AKAP9, OTX2, NTM 
    	   PDC, MYO3B, OTX2-AS1, PDE1C, ZNF385B, EGFLAM, PEX5L, TAFA4, PCBP3, STXBP5L 
    	   MPP4, MIR124-1HG, AGAP1, TMEM108, RIMS2, SEMA6D, NEGR1, ALK, NLK, SLC8A1 
    PC_ 2 
    Positive:  DCBLD2, CDH6, FSTL5, EYS, GEM, AKAP12, ARPP21, MT-RNR2, NRCAM, TRIO 
    	   TENM3, SERPINE1, LRP1B, HIVEP3, DPP10, RASAL2, TCERG1L, HES1, MAML2, HS6ST3 
    	   CADM2, ROBO2, HMGA2, WEE1, XKR4, RORB, RP11-152K4.2, PCDH11X, CADPS, COL25A1 
    Negative:  TRPM3, COL8A1, TSHZ2, TRPM1, DCT, CD96, PTPRT, CNGB3, LRMDA, IMMP2L 
    	   BEST1, NEAT1, RP11-779P15.2, CHSY3, ENPP2, RP4-678D15.1, CPEB4, TIMP3, RLBP1, PCAT1 
    	   TYR, LRRTM4, IGFBP5, KCNQ1OT1, PLD5, ELN, LIMD1, SLC26A7, SILC1, PRKG1 
    PC_ 3 
    Positive:  RORB, NRXN3, RP11-78L19.1, MIAT, KCTD16, FSTL4, FAM155A, CNTNAP2, NCAM1, LINC00461 
    	   PTPRD, ADGRL2, ZNF385D, CUX2, MOB3B, CKB, STMN2, PTPRZ1, MDGA2, HES5 
    	   LINGO2, PHYHIPL, KLHL1, GRIA4, MIR503HG, LINC02334, ZNF804A, EPHA5, EBF1, RUNX1T1 
    Negative:  DCT, COL8A1, TRPM3, TRPM1, TSHZ2, CD96, EYS, GEM, FSTL5, DCBLD2 
    	   PTPRT, SGCD, LRMDA, CNGB3, DLGAP1, GPC6, CDH6, OTX2-AS1, KCNB2, ENPP2 
    	   RP11-779P15.2, CPEB4, BEST1, TYR, MITF, PHACTR2, TIMP3, HMGA2, LHFPL6, PLD5 
    PC_ 4 
    Positive:  NR2F1, RP11-78L19.1, MOB3B, MIR503HG, CRB1, CNTLN, CHSY3, ADGRL2, PCDH11X, PTPRZ1 
    	   LINC02334, CKB, SORCS1, SLC26A7, CLEC2A, CLVS1, LINC02232, ADAMTS6, ARRDC3-AS1, TLL1 
    	   LINC00461, COL25A1, LINC01697, PRKG1, SBF2, RTN4, ZFP36L2, HES6, NFIB, CUX2 
    Negative:  CNTNAP2, FSTL4, KCNIP4, ZNF804A, STMN2, EBF1, KLHL1, CTNNA2, AKAP6, LRFN5 
    	   CSMD1, FGF14, DCC, GRIA2, PTPRD, MYT1L, CNTN4, RUNX1T1, GRIP1, NRXN1 
    	   ELAVL2, CADM2, DPP6, ISL1, SLIT2, OLFM3, CELF4, CTC-340A15.2, ERBB4, ONECUT1 
    PC_ 5 
    Positive:  COL8A1, TSHZ2, DCT, RORB, PCDH11X, TRPM1, CD96, CNGB3, BEST1, RP4-678D15.1 
    	   RP11-779P15.2, PTPRZ1, IMMP2L, RLBP1, LRMDA, DLGAP1, ITGB8, CKB, COL25A1, ZNF385D 
    	   LINC00511, MOB3B, GRB10, LIMD1, LINC00461, PCAT1, HES1, PLD5, CPEB4, FRMD5 
    Negative:  CHSY3, PCDH9, MECOM, EFNA5, RMST, GPC6, CTC-575N7.1, RP11-120I21.2, SLC26A7, CCSER1 
    	   HPSE2, LHFPL3, ADAMTS19, COL9A1, HAS2, UNC5C, RP11-305F18.1, CDH12, MYOCD, NR2F1 
    	   SGCZ, SEMA3E, RELN, SYNPR, PCDH7, GRIA4, CACNA1D, NRG1, NEAT1, LINC01122 
    
    PC_ 1 
    Positive:  NRXN3, RALYL, DSCAM, NRG3, DCC, MDGA2, GRM8, ASIC2, DLGAP1, ZFHX3 
    	   MIR181A1HG, PTPRD, NEGR1, PCDH9, NRP1, ATRNL1, RUNX1T1, FGF14, MYT1L, FRMD5 
    	   CTNNA2, CSMD1, AFF3, MEIS2, PDE3A, ENOX1, KCNIP4, GRIA4, SNTG1, PTPRT 
    Negative:  EYS, NRG1, MARCHF1, FSTL5, ZNF385B, TMEM108, IQCJ-SCHIP1, PRDM1, EPB41L2, RP1 
    	   ADGRV1, GPC6, ANKRD33B, ANO2, EGFLAM, IMPG2, EPS8, PDE1C, USH2A, DPP10 
    	   OTX2, HTR1F, OTX2-AS1, PDC, NLK, TMEM244, PIK3R1, MPPED2, SLC35F3, DLEU2 
    PC_ 2 
    Positive:  TRPM3, SLC38A11, ADGRL2, SORCS1, FOXP2, GABRG3, SLC4A4, LRMDA, LGR4, PARD3B 
    	   SGCD, PRKG1, ROBO1, PCDH11X, NAV2, TF, GLI3, WIF1, PIP5K1B, TYR 
    	   PCDH11Y, CLVS1, MIR924HG, SPP1, ST6GALNAC3, CHSY3, NFIA, GULP1, MIR99AHG, SFRP2 
    Negative:  EYS, NRXN3, NRG1, DSCAM, RALYL, DCC, NRXN1, GRIP1, CADM2, MDGA2 
    	   MARCHF1, MEIS2, FSTL5, ZNF385B, ASIC2, CADPS, SNTG1, NEGR1, NRG3, CNTNAP2 
    	   PPFIA2, MIR181A1HG, TMEM108, RIMS2, ROBO2, ZFHX3, PRDM1, IQCJ-SCHIP1, RP1, RUNX1T1 
    PC_ 3 
    Positive:  ADGRL2, SORCS1, RORB, LGR4, LSAMP, ZNF385D, ST6GALNAC3, PCDH11X, PCDH11Y, TF 
    	   RP11-384F7.2, NAV2, WIF1, SPP1, CRB1, LINC00461, RTN4, NCAM1, LINC00511, NRXN3 
    	   RP11-436K8.1, LINC01697, COL25A1, FOXP2, ZHX2, MOB3B, SLC1A3, KCTD8, CREB5, CUX2 
    Negative:  SLC38A11, SGCD, LRMDA, SLC4A4, ROBO1, TRPM3, KCND2, TYR, DCT, AC013463.2 
    	   SGK1, CCBE1, TIMP3, DLGAP1, PLD5, TMEM132C, CPA6, LHFPL6, NIBAN1, MITF 
    	   ATP6V1C2, PLD1, RP11-141F7.1, BNC2, RBM20, IGF1, OTX2-AS1, IQGAP2, SPOCK3, EYS 
    PC_ 4 
    Positive:  NKAIN2, ONECUT1, TNR, TENM2, CNTN4, RYR3, MEGF10, CSMD1, MALT1, PCDH9 
    	   CPNE4, TMEFF2, TANC1, PROX1, FRMPD4, GRIA3, KCNIP1, ONECUT2, PDE4B, CACNG3 
    	   NRXN3, DRD2, ONECUT3, RP11-126O1.4, NTRK2, TPM3, ZNF804A, MAGI2, RP11-209K10.2, SYN3 
    Negative:  DSCAM, DCC, RALYL, NRG3, MDGA2, MEIS2, SNTG1, GRM8, ZFHX3, DLGAP1 
    	   ADGRL3, FGF14, PTPRT, NEGR1, ATRNL1, POU6F2, IL1RAPL2, CADM2, GALNTL6, PTCHD4 
    	   PBX1, NRXN1, GALNT17, IL1RAPL1, RBFOX1, FAT3, FRMD5, KCNIP4, PCDH7, SGCZ 
    PC_ 5 
    Positive:  IQCJ-SCHIP1, GRIK1, NFIB, NLK, RP11-120J1.2, NRXN3, FAM155A, SOX6, SEMA5A, NCAM2 
    	   FABP7, NOL4, PIK3R1, DPP10, PEX5L, KCTD16, PAG1, RORB, WDR72, PCSK2 
    	   CD44, GRID2, ST18, CTNND2, LRRC4C, KCNQ5, NWD2, TMTC2, SLC38A11, PDE4D 
    Negative:  SEMA6D, GPC5, RP11-13N12.1, EYS, MSR1, GRIP1, GPC6, PRKN, THRB, EGFLAM 
    	   TTN, OTX2-AS1, RP11-152L20.3, SLC8A1, MCC, PCDH9, ROBO2, FMN1, AGBL4, CHSY3 
    	   GRIK4, GRAMD1C, LHFPL3, CCDC141, IMPG1, MPP4, WWC1, COBL, MECOM, FSTL5 
    



```R
#Merge by integration anchors
immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, normalization.method = "SCT",
    anchor.features = features, dims = 1:50, reduction = "rpca", k.anchor = 20)
immune.combined.sct <- IntegrateData(anchorset = immune.anchors, normalization.method = "SCT", dims = 1:50)
immune.combined.sct <- RunPCA(immune.combined.sct, verbose = FALSE)
immune.combined.sct <- RunUMAP(immune.combined.sct, reduction = "pca", dims = 1:50)
```

    Computing within dataset neighborhoods
    
    Finding all pairwise anchors
    
    Projecting new data onto SVD
    
    Projecting new data onto SVD
    
    Finding neighborhoods
    
    Finding anchors
    
    	Found 12494 anchors
    
    Merging dataset 1 into 2
    
    Extracting anchors for merged samples
    
    Finding integration vectors
    
    Finding integration vector weights
    
    Integrating data
    
    Warning message:
    “The default method for RunUMAP has changed from calling Python UMAP via reticulate to the R-native UWOT using the cosine metric
    To use Python UMAP via reticulate, set umap.method to 'umap-learn' and metric to 'correlation'
    This message will be shown once per session”
    15:49:02 UMAP embedding parameters a = 0.9922 b = 1.112
    
    15:49:02 Read 14789 rows and found 50 numeric columns
    
    15:49:02 Using Annoy for neighbor search, n_neighbors = 30
    
    15:49:02 Building Annoy index with metric = cosine, n_trees = 50
    
    0%   10   20   30   40   50   60   70   80   90   100%
    
    [----|----|----|----|----|----|----|----|----|----|
    
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    |
    
    15:49:05 Writing NN index file to temp file /tmp/RtmptMrHpw/file2824253ecc524
    
    15:49:05 Searching Annoy index using 1 thread, search_k = 3000
    
    15:49:10 Annoy recall = 100%
    
    15:49:11 Commencing smooth kNN distance calibration using 1 thread
    
    15:49:13 Initializing from normalized Laplacian + noise
    
    15:49:13 Commencing optimization for 200 epochs, with 676510 positive edges
    
    15:49:32 Optimization finished
    



```R
#Check umap
DimPlot(immune.combined.sct, reduction = "umap", group.by = "final_annotated_0807_2", label = TRUE, repel = TRUE)
```


    
![png](Merge_fetal_organoid_d60_d90_files/Merge_fetal_organoid_d60_d90_9_0.png)
    



```R
#Check batch effects
DimPlot(immune.combined.sct, reduction = "umap", group.by = "tissue", label = TRUE, repel = TRUE)
```


    
![png](Merge_fetal_organoid_d60_d90_files/Merge_fetal_organoid_d60_d90_10_0.png)
    



```R
#Run harmony to remove batch effects
library(harmony)
```

    Loading required package: Rcpp
    



```R
pancreas_merged <- RunHarmony(object = immune.combined.sct,
                                  assay.use = immune.combined.sct@active.assay,
                                  reduction = "pca",
                                  dims.use = 1:50,
                                  group.by.vars = "treatment",
                                  plot_convergence = TRUE)
reduction_name="harmony"

pancreas_merged <- RunTSNE(object = pancreas_merged, assay = pancreas_merged@active.assay, reduction = reduction_name, dims = 1:50)
  
############## UMAP
pancreas_merged <- RunUMAP(object = pancreas_merged, assay = pancreas_merged@active.assay, reduction = reduction_name, dims = 1:50)
  
############## cluster cells
pancreas_merged <- FindNeighbors(object = pancreas_merged, assay = pancreas_merged@active.assay, reduction = reduction_name, dims = 1:50)
pancreas_merged <- FindClusters(object = pancreas_merged, resolution = 0.4)

```

    Harmony 1/10
    
    Harmony 2/10
    
    Harmony 3/10
    
    Harmony 4/10
    
    Harmony 5/10
    
    Harmony 6/10
    
    Harmony 7/10
    
    Harmony 8/10
    
    Harmony 9/10
    
    Harmony 10/10
    
    Warning message:
    “Invalid name supplied, making object name syntactically valid. New object name is Seurat..ProjectDim.integrated.harmony; see ?make.names for more details on syntax validity”
    15:51:07 UMAP embedding parameters a = 0.9922 b = 1.112
    
    15:51:07 Read 14789 rows and found 50 numeric columns
    
    15:51:07 Using Annoy for neighbor search, n_neighbors = 30
    
    15:51:07 Building Annoy index with metric = cosine, n_trees = 50
    
    0%   10   20   30   40   50   60   70   80   90   100%
    
    [----|----|----|----|----|----|----|----|----|----|
    
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    |
    
    15:51:09 Writing NN index file to temp file /tmp/RtmptMrHpw/file282427bf0f7c
    
    15:51:09 Searching Annoy index using 1 thread, search_k = 3000
    
    15:51:13 Annoy recall = 100%
    
    15:51:14 Commencing smooth kNN distance calibration using 1 thread
    
    15:51:16 Initializing from normalized Laplacian + noise
    
    15:51:17 Commencing optimization for 200 epochs, with 686100 positive edges
    
    15:51:36 Optimization finished
    
    Computing nearest neighbor graph
    
    Computing SNN
    


    Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    
    Number of nodes: 14789
    Number of edges: 683835
    
    Running Louvain algorithm...
    Maximum modularity in 10 random starts: 0.9251
    Number of communities: 13
    Elapsed time: 2 seconds



    
![png](Merge_fetal_organoid_d60_d90_files/Merge_fetal_organoid_d60_d90_12_2.png)
    



```R
DimPlot(pancreas_merged, reduction = "umap", group.by = "final_annotated_0807_2", label = TRUE, repel = TRUE)
```


    
![png](Merge_fetal_organoid_d60_d90_files/Merge_fetal_organoid_d60_d90_13_0.png)
    



```R
DimPlot(pancreas_merged, reduction = "umap", group.by = "tissue", label = TRUE, repel = TRUE)
```


    
![png](Merge_fetal_organoid_d60_d90_files/Merge_fetal_organoid_d60_d90_14_0.png)
    



```R
#Save data
p <- DimPlot(pancreas_merged, reduction = "umap", group.by = "final_annotated_0807_2", label = TRUE, repel = TRUE)
ggsave('../figures/fetal_cmz_60_90.pdf',p)
saveRDS(pancreas_merged,'./fetal_cmz_60_90_harmnoy.rds')
```

    Saving 6.67 x 6.67 in image
    

