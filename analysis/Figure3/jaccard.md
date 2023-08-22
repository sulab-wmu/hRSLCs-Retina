```R
library(Seurat)
library(dplyr)
library(ggplot2)
library(Matrix)
library(scales)
library(reshape)
```

    Attaching SeuratObject
    
    Attaching sp
    
    
    Attaching package: ‘dplyr’
    
    
    The following objects are masked from ‘package:stats’:
    
        filter, lag
    
    
    The following objects are masked from ‘package:base’:
    
        intersect, setdiff, setequal, union
    
    
    
    Attaching package: ‘reshape’
    
    
    The following object is masked from ‘package:Matrix’:
    
        expand
    
    
    The following object is masked from ‘package:dplyr’:
    
        rename
    
    



```R
setwd('/share2/pub/zhouyj/zhouyj/Liu/analysis/Reg/drived_data')
```


```R
#Load development organoid snRMA-seq data, select  D0(60+0), D90
Reg <- readRDS('/share2/pub/zhouyj/zhouyj/Liu/analysis/Dev/drived_data/Final_Dev.rds') %>% subset(treatment %in% c('CMZ-0','CMZ-90'))
Reg$final_annotated_0807_2 <- Reg$raw_annotated6
Reg$tissue <- 'Organoid'
```


```R
#Load fetal scRNA-seq data
Fetal <- readRDS('/share2/pub/chenchg/chenchg/SingleCell/LiHui/final/final_annotation.rds') %>% subset(treatment == 'CMZ')
Fetal$final_annotated_0807_2 <- Fetal$wsnn_res_0.5_cell_type
Fetal$tissue <- 'Fetal'
```

    Loading required package: Signac
    



```R
Idents(Reg) <- Reg$final_annotated_0807_2
```


```R
Idents(Fetal) <- Fetal$wsnn_res_0.5_cell_type
```


```R
#Print matrix informations
Reg
```


    An object of class Seurat 
    99676 features across 7227 samples within 3 assays 
    Active assay: integrated (3000 features, 3000 variable features)
     2 other assays present: RNA, SCT
     4 dimensional reductions calculated: pca, umap, harmony, tsne



```R
Fetal
```


    An object of class Seurat 
    308828 features across 7562 samples within 4 assays 
    Active assay: RNA (60649 features, 0 variable features)
     3 other assays present: ATAC, peaks, SCT
     5 dimensional reductions calculated: pca, umap.rna, lsi, umap.atac, umap.biMod


# Compute marker genes


```R
#Set defaultAssay
DefaultAssay(Reg) <- 'SCT'
DefaultAssay(Fetal) <- 'SCT'
```


```R
Reg <- PrepSCTFindMarkers(Reg)
```

    Found 2 SCT models. Recorrecting SCT counts using minimum median counts: 1839
    



```R
Fetal <- PrepSCTFindMarkers(Fetal)
```

    Only one SCT model is stored - skipping recalculating corrected counts
    



```R
#Compute marker genes of fetal and organoid
Reg.sct.markers <- FindAllMarkers(Reg)
Fetal.sct.markers <- FindAllMarkers(Fetal)
```

    Calculating cluster hRSLCs
    
    Calculating cluster RPCs
    
    Calculating cluster PCs
    
    Calculating cluster BCs
    
    Calculating cluster RGCs
    
    Calculating cluster RPE_progenitors
    
    Calculating cluster ACs
    
    Calculating cluster RPE
    
    Calculating cluster PC_precursors
    
    Calculating cluster HCs
    
    Calculating cluster RPE
    
    Calculating cluster RPE_progenitors
    
    Calculating cluster hRSLCs
    
    Calculating cluster RPCs
    
    Calculating cluster PC_precursors
    
    Calculating cluster PCs
    
    Calculating cluster RGCs
    
    Calculating cluster ACs
    
    Calculating cluster HCs
    
    Calculating cluster BCs
    



```R
#Set cell type levels
my_levels2 <- c('RPE','RPE_progenitors','hRSLCs','RPCs','PC_precursors','PCs','RGCs','ACs','HCs','BCs')
```


```R
#Create empty matrix 
jacd_min <- matrix(,nrow=10,ncol=10)
rownames(jacd_min) <- paste('Organoid',my_levels2,sep = '_')
colnames(jacd_min) <- paste('Fetal',my_levels2,sep = '_')
```


```R
#Calculate Jaccard index
for(i in c(1:10)){
    reg <- Reg.sct.markers %>% subset(cluster == my_levels2[[i]]) %>% filter(avg_log2FC > 0) %>% filter(p_val_adj < 0.05)
    print('organoid')
    print(length(reg$gene))
    for(j in c(1:10)){
        fet <- Fetal.sct.markers %>% subset(cluster == my_levels2[[j]])%>% filter(avg_log2FC > 0)  %>% filter(p_val_adj < 0.05)
        print('fetal')
        print(length(fet$gene))
        int_min <- min(length(reg$gene),length(fet$gene))
        reg <- reg[c(1:int_min),]
        fet <- fet[c(1:int_min),]
        int <- length(intersect(reg$gene,fet$gene))
        jacd_min[i,j] = int/(length(reg$gene)+length(fet$gene)-int)   
    }
}
```

    [1] "organoid"
    [1] 326
    [1] "fetal"
    [1] 239
    [1] "fetal"
    [1] 264
    [1] "fetal"
    [1] 271
    [1] "fetal"
    [1] 302
    [1] "fetal"
    [1] 210
    [1] "fetal"
    [1] 285
    [1] "fetal"
    [1] 242
    [1] "fetal"
    [1] 349
    [1] "fetal"
    [1] 315
    [1] "fetal"
    [1] 43
    [1] "organoid"
    [1] 179
    [1] "fetal"
    [1] 239
    [1] "fetal"
    [1] 264
    [1] "fetal"
    [1] 271
    [1] "fetal"
    [1] 302
    [1] "fetal"
    [1] 210
    [1] "fetal"
    [1] 285
    [1] "fetal"
    [1] 242
    [1] "fetal"
    [1] 349
    [1] "fetal"
    [1] 315
    [1] "fetal"
    [1] 43
    [1] "organoid"
    [1] 213
    [1] "fetal"
    [1] 239
    [1] "fetal"
    [1] 264
    [1] "fetal"
    [1] 271
    [1] "fetal"
    [1] 302
    [1] "fetal"
    [1] 210
    [1] "fetal"
    [1] 285
    [1] "fetal"
    [1] 242
    [1] "fetal"
    [1] 349
    [1] "fetal"
    [1] 315
    [1] "fetal"
    [1] 43
    [1] "organoid"
    [1] 148
    [1] "fetal"
    [1] 239
    [1] "fetal"
    [1] 264
    [1] "fetal"
    [1] 271
    [1] "fetal"
    [1] 302
    [1] "fetal"
    [1] 210
    [1] "fetal"
    [1] 285
    [1] "fetal"
    [1] 242
    [1] "fetal"
    [1] 349
    [1] "fetal"
    [1] 315
    [1] "fetal"
    [1] 43
    [1] "organoid"
    [1] 115
    [1] "fetal"
    [1] 239
    [1] "fetal"
    [1] 264
    [1] "fetal"
    [1] 271
    [1] "fetal"
    [1] 302
    [1] "fetal"
    [1] 210
    [1] "fetal"
    [1] 285
    [1] "fetal"
    [1] 242
    [1] "fetal"
    [1] 349
    [1] "fetal"
    [1] 315
    [1] "fetal"
    [1] 43
    [1] "organoid"
    [1] 214
    [1] "fetal"
    [1] 239
    [1] "fetal"
    [1] 264
    [1] "fetal"
    [1] 271
    [1] "fetal"
    [1] 302
    [1] "fetal"
    [1] 210
    [1] "fetal"
    [1] 285
    [1] "fetal"
    [1] 242
    [1] "fetal"
    [1] 349
    [1] "fetal"
    [1] 315
    [1] "fetal"
    [1] 43
    [1] "organoid"
    [1] 203
    [1] "fetal"
    [1] 239
    [1] "fetal"
    [1] 264
    [1] "fetal"
    [1] 271
    [1] "fetal"
    [1] 302
    [1] "fetal"
    [1] 210
    [1] "fetal"
    [1] 285
    [1] "fetal"
    [1] 242
    [1] "fetal"
    [1] 349
    [1] "fetal"
    [1] 315
    [1] "fetal"
    [1] 43
    [1] "organoid"
    [1] 137
    [1] "fetal"
    [1] 239
    [1] "fetal"
    [1] 264
    [1] "fetal"
    [1] 271
    [1] "fetal"
    [1] 302
    [1] "fetal"
    [1] 210
    [1] "fetal"
    [1] 285
    [1] "fetal"
    [1] 242
    [1] "fetal"
    [1] 349
    [1] "fetal"
    [1] 315
    [1] "fetal"
    [1] 43
    [1] "organoid"
    [1] 134
    [1] "fetal"
    [1] 239
    [1] "fetal"
    [1] 264
    [1] "fetal"
    [1] 271
    [1] "fetal"
    [1] 302
    [1] "fetal"
    [1] 210
    [1] "fetal"
    [1] 285
    [1] "fetal"
    [1] 242
    [1] "fetal"
    [1] 349
    [1] "fetal"
    [1] 315
    [1] "fetal"
    [1] 43
    [1] "organoid"
    [1] 113
    [1] "fetal"
    [1] 239
    [1] "fetal"
    [1] 264
    [1] "fetal"
    [1] 271
    [1] "fetal"
    [1] 302
    [1] "fetal"
    [1] 210
    [1] "fetal"
    [1] 285
    [1] "fetal"
    [1] 242
    [1] "fetal"
    [1] 349
    [1] "fetal"
    [1] 315
    [1] "fetal"
    [1] 43



```R
#Check Jaccard matrix
jacd_min
```


<table class="dataframe">
<caption>A matrix: 10 × 10 of type dbl</caption>
<thead>
	<tr><th></th><th scope=col>Fetal_RPE</th><th scope=col>Fetal_RPE_progenitors</th><th scope=col>Fetal_hRSLCs</th><th scope=col>Fetal_RPCs</th><th scope=col>Fetal_PC_precursors</th><th scope=col>Fetal_PCs</th><th scope=col>Fetal_RGCs</th><th scope=col>Fetal_ACs</th><th scope=col>Fetal_HCs</th><th scope=col>Fetal_BCs</th></tr>
</thead>
<tbody>
	<tr><th scope=row>Organoid_RPE</th><td>0.33519553</td><td>0.154589372</td><td>0.081447964</td><td>0.012711864</td><td>0.016949153</td><td>0.012048193</td><td>0.02689487</td><td>0.01941748</td><td>0.02189781</td><td>0.00000000</td></tr>
	<tr><th scope=row>Organoid_RPE_progenitors</th><td>0.17763158</td><td>0.221843003</td><td>0.166123779</td><td>0.034682081</td><td>0.008450704</td><td>0.005617978</td><td>0.05604720</td><td>0.03468208</td><td>0.02285714</td><td>0.00000000</td></tr>
	<tr><th scope=row>Organoid_hRSLCs</th><td>0.08121827</td><td>0.121052632</td><td>0.260355030</td><td>0.073047859</td><td>0.021897810</td><td>0.007194245</td><td>0.07142857</td><td>0.06060606</td><td>0.03703704</td><td>0.00000000</td></tr>
	<tr><th scope=row>Organoid_RPCs</th><td>0.01718213</td><td>0.020689655</td><td>0.038596491</td><td>0.228215768</td><td>0.017182131</td><td>0.006802721</td><td>0.02422145</td><td>0.01369863</td><td>0.01369863</td><td>0.00000000</td></tr>
	<tr><th scope=row>Organoid_PC_precursors</th><td>0.00000000</td><td>0.000000000</td><td>0.004366812</td><td>0.133004926</td><td>0.074766355</td><td>0.017699115</td><td>0.00000000</td><td>0.00877193</td><td>0.01321586</td><td>0.00000000</td></tr>
	<tr><th scope=row>Organoid_PCs</th><td>0.01904762</td><td>0.019047619</td><td>0.011820331</td><td>0.004694836</td><td>0.126005362</td><td>0.304347826</td><td>0.02689487</td><td>0.01941748</td><td>0.01694915</td><td>0.01176471</td></tr>
	<tr><th scope=row>Organoid_RGCs</th><td>0.01500000</td><td>0.015000000</td><td>0.022670025</td><td>0.025252525</td><td>0.054545455</td><td>0.027848101</td><td>0.17002882</td><td>0.17341040</td><td>0.12777778</td><td>0.00000000</td></tr>
	<tr><th scope=row>Organoid_ACs</th><td>0.01481481</td><td>0.014814815</td><td>0.022388060</td><td>0.022388060</td><td>0.033962264</td><td>0.022388060</td><td>0.13223140</td><td>0.27441860</td><td>0.10040161</td><td>0.00000000</td></tr>
	<tr><th scope=row>Organoid_HCs</th><td>0.01901141</td><td>0.015151515</td><td>0.019011407</td><td>0.019011407</td><td>0.019011407</td><td>0.011320755</td><td>0.11203320</td><td>0.16017316</td><td>0.25821596</td><td>0.00000000</td></tr>
	<tr><th scope=row>Organoid_BCs</th><td>0.01345291</td><td>0.008928571</td><td>0.000000000</td><td>0.004444444</td><td>0.130000000</td><td>0.113300493</td><td>0.01801802</td><td>0.02727273</td><td>0.01801802</td><td>0.04878049</td></tr>
</tbody>
</table>




```R
#Convert matrix to list
jacd_mat_min <- melt(jacd_min)
```

    Warning message in type.convert.default(X[[i]], ...):
    “'as.is' should be specified by the caller; using TRUE”
    Warning message in type.convert.default(X[[i]], ...):
    “'as.is' should be specified by the caller; using TRUE”



```R
head(jacd_mat_min)
```


<table class="dataframe">
<caption>A data.frame: 6 × 3</caption>
<thead>
	<tr><th></th><th scope=col>X1</th><th scope=col>X2</th><th scope=col>value</th></tr>
	<tr><th></th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>1</th><td>Organoid_RPE            </td><td>Fetal_RPE</td><td>0.33519553</td></tr>
	<tr><th scope=row>2</th><td>Organoid_RPE_progenitors</td><td>Fetal_RPE</td><td>0.17763158</td></tr>
	<tr><th scope=row>3</th><td>Organoid_hRSLCs         </td><td>Fetal_RPE</td><td>0.08121827</td></tr>
	<tr><th scope=row>4</th><td>Organoid_RPCs           </td><td>Fetal_RPE</td><td>0.01718213</td></tr>
	<tr><th scope=row>5</th><td>Organoid_PC_precursors  </td><td>Fetal_RPE</td><td>0.00000000</td></tr>
	<tr><th scope=row>6</th><td>Organoid_PCs            </td><td>Fetal_RPE</td><td>0.01904762</td></tr>
</tbody>
</table>




```R
#Plot Jaccard similarity index
p <- ggplot(jacd_mat_min,aes(x=X2,y=X1,fill= value))+geom_tile()+scale_fill_gradient2(low="#67ADB7",high="#af2157")+ggtitle('min')+theme(axis.text.x = element_text(angle=90,hjust = 1,vjust = 1))
p
```


    
![png](jaccard_files/jaccard_19_0.png)
    



```R
#Save files
ggsave('../figures/jacad_heatmap.pdf',p)
write.csv(as.data.frame(jacd_min),'./jacad_mtx_0822.csv')
write.csv(as.data.frame(jacd_mat_min),'./jacad_list_0822.csv')
write.csv(Reg.sct.markers,'./Reg.sct.markers_cmz.csv')
write.csv(Fetal.sct.markers,'./Fetal.sct.markers_cmz.csv')
```

    Saving 6.67 x 6.67 in image
    



```R

```
