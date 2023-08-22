```R
library(clusterProfiler)
library(DESeq2)
library(ggplot2)
library(tximport)
#library(tximportData)
library(pheatmap)
library(RUVSeq)
library(dplyr)
library(ggplot2)
```

    
    
    Registered S3 method overwritten by 'ggtree':
      method      from 
      identify.gg ggfun
    
    clusterProfiler v4.2.0  For help: https://yulab-smu.top/biomedical-knowledge-mining-book/
    
    If you use clusterProfiler in published research, please cite:
    T Wu, E Hu, S Xu, M Chen, P Guo, Z Dai, T Feng, L Zhou, W Tang, L Zhan, X Fu, S Liu, X Bo, and G Yu. clusterProfiler 4.0: A universal enrichment tool for interpreting omics data. The Innovation. 2021, 2(3):100141
    
    
    Attaching package: ‘clusterProfiler’
    
    
    The following object is masked from ‘package:stats’:
    
        filter
    
    
    Loading required package: S4Vectors
    
    Loading required package: stats4
    
    Loading required package: BiocGenerics
    
    
    Attaching package: ‘BiocGenerics’
    
    
    The following objects are masked from ‘package:stats’:
    
        IQR, mad, sd, var, xtabs
    
    
    The following objects are masked from ‘package:base’:
    
        anyDuplicated, append, as.data.frame, basename, cbind, colnames,
        dirname, do.call, duplicated, eval, evalq, Filter, Find, get, grep,
        grepl, intersect, is.unsorted, lapply, Map, mapply, match, mget,
        order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
        rbind, Reduce, rownames, sapply, setdiff, sort, table, tapply,
        union, unique, unsplit, which.max, which.min
    
    
    
    Attaching package: ‘S4Vectors’
    
    
    The following object is masked from ‘package:clusterProfiler’:
    
        rename
    
    
    The following objects are masked from ‘package:base’:
    
        expand.grid, I, unname
    
    
    Loading required package: IRanges
    
    
    Attaching package: ‘IRanges’
    
    
    The following object is masked from ‘package:clusterProfiler’:
    
        slice
    
    
    Loading required package: GenomicRanges
    
    Loading required package: GenomeInfoDb
    
    Loading required package: SummarizedExperiment
    
    Loading required package: MatrixGenerics
    
    Loading required package: matrixStats
    
    
    Attaching package: ‘MatrixGenerics’
    
    
    The following objects are masked from ‘package:matrixStats’:
    
        colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
        colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
        colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
        colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
        colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
        colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
        colWeightedMeans, colWeightedMedians, colWeightedSds,
        colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
        rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
        rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
        rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
        rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
        rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
        rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
        rowWeightedSds, rowWeightedVars
    
    
    Loading required package: Biobase
    
    Welcome to Bioconductor
    
        Vignettes contain introductory material; view with
        'browseVignettes()'. To cite Bioconductor, see
        'citation("Biobase")', and for packages 'citation("pkgname")'.
    
    
    
    Attaching package: ‘Biobase’
    
    
    The following object is masked from ‘package:MatrixGenerics’:
    
        rowMedians
    
    
    The following objects are masked from ‘package:matrixStats’:
    
        anyMissing, rowMedians
    
    
    Loading required package: EDASeq
    
    Loading required package: ShortRead
    
    Loading required package: BiocParallel
    
    Loading required package: Biostrings
    
    Loading required package: XVector
    
    
    Attaching package: ‘Biostrings’
    
    
    The following object is masked from ‘package:base’:
    
        strsplit
    
    
    Loading required package: Rsamtools
    
    Loading required package: GenomicAlignments
    
    Loading required package: edgeR
    
    Loading required package: limma
    
    
    Attaching package: ‘limma’
    
    
    The following object is masked from ‘package:DESeq2’:
    
        plotMA
    
    
    The following object is masked from ‘package:BiocGenerics’:
    
        plotMA
    
    
    
    Attaching package: ‘dplyr’
    
    
    The following object is masked from ‘package:ShortRead’:
    
        id
    
    
    The following objects are masked from ‘package:GenomicAlignments’:
    
        first, last
    
    
    The following objects are masked from ‘package:Biostrings’:
    
        collapse, intersect, setdiff, setequal, union
    
    
    The following object is masked from ‘package:XVector’:
    
        slice
    
    
    The following object is masked from ‘package:Biobase’:
    
        combine
    
    
    The following object is masked from ‘package:matrixStats’:
    
        count
    
    
    The following objects are masked from ‘package:GenomicRanges’:
    
        intersect, setdiff, union
    
    
    The following object is masked from ‘package:GenomeInfoDb’:
    
        intersect
    
    
    The following objects are masked from ‘package:IRanges’:
    
        collapse, desc, intersect, setdiff, slice, union
    
    
    The following objects are masked from ‘package:S4Vectors’:
    
        first, intersect, rename, setdiff, setequal, union
    
    
    The following objects are masked from ‘package:BiocGenerics’:
    
        combine, intersect, setdiff, union
    
    
    The following objects are masked from ‘package:stats’:
    
        filter, lag
    
    
    The following objects are masked from ‘package:base’:
    
        intersect, setdiff, setequal, union
    
    



```R
#Set file path
file_pre="/share2/pub/zhouyj/zhouyj/Liu/20230810_bulk/bulk_fatsq/03-mapping/gene_expression"
output_pre=paste0(file_pre,"/", "output-seperate/")
if(!file.exists(output_pre)){
  dir.create(output_pre)
}
```


```R
#Create meta data for bulk RNA data sets
samples=data.frame(id=c("S01T0001","S01T0002","S01T0003","S01T0004","S01T0005","S01T0006"), Treatment=c("C-D20-NR-1","C-D20-NR-2","C-D20-NR-3","D80-CMZ-1","D80-CMZ-2","D80-CMZ-3"), batch=c("batch7","batch7","batch7","batch7","batch7","batch7"),tissue = c('NR_D20','NR_D20','NR_D20','CMZ_D80','CMZ_D80','CMZ_D80') )
samples$merge=paste0(samples$Treatment,"_",samples$id)
row.names(samples) = samples$merge
```


```R
#Load data
files2 <- file.path(file_pre, paste0("03_",samples$id, ".rsem.genes.results"))
names(files2) <- paste0(samples$Treatment,"_",samples$id)
txi.rsem2 <- tximport(files2, type = "rsem", txIn = FALSE, txOut = FALSE)
head(txi.rsem2$counts)
```

    reading in files with read_tsv
    
    1 
    2 
    3 
    4 
    5 
    6 
    
    



<table class="dataframe">
<caption>A matrix: 6 × 6 of type dbl</caption>
<thead>
	<tr><th></th><th scope=col>C-D20-NR-1_S01T0001</th><th scope=col>C-D20-NR-2_S01T0002</th><th scope=col>C-D20-NR-3_S01T0003</th><th scope=col>D80-CMZ-1_S01T0004</th><th scope=col>D80-CMZ-2_S01T0005</th><th scope=col>D80-CMZ-3_S01T0006</th></tr>
</thead>
<tbody>
	<tr><th scope=row>ENSG00000000003.15_TSPAN6</th><td>1684.00</td><td>1614.00</td><td>1719.00</td><td>1120.00</td><td>1228.00</td><td>1396.00</td></tr>
	<tr><th scope=row>ENSG00000000005.6_TNMD</th><td>   0.00</td><td>   0.00</td><td>   0.00</td><td>   1.00</td><td>   3.00</td><td>   0.00</td></tr>
	<tr><th scope=row>ENSG00000000419.14_DPM1</th><td> 869.34</td><td> 860.00</td><td>1054.55</td><td> 806.06</td><td> 947.93</td><td> 853.13</td></tr>
	<tr><th scope=row>ENSG00000000457.14_SCYL3</th><td> 309.71</td><td> 361.33</td><td> 392.09</td><td> 319.55</td><td> 311.80</td><td> 321.15</td></tr>
	<tr><th scope=row>ENSG00000000460.17_C1orf112</th><td> 272.29</td><td> 155.67</td><td> 223.91</td><td> 153.45</td><td> 142.20</td><td> 176.85</td></tr>
	<tr><th scope=row>ENSG00000000938.13_FGR</th><td>   4.00</td><td>   6.00</td><td>   1.00</td><td>   0.00</td><td>   0.00</td><td>   0.00</td></tr>
</tbody>
</table>




```R
#Set tximport for deseq
txi.rsem2$length[txi.rsem2$length == 0] <- 1
dds <- DESeqDataSetFromTximport(txi.rsem2, colData = samples, design = ~ tissue)
```

    Warning message in DESeqDataSet(se, design = design, ignoreRank):
    “some variables in design formula are characters, converting to factors”
    using counts and average transcript lengths from tximport
    



```R
dds
```


    class: DESeqDataSet 
    dim: 60649 6 
    metadata(1): version
    assays(2): counts avgTxLength
    rownames(60649): ENSG00000000003.15_TSPAN6 ENSG00000000005.6_TNMD ...
      ENSG00000288724.1_RP13-546I2.2 ENSG00000288725.1_RP11-413H22.3
    rowData names(0):
    colnames(6): C-D20-NR-1_S01T0001 C-D20-NR-2_S01T0002 ...
      D80-CMZ-2_S01T0005 D80-CMZ-3_S01T0006
    colData names(5): id Treatment batch tissue merge



```R
#Filter counts and perform deseq
dds <- dds[rowSums(counts(dds)) > 1,]
dds <- DESeq(dds)
```

    estimating size factors
    
    using 'avgTxLength' from assays(dds), correcting for library size
    
    estimating dispersions
    
    gene-wise dispersion estimates
    
    mean-dispersion relationship
    
    final dispersion estimates
    
    fitting model and testing
    



```R
dds
```


    class: DESeqDataSet 
    dim: 27938 6 
    metadata(1): version
    assays(6): counts avgTxLength ... H cooks
    rownames(27938): ENSG00000000003.15_TSPAN6 ENSG00000000005.6_TNMD ...
      ENSG00000288722.1_F8A1 ENSG00000288725.1_RP11-413H22.3
    rowData names(22): baseMean baseVar ... deviance maxCooks
    colnames(6): C-D20-NR-1_S01T0001 C-D20-NR-2_S01T0002 ...
      D80-CMZ-2_S01T0005 D80-CMZ-3_S01T0006
    colData names(5): id Treatment batch tissue merge



```R
#Filter results, group by tissue column
result = results(dds,alpha = 0.05, contrast=c("tissue", "NR_D20", "CMZ_D80"))
```


```R
result_data_filter <- subset(result,result$padj<0.05 & abs(result$log2FoldChange) >0.25)
```


```R
result_data_filter
```


    log2 fold change (MLE): tissue NR_D20 vs CMZ_D80 
    Wald test p-value: tissue NR D20 vs CMZ D80 
    DataFrame with 11693 rows and 6 columns
                                     baseMean log2FoldChange     lfcSE      stat
                                    <numeric>      <numeric> <numeric> <numeric>
    ENSG00000000971.16_CFH            104.401      -9.139562  1.211958  -7.54115
    ENSG00000001036.14_FUCA2         1156.438      -0.353381  0.108117  -3.26852
    ENSG00000001084.13_GCLC          1004.737      -1.866872  0.119846 -15.57724
    ENSG00000001461.17_NIPAL3         332.769      -0.688031  0.162987  -4.22138
    ENSG00000001497.18_LAS1L         1258.338       0.418763  0.108207   3.87003
    ...                                   ...            ...       ...       ...
    ENSG00000288637.1_RP11-559F19.4   8.04737       -6.67780  1.658633  -4.02609
    ENSG00000288656.1_RP11-666O2.7   26.94384       -3.59118  0.810595  -4.43030
    ENSG00000288663.1_RP11-680A11.7 484.15505        0.80981  0.188473   4.29669
    ENSG00000288674.1_RP5-1087E8.6   23.06754        2.26973  0.872597   2.60112
    ENSG00000288721.1_RP5-973N23.5   16.64926        1.84658  0.808738   2.28329
                                         pvalue        padj
                                      <numeric>   <numeric>
    ENSG00000000971.16_CFH          4.65840e-14 3.70918e-13
    ENSG00000001036.14_FUCA2        1.08111e-03 2.82643e-03
    ENSG00000001084.13_GCLC         1.03944e-54 5.14167e-53
    ENSG00000001461.17_NIPAL3       2.42809e-05 8.24179e-05
    ENSG00000001497.18_LAS1L        1.08823e-04 3.35626e-04
    ...                                     ...         ...
    ENSG00000288637.1_RP11-559F19.4 5.67124e-05 1.82632e-04
    ENSG00000288656.1_RP11-666O2.7  9.41009e-06 3.38291e-05
    ENSG00000288663.1_RP11-680A11.7 1.73371e-05 6.00958e-05
    ENSG00000288674.1_RP5-1087E8.6  9.29191e-03 2.01475e-02
    ENSG00000288721.1_RP5-973N23.5  2.24132e-02 4.42790e-02



```R
#Save data
write.csv(result_data_filter,'./output-seperate/result_data_filter.csv')
```
