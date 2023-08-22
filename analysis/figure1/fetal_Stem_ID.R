##run Stem_ID
setwd("C:/Users/Lenovo/Desktop/final_version/Stem_ID/")
##devtools::install_github("dgrun/RaceID3_StemID2_package")
library(RaceID)
library(Seurat)

##use the compentropy function
edit(compentropy)

##import singlecell data from fetal_celltypes_annotation.R
single_data<-readRDS("C:/Users/Lenovo/Desktop/final_version/annotation/final_annotation.rds")

count<-as.matrix(single_data@assays$RNA@counts)
probs <- t(t(count)/apply(count, 2, sum))
entropy <- -apply(probs * log(probs + 1e-10)/log(nrow(count)), 2, sum)

##save the result
save(entropy,file="entropy.Rdata")