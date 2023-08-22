##run Stem_ID
##organoids
library(RaceID)
library(Seurat)
setwd("C:/Users/Lenovo/Desktop/regeneration/Stem_ID")

##organoids singlecell data
single_data<-readRDS("C:/Users/Lenovo/Desktop/regeneration/Final_Reg_CMZ_0807.rds")

##subset the hRSLCs
Idents(single_data)<-single_data$final_annotated_0807_2
single_data<-subset(single_data,idents = "hRSLCs")

count<-as.matrix(single_data@assays$RNA@counts)
probs <- t(t(count)/apply(count, 2, sum))
entropy <- -apply(probs * log(probs + 1e-10)/log(nrow(count)), 2, sum)

##save the result
save(entropy,file="entropy.Rdata")