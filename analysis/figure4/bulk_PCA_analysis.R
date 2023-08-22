##bulk PCA analysis

##import the PCA data
library(openxlsx)
PCA<-read.xlsx("PCA.xlsx",sheet = "Sheet1")
row.names(PCA)<-PCA[,1]
PCA<-PCA[,-1]

##scale
PCA<-t(scale(t(PCA),center = T,scale = T))

##transpose
PCA<-t(PCA)

##load packages
#install.packages("FactoMineR")
library(FactoMineR)
library(ggpubr)
library(ggplot2)
library(ggthemes)

gene.pca <- PCA(PCA, ncp = 2, scale.unit = T, graph = FALSE)

##get the embedding
pca_sample <- data.frame(gene.pca$ind$coord[ ,1:2])
pca_sample$Sample=row.names(pca_sample)

##get the contribution of each pc
pca_eig1 <- round(gene.pca$eig[1,2], 2)
pca_eig2 <- round(gene.pca$eig[2,2],2 )

##add a column as label
pca_sample$Group<-c(rep("Normal retinogenesis",6),rep("Regenerative retinogenesis",5))

pca_sample<-pca_sample[,c(3,1,2,4)]
pca_sample$Group<-as.factor(pca_sample$Group)

##save the result
save(pca_sample,file = "PCA_result.Rdata")

