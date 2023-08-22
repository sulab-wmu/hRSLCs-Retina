##Figure 6B
##scATAC peaks heatmap plot
library(Seurat)

##import the singlecell data from fetal_celltypes_annotation.R
single_data<-readRDS("C:/Users/Lenovo/Desktop/final_version/annotation/final_annotation.rds")

Idents(single_data)<-single_data$wsnn_res_0.5_cell_type
DefaultAssay(single_data)<-"peaks"
table(single_data$wsnn_res_0.5_cell_type)
head(rownames(single_data))

##subset the peaks
library(stringr)
peaks<-as.data.frame(str_split_fixed(rownames(single_data),pattern = "-",n=3))
colnames(peaks)<-c("chr","start","end")
peaks[,2]<-as.numeric(peaks[,2])
peaks[,3]<-as.numeric(peaks[,3])

##import some import TFs
library(openxlsx)
marker<-read.xlsx("fetal TF list.xlsx")

##import gtf_data from Gencode
library(rtracklayer)
gtf_data=import("C:/Users/Lenovo/Desktop/gencode.v44.annotation.gtf")
gtf_data = as.data.frame(gtf_data)

##subset the gtf_data which type belongs to "gene"
gtf_data <-gtf_data[gtf_data$type=="gene",]
##match marker and gtf_data
gtf_data<-gtf_data[match(marker$GENES,gtf_data$gene_name),]
#gtf_data<-gtf_data[match(c("MECOM","MEIS1","TBX20"),gtf_data$gene_name),]
gtf_data<-gtf_data[,c(1:5,12)]

##define the location of the promoter
for (i in 1:23) {
  if(gtf_data[i,5]=="+"){
    gtf_data[i,7]<-gtf_data[i,2]-500
    gtf_data[i,8]<-gtf_data[i,2]+500
  }
  else{
    gtf_data[i,7]<-gtf_data[i,3]-500
    gtf_data[i,8]<-gtf_data[i,3]+500 
  }
}


##match the peak which is in or within the promoter(0.5kb)
df = as.data.frame(matrix(nrow=0,ncol=4)) 
for (i in 1:23) {
  peaks1<-peaks[peaks$chr==gtf_data[i,1],]
  peaks1<-peaks1[(peaks1$start<gtf_data[i,7]&peaks1$end>gtf_data[i,7])|(peaks1$start>gtf_data[i,7]&peaks1$start<gtf_data[i,8]),]
  peaks1$gene<-rep(gtf_data[i,6],nrow(peaks1))
  df<-rbind(df,peaks1)
}

df$peaks<-paste(df$chr,df$start,df$end,sep = "-")

##AverageExpression
single_data1<-AverageExpression(single_data,assays = "peaks",group.by = "wsnn_res_0.5_cell_type",slot = "data")
single_data1<-as.data.frame(single_data1)

##add gene name, as label to sum
single_data1<-single_data1[match(df$peaks,row.names(single_data1)),]
single_data1<-cbind(df$gene,single_data1)
single_data1<-as.data.frame(single_data1)

for (i in 2:12) {
  single_data1[,i]<-as.numeric(single_data1[,i])
}

single_data2<-aggregate(single_data1[,2:12],by=list(single_data1$`df$gene`),FUN = sum)

single_data2<-single_data2[match(marker$GENES,single_data2$Group.1),]
single_data2<-na.omit(single_data2)
row.names(single_data2)<-single_data2$Group.1                    
single_data2<-single_data2[,-1]      

##scale
single_data2<-t(scale(t(single_data2),center = T,scale = T))
colnames(single_data2)<-gsub("peaks.","",colnames(single_data2))
single_data2<-single_data2[,c(3:11,1:2)]
single_data2<-single_data2[c(4:21,1:3),]

##pheatmap
library(pheatmap)
library(ggplot2)
library(grid)
library(gtable)
range(single_data2)
bk<-c(seq(-2,-0.1,by=0.01),seq(0,2,by=0.01))

p<-pheatmap(t(single_data2),cluster_cols = F,cluster_rows = F,
            scale = "row",
            show_rownames = T,show_colnames = T,
            color = c(colorRampPalette(colors = c("#3A71AA","white"))(length(bk)/2),colorRampPalette(colors = c("white","#B22028"))(length(bk)/2)),
            breaks = bk,
            cellwidth = 15,cellheight = 15,
            border_color="black"
)
ggsave("fetal_TF_peaks_heatmap_0.5kb.pdf",p)






