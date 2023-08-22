##scATAC peaks heatmap plot
##Supplemental_Figure S2B
##input the singlecell data from fetal_celltypes_annotation.R
single_data<-readRDS("C:/Users/Lenovo/Desktop/final_version/annotation/final_annotation.rds")

single_data$wsnn_res_0.5_cell_type1<-single_data$wsnn_res_0.5_cell_type

##add another column
##hRSLCs is equal to Stem1, RPE_progenitors is equal to Stem2
single_data$wsnn_res_0.5_cell_type1<-as.character(single_data$wsnn_res_0.5_cell_type1)
single_data$wsnn_res_0.5_cell_type1[WhichCells(single_data,idents = c("hRSLCs"))]<-"Stem1"
single_data$wsnn_res_0.5_cell_type1[WhichCells(single_data,idents = c("RPE_progenitors"))]<-"Stem2"

##change the character to factor
single_data$wsnn_res_0.5_cell_type1<-factor(single_data$wsnn_res_0.5_cell_type1,levels=c("RPE","Stem2","Stem1","RPCs","PC_precursors","PCs","RGCs","ACs",
                                                                                         "HCs","BCs","MCs"))
Idents(single_data)<-single_data$wsnn_res_0.5_cell_type1
DefaultAssay(single_data)<-"peaks"
table(single_data$wsnn_res_0.5_cell_type1)
head(rownames(single_data))

##subset the peaks
library(stringr)
peaks<-as.data.frame(str_split_fixed(rownames(single_data),pattern = "-",n=3))
colnames(peaks)<-c("chr","start","end")
peaks[,2]<-as.numeric(peaks[,2])
peaks[,3]<-as.numeric(peaks[,3])

##input the markers for each celltype
library(openxlsx)
marker<-read.xlsx("C:/Users/Lenovo/Desktop/fetal Marker Genes.xlsx")


##import gtf from Gencode
library(rtracklayer)
gtf_data=import("C:/Users/Lenovo/Desktop/gencode.v44.annotation.gtf")
gtf_data = as.data.frame(gtf_data)

##subset the gtf_data which type belongs to "gene"
gtf_data <-gtf_data[gtf_data$type=="gene",]
##match marker and gtf_data
gtf_data<-gtf_data[match(marker$GENES,gtf_data$gene_name),]
gtf_data<-gtf_data[,c(1:5,12)]

##define the location of the promoter
for (i in 1:62) {
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
for (i in 1:62) {
  peaks1<-peaks[peaks$chr==gtf_data[i,1],]
  peaks1<-peaks1[(peaks1$start<gtf_data[i,7]&peaks1$end>gtf_data[i,7])|(peaks1$start>gtf_data[i,7]&peaks1$start<gtf_data[i,8]),]
  peaks1$gene<-rep(gtf_data[i,6],nrow(peaks1))
  df<-rbind(df,peaks1)
}

df$peaks<-paste(df$chr,df$start,df$end,sep = "-")

##AverageExpression
single_data1<-AverageExpression(single_data,assays = "peaks",slot = "data")
single_data1<-as.data.frame(single_data1)
single_data1<-single_data1[match(df$peaks,row.names(single_data1)),]

##add gene name, as label to sum
single_data1<-cbind(df$gene,single_data1)
single_data1<-as.data.frame(single_data1)

for (i in 2:12) {
  single_data1[,i]<-as.numeric(single_data1[,i])
}

##according to gene name, sum the peaks belong to same gene
single_data2<-aggregate(single_data1[,2:12],by=list(single_data1$`df$gene`),FUN = sum)

##adjust the ordering
single_data2<-single_data2[match(marker$GENES,single_data2$Group.1),]
single_data2<-na.omit(single_data2)
row.names(single_data2)<-single_data2$Group.1                    
single_data2<-single_data2[,-1]      

##scale
single_data2<-t(scale(t(single_data2),center = T,scale = T))

##phaetmap
library(pheatmap)
library(ggplot2)
range(single_data2)
bk<-c(seq(-3,-0.1,by=0.01),seq(0,3,by=0.01))

p<-pheatmap(t(single_data2),cluster_cols = F,cluster_rows = F,
            scale = "row",
            show_rownames = T,show_colnames = T,
            color = c(colorRampPalette(colors = c("#3A71AA","white"))(length(bk)/2),colorRampPalette(colors = c("white","#B22028"))(length(bk)/2)),
            breaks = bk,
            cellwidth = 10,cellheight = 8,
            border_color="black"
)
ggsave("peaks_heatmap_0.5kb.pdf",p,height = 10)
