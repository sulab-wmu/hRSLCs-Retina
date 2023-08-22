library("Seurat")
library("ggplot2")
library("patchwork")
library("dplyr")
#load single cell data
H9_D60<- readRDS("/share/pub/dengcy/STanalysis/H9_D60_2_stdata.rds")

H9_D60_Allmarkers <- Seurat::FindAllMarkers(object = H9_D60,
                                              assay = "SCT",
                                              slot = "data",
                                              verbose = TRUE,
                                              only.pos = TRUE)
write.csv(H9_D60_Allmarkers,file="D60_20_40_remove_Allmarkers.csv")


#library(SPOTlight)
sc_RDS_merged <- readRDS("/share/pub/dengcy/STanalysis/scdata/NEW_RDS_merged_clustered_v3.rds")

Idents(sc_RDS_merged) <- "cell_type_annotated"
sc_RDS_merged_markers <- Seurat::FindAllMarkers(object = sc_RDS_merged,
                                              assay = "SCT",
                                              slot = "data",
                                              verbose = TRUE,
                                              only.pos = TRUE)
saveRDS(object = sc_RDS_merged_markers, file = here::here("sc_RDS_merged_markers.rds"))