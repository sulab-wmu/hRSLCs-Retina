import stlearn as st
import scanpy as sc
from pathlib import Path
import matplotlib.pyplot as plt
import os
import scanpy as sc
import pandas as pd
os.chdir("/share/pub/dengcy/STanalysis/")
st.settings.set_figure_params(dpi=180)
BASE_PATH = Path("/share/pub/dengcy/STanalysis/H9_D60_2/")
# spot tile is the intermediate result of image pre-processing
TILE_PATH = Path("/share/pub/dengcy/STanalysis/H9_D60_2/")
TILE_PATH.mkdir(parents=True, exist_ok=True)
# output path
OUT_PATH = Path("/share/pub/dengcy/STanalysis/H9_D60_2/")
OUT_PATH.mkdir(parents=True, exist_ok=True)
adata = st.Read10X(BASE_PATH)
# pre-processing for gene count table
st.pp.filter_genes(adata,min_cells=1)
st.pp.normalize_total(adata)
st.pp.log1p(adata)
# pre-processing for spot image
st.pp.tiling(adata, TILE_PATH)
# this step uses deep learning model to extract high-level features from tile images
# may need few minutes to be completed
st.pp.extract_feature(adata)
#2. run stSME clustering
# run PCA for gene expression data
st.em.run_pca(adata,n_comps=50)
adata.write_h5ad(filename='H9_D60_2.h5ad')
#metadata frome seurat analysis
metadata = pd.read_csv('H9_D60_2_final_metadata.csv')
metadata.index = metadata['Unnamed: 0']
metadata = metadata.drop(['Unnamed: 0'],axis=1)
adata.obs = adata.obs.join(metadata, how="left")
adata.obs.New_celltype = adata.obs.New_celltype.fillna("remove")
adata.obs.New_celltype = adata.obs.New_celltype.astype("category")
#Display the images of each cell type separately
celltype = adata.obs.New_celltype.unique()
for i in celltype:
    adata.obs['New_celltype_2'] = adata.obs['New_celltype'].apply(lambda x: i if x == i else 'remove')
    adata.obs['New_celltype_2'] = adata.obs['New_celltype_2'].astype('category')
    adata.uns['New_celltype_2_colors'] = ['#FF0000', '#808080']
    st.pl.cluster_plot(adata, use_label="New_celltype_2",size=13,image_alpha=0.2)
    i = i.replace('/','_')
    plt.savefig('/share/pub/dengcy/STanalysis/H9_D60_2_New_celltype_'+i+'.pdf',dpi=300)

adata.write_h5ad(filename='H9_D60_2.h5ad')
adata = sc.read_h5ad('H9_D60_2.h5ad')
########################################
#H9_D60_2 trajectory
#####################
# Save raw_count
adata = st.Read10X(BASE_PATH)
adata.layers["raw_count"] = adata.X
# Preprocessing
st.pp.filter_genes(adata,min_cells=3)
st.pp.normalize_total(adata)
st.pp.log1p(adata)
# Keep raw data
adata.raw = adata
st.pp.scale(adata)
st.em.run_pca(adata,n_comps=50,random_state=0)
# Tiling image
st.pp.tiling(adata,out_path="tiling",crop_size = 40)
# Using Deep Learning to extract feature
st.pp.extract_feature(adata)
# Apply stSME spatial-PCA option
st.spatial.morphology.adjust(adata,use_data="X_pca",radius=50,method="mean")
st.pp.neighbors(adata,n_neighbors=25,use_rep='X_pca_morphology',random_state=0)
st.tl.clustering.louvain(adata,random_state=0,resolution=2)
st.pl.cluster_plot(adata,use_label="louvain",size=12)
sc.pl.spatial(adata, color="louvain",size=1)

adata.obs.seurat_clusters = adata.obs.seurat_clusters.fillna(7)
adata.obs.seurat_clusters = adata.obs.seurat_clusters.astype(int)
adata.obs.seurat_clusters = adata.obs.seurat_clusters.astype(str)
adata.obs.seurat_clusters = adata.obs.seurat_clusters.astype("category")

adata.uns["iroot"] = st.spatial.trajectory.set_root(adata,use_label="seurat_clusters",cluster="4",use_raw=True)
st.spatial.trajectory.pseudotime(adata,eps=50,use_rep="X_pca",use_label="seurat_clusters")

adata.write_h5ad(filename='H9_D60_2.h5ad')
