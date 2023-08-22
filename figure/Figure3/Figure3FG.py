import stlearn as st
import scanpy as sc
import matplotlib.pyplot as plt
import scanpy as sc
adata = sc.read_h5ad('H9_D60_2.h5ad')
st.pl.trajectory.pseudotime_plot(adata, use_label="seurat_clusters", pseudotime_key="dpt_pseudotime", list_clusters=['4'], show_node=True,  image_alpha=0.2,spot_size=17,edge_alpha=0,node_alpha=0,margin=50)
plt.savefig('/share/pub/dengcy/STanalysis/Figure3G_merge.pdf',dpi=300)
#gene_symbols is the marker genes for hRSLCs cells
data2 = st.convert_scanpy(adata)
st.pl.gene_plot(adata2, gene_symbols=['RMST','PCDH9','PCDH7','MECOM','RELN','COL9A1','CPAMD8','ZIC1','SYNPR','UNC5C','RP11-141M1.3','FGF19'], method="CumSum",size=12,image_alpha=0.2,title='H9_D60_2: hRSLCs marker genes')
plt.savefig('Figure3F.pdf',dpi=300)


#################
#Perform image segmentation on adata based on spatial coordinates, and divide it into 4 regions for calculation.
adata.obs['array_row'].describe()
#min        5.000000
#25%       18.000000
#50%       35.500000
#75%       43.000000
#max       67.000000
adata.obs['array_col'].describe()
#min       28.000000
#25%       46.000000
#50%       70.000000
#75%       98.000000
#max      118.000000
#Top left corner
adata_1 = adata[adata.obs['array_row'] < 26,:]
adata_1 = adata_1[adata_1.obs['array_col'] < 59,:]
#Check whether the partition is successful
adata_1.obs['array_row'].describe()
#array_row&lt; 12 and 45 points off
test = adata_1[(adata_1.obs['array_row'] < 12) & (adata_1.obs['array_col'] > 50),:]
adata_1 = adata_1[~adata_1.obs.index.isin(test.obs.index),:]

st.pl.cluster_plot(adata_1,use_label="seurat_clusters",size=40)
plt.savefig('/share/pub/dengcy/STanalysis/Figure3G_split1.pdf',dpi=300)
adata_1.obs.seurat_clusters = adata_1.obs.seurat_clusters.astype('category')
adata_1.uns["iroot"] = st.spatial.trajectory.set_root(adata_1,use_label="seurat_clusters",cluster="3",use_raw=True)
st.spatial.trajectory.pseudotime(adata_1,eps=50,use_rep="X_pca",use_label="seurat_clusters")
st.spatial.trajectory.pseudotimespace_global(adata_1,use_label="seurat_clusters",list_clusters=["5","4"])
st.pl.cluster_plot(adata_1,use_label="seurat_clusters",show_trajectories=True,list_clusters=["5","4"],show_subcluster=True)
plt.savefig('/share/pub/dengcy/STanalysis/Figure3G_4_5_sub1.pdf',dpi=300)