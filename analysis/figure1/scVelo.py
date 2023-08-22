##run scVelo
import anndata
import scvelo as scv
import pandas as pd
import numpy as np
import matplotlib as plt
plt.use('pdf')
import scanpy as sc

##input the h5ad from scVelo_input.R
adata_h5ad=sc.read_h5ad("/public/ojsys/eye/sujianzhong/chencheng/Singlecell/LiHui/final/scVelo/data.h5ad")
adata_h5ad.obs["wsnn_res_0.5_cell_type"]

##input metadata from scVelo_input.R
meta=pd.read_csv("/public/ojsys/eye/sujianzhong/chencheng/Singlecell/LiHui/final/scVelo/metadata.csv",index_col=0)
adata_h5ad.obs=meta

##input the loom
adata_loom = anndata.read_loom("/public/ojsys/eye/sujianzhong/chencheng/Singlecell/LiHui/final_version/scVelo/total_fetus.loom")
adata_loom.var_names_make_unique()

##change the cell barcodes in loom
barcodes=adata_loom.obs.index.tolist()
for i in range(len(barcodes)):
    if barcodes[i].startswith('NR'):
        barcodes[i]=barcodes[i].replace('NR:','').replace('x','-1_1')
    elif barcodes[i].startswith('CMZ'):
        barcodes[i]=barcodes[i].replace('CMZ:','').replace('x','-1_2')
adata_loom.obs.index=barcodes

##subset the cells both in adata_h5ad and adata_loom
cell_names=adata_h5ad.obs.index 
adata_loom_subset=adata_loom[adata_loom.obs_names.isin(cell_names), :]

##merge adata_h5ad and adata_loom
adata_merge=scv.utils.merge(adata_h5ad,adata_loom_subset)
#sc.write('merged.h5ad',adata_merge,compression='gzip',compression_opts=1)

##figure setting
scv.settings.verbosity = 3
scv.settings.set_figure_params('scvelo', facecolor='white', dpi=100, frameon=False)
scv.settings.rcParams['font.size']=4  ##font.size

##data preprocess
scv.pp.filter_and_normalize(adata_merge, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(adata_merge, n_pcs=30, n_neighbors=30)

##RNA Velocity
scv.tl.recover_dynamics(adata_merge)
scv.tl.velocity(adata_merge, mode='dynamical')
#scv.tl.velocity(adata_merge, mode='stochastic')
scv.tl.velocity_graph(adata_merge)

##save the result
adata_merge.write('scVelo.h5ad')
