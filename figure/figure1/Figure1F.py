##scVelo plot
import anndata
import scvelo as scv
import pandas as pd
import numpy as np
import matplotlib as plt
plt.use('pdf')
import scanpy as sc

##read the result from scVelo.py
adata_merge=sc.read_h5ad("/public/ojsys/eye/sujianzhong/chencheng/Singlecell/LiHui/final/scVelo/scVelo.h5ad")

##embedding_stream
mycolors=["#ECBA84","#A13B46","#CA8C74","#9569AB","#C0BFDF","#E77A77","#7B6148","#6A8473","#71C89C","#67ADB7","#36600E"]
scv.pl.velocity_embedding_stream(adata_merge,basis='X_umap.biMod',color='wsnn_res_0.5_cell_type',frameon=False, dpi=100,save='scVelo_embedding_stream.svg',title='',arrow_size=0.5,linewidth=1, palette=mycolors)
scv.pl.velocity_embedding_stream(adata_merge,basis='X_umap.biMod',color='wsnn_res_0.5_cell_type',frameon=False, dpi=100,save='scVelo_embedding_stream.pdf',title='',arrow_size=0.5,linewidth=1, palette=mycolors)

##embedding
scv.pl.velocity_embedding(adata_merge,basis='X_umap.biMod',color='wsnn_res_0.5_cell_type', arrow_length=3, arrow_size=2, dpi=120,save='scVelo_embedding.pdf',title='',palette=mycolors)

##embedding_grid
scv.pl.velocity_embedding_grid(adata_merge, basis='X_umap.biMod', color='wsnn_res_0.5_cell_type', save='scVelo_embedding_grid.pdf', title='', scale=0.25,palette=mycolors)