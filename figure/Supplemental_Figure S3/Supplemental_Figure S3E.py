import numpy as np
import pandas as pd
import scanpy as sc
from matplotlib import rcParams
import seaborn as sns 

sc.settings.verbosity = 3             
sc.logging.print_versions()
sc.settings.set_figure_params(dpi=100,frameon=False)

adata = sc.read_h5ad("/public/ojsys/eye/sujianzhong/chencheng/Singlecell/LiHui/final/scVelo/data.h5ad")
adata.obs["wsnn_res_0.5_cell_type"]

meta=pd.read_csv("/public/ojsys/eye/sujianzhong/chencheng/Singlecell/LiHui/final/scVelo/metadata.csv",index_col=0)
adata.obs=meta

sc.pp.neighbors(adata)

sc.tl.paga(adata, groups="wsnn_res_0.5_cell_type")

adata.uns['wsnn_res_0.5_cell_type_colors']=["#ECBA84","#A13B46","#CA8C74","#9569AB","#C0BFDF","#E77A77","#7B6148","#6A8473","#71C89C","#67ADB7","#36600E"]
sc.pl.paga(adata, threshold=0.05,color="wsnn_res_0.5_cell_type",fontsize=10,edge_width_scale=0.5,frameon=False,show=True,save=True)