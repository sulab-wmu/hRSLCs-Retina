import numpy as np
import scanpy as sc
import pandas as pd
import os
import squidpy as sq
import matplotlib.pyplot as plt
import seaborn as sns
sc.logging.print_versions()
sc.set_figure_params(facecolor="white", figsize=(8, 8))
sc.settings.verbosity = 3


adata=sc.read_h5ad('/share/pub/dengcy/STanalysis/fetusdata/report/upload/1_SAW_analysis/2_expression/W19-fetal.spatial.cluster.h5ad')

resolution=0.5
#sc.pp.neighbors(adata, n_neighbors=30, n_pcs=30)
#sc.tl.umap(adata)
sc.tl.louvain(adata, resolution=2)

os.listdir("/share/pub/dengcy/STanalysis/fetusdata/analysis/")
os.chdir("/share/pub/dengcy/STanalysis/fetusdata/analysis/")
sc.pl.umap(adata, color='leiden', legend_fontsize=12, legend_fontoutline=2,frameon=False,title='clustering of spots', palette='Set1',save="fetal_leiden_scatter")

sc.pl.umap(adata, color='louvain', legend_fontsize=12, legend_fontoutline=2,frameon=False,title='clustering of spots', palette='Set1',save="fetal_louvain_scatter")

import matplotlib.pyplot as plt

# 假设你的AnnData对象名为adata

# 自定义离散颜色列表
#custom_colors = ["#71C89C","#67ADB7","#36600E", "#6A8473", "#C0BFDF","#E77A77","#7B6148","#ECBA84","#CA8C74","#A13B46","#9569AB","#71C89C","#67ADB7","#36600E", "#6A8473", "#C0BFDF","#E77A77","#7B6148","#ECBA84","#CA8C74","#A76F6F","#A13B46","#9569AB","#AAC8A7","#FFD6A5","#967E76","#9D5353","#C37B89","#5B6D5B","#DE8971"]

adata.obs['louvain2'] = adata.obs['louvain']
sc.pl.embedding(adata,basis='spatial',color ='louvain2',save='Spatial_louvain_2_plot',title="resulotion=2")

sc.tl.louvain(adata, resolution=3)
adata.obs['louvain3'] = adata.obs['louvain']
sc.pl.embedding(adata,basis='spatial',color ='louvain3',save='Spatial_louvain_3_plot',title="resulotion=3")

sc.tl.louvain(adata, resolution=4)
adata.obs['louvain4'] = adata.obs['louvain']
sc.pl.embedding(adata,basis='spatial',color ='louvain4',save='Spatial_louvain_4_plot',title="resulotion=4")

sc.tl.score_genes(adata, gene_list=['RMST','PCDH9','PCDH7','MECOM','RELN','COL9A1','CPAMD8','ZIC1','SYNPR','UNC5C','RP11-141M1.3','FGF19'],score_name='fetal_cmz_score')

import pandas as pd
# 选择大于10的行，并将它们的值设置为10
adata.obs['fetal_cmz_score'] = adata.obs['fetal_cmz_score'].apply(lambda x: 10 if x > 10 else x)

sc.pl.embedding(adata,basis='spatial',color ='fetal_cmz_score',save='Spatial_fetal_cmz_score_plot')


###############
#绘制聚类图
###############
method="wilcoxon"
sc.tl.rank_genes_groups(adata, groupby='louvain4', method="wilcoxon", key_added="dea_1_louvain4")
sc.tl.filter_rank_genes_groups(adata,min_in_group_fraction=0.1,max_out_group_fraction=0.1, key="dea_1_louvain4",key_added="dea_1_louvain4_filtered")
n_genes=500
result = adata.uns['dea_1_louvain4']
groups = result['names'].dtype.names
pval_table = pd.DataFrame({group + '_' + key[:2]: result[key][group] for group in groups for key in ['names', 'pvals_adj']}).head(n_genes)
pval_table.to_excel(os.path.join('_pval_table_' + method + '_louvain4_' + str(n_genes) + '_annotation.xlsx'), engine='openpyxl')
del adata.uns['dea_1_leiden_filtered']
adata.write("adata.h5ad")
#输出结果
n_genes=500
result = adata.uns['dea_1_louvain4_filtered']
groups = result['names'].dtype.names
pval_table = pd.DataFrame({group + '_' + key[:2]: result[key][group] for group in groups for key in ['names', 'pvals_adj']}).head(n_genes)
pval_table.to_excel(os.path.join('_pval_table_' + method + '_louvain4_' + str(n_genes) + '_filter_annotation.xlsx'), engine='openpyxl')

sc.tl.dendrogram(adata,groupby="louvain")
sc.pl.rank_genes_groups_dotplot(adata,groupby="louvain",standard_scale="var",n_genes=5,key="dea_1_louvain_filtered",save='louvain_filter_gene_dotplot')
obs=adata.obs
obs['spotnames']= "W19-fetal:" + adata.obs_names
obs.to_csv("/share/pub/dengcy/STanalysis/fetusdata/analysis/obs.csv")

####################
#分割图形数据
#################
os.chdir("/share/pub/dengcy/STanalysis/fetusdata/analysis/")
adata=sc.read_h5ad('/share/pub/dengcy/STanalysis/fetusdata/analysis/adata.h5ad')
custom_colors = ["#71C89C","#67ADB7","#36600E", "#6A8473", "#C0BFDF","#E77A77","#7B6148","#ECBA84","#CA8C74","#A13B46","#9569AB","#71C89C","#67ADB7","#36600E", "#6A8473", "#C0BFDF","#E77A77","#7B6148","#ECBA84","#CA8C74","#A76F6F","#A13B46","#9569AB","#AAC8A7","#FFD6A5","#967E76","#9D5353","#C37B89","#5B6D5B","#DE8971"]

adata.obs['x'].describe()
adata.obs['y'].describe()
adata1 = adata[adata.obs['x'] <5000,:]
adata1 = adata1[adata1.obs['y'] < 11000,:]
adata1 = adata1[adata1.obs['y'] > 3000,:]
sc.set_figure_params(facecolor="white", figsize=(4, 6))

sc.pl.embedding(adata1,basis='spatial',color ='louvain4',save='Sub_Spatial_rawlouvain4_plot')

sc.tl.louvain(adata1, resolution=4)
#adata.obs['louvain3'] = adata.obs['louvain']
sc.pl.embedding(adata1,basis='spatial',color ='louvain',save='sub_Spatial_louvain_4_plot',title="resulotion=4")

sc.tl.score_genes(adata1, gene_list=['RMST','PCDH9','PCDH7','MECOM','RELN','COL9A1','CPAMD8','ZIC1','SYNPR','UNC5C','RP11-141M1.3','FGF19'],score_name='fetal_cmz_score')

import pandas as pd

# 选择大于10的行，并将它们的值设置为10
#adata1.obs['fetal_cmz_score'] = adata1.obs['fetal_cmz_score'].apply(lambda x: 6 if x > 6 else x)

sc.pl.embedding(adata1,basis='spatial',color ='fetal_cmz_score',save='Sub_Spatial_fetal_cmz_score_plot')

sc.pl.embedding(adata1,basis='spatial',color =['PCDH9','PCDH7','MECOM','RELN','COL9A1','CPAMD8','ZIC1','SYNPR','UNC5C','FGF19'],save='Sub_Spatial_fetal_cmz_gene_plot')
#############
#差异基因计算
###########

method="wilcoxon"
sc.tl.rank_genes_groups(adata1, groupby='louvain', method="wilcoxon", key_added="dea_1_louvain")
sc.tl.filter_rank_genes_groups(adata1,min_in_group_fraction=0.1,max_out_group_fraction=0.1, key="dea_1_louvain",key_added="dea_1_louvain_filtered")
n_genes=100
result = adata1.uns['dea_1_louvain']
groups = result['names'].dtype.names
pval_table = pd.DataFrame({group + '_' + key[:2]: result[key][group] for group in groups for key in ['names', 'pvals_adj']}).head(n_genes)
pval_table.to_excel(os.path.join('_subspatial_pval_table_' + method + '_louvain_' + str(n_genes) + '_annotation.xlsx'), engine='openpyxl')


sc.tl.dendrogram(adata1,groupby="louvain")
sc.pl.rank_genes_groups_dotplot(adata1,groupby="louvain",standard_scale="var",n_genes=5,key="dea_1_louvain_filtered",save='subspatial_louvain_filter_gene_dotplot')


####################################################
#新的空间数据的分析
##################################

import anndata
import pyarrow.feather as feather
#import stlearn as st
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as pl
from matplotlib import rcParams
sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_versions()
results_file = '/share/pub/dengcy/STanalysis/fetusdata/analysis/fetal_data_adata.h5ad'
sc.settings.set_figure_params(dpi=80, frameon=False, figsize=(3, 3), facecolor='white')  # low dpi (dots per inch) yields small inline figures

from pathlib import Path
#import matplotlib.pyplot as plt
os.chdir("/share/pub/dengcy/STanalysis/fetusdata/analysis/")
sct_matrix = feather.read_feather("sct_matrix.feather")
sct_matrix = sct_matrix.T
# 读取SCT矩阵文件
metadata = pd.read_csv("metadata.csv", index_col=0)
var_names = pd.read_csv("genenames.csv", index_col=0)
adata = anndata.AnnData(X=sct_matrix, obs=metadata)
adata.var_names = var_names['rownames.fetal_data.assays.SCT.counts.']
# 打印AnnData对象的摘要信息
print(adata)
sc.set_figure_params(fontsize=14, dpi=80, dpi_save=300, format='svg')
adata.X = adata.X.astype('float64')
sc.pp.recipe_zheng17(adata,n_top_genes=5000)
sc.tl.pca(adata, svd_solver='arpack')
sc.pp.neighbors(adata, n_neighbors=4, n_pcs=20)

sc.tl.louvain(adata, resolution=2.0)
sc.tl.paga(adata, groups='louvain')
sc.pl.paga(adata, color=['louvain'],save='fetal_paga_cluster',fontsize =6,fontoutline =0.1,node_size_scale=0.5,edge_width_scale=0.5)
sc.pl.paga(adata, threshold=0.05,fontsize =6,fontoutline =0.1,node_size_scale=0.5,edge_width_scale=0.5, show=False,save='_louvain_paga_filter')


gene_marker=pd.read_csv("singlecell_fetal_genemarkers.csv")
genes1 = gene_marker[gene_marker['CELL TYPLE'] == 'CMZ_RSCs'][' GENES'].values
genes2 = gene_marker[gene_marker['CELL TYPLE'] == 'RPR_progenitors'][' GENES'].values

genes1=set(genes1)
genes2=set(genes2)

var_names=set(adata.var_names)
genes1=var_names.intersection(genes1)
genes2=var_names.intersection(genes2)

sc.tl.draw_graph(adata, init_pos='paga')
genes=['OPTC','AOC2','COL9A1','RELN','ST13','CHSY3','RAX','GPC6','COL9A2','GPX3','NBL1','COL9A3','SFRP2','SERPINF1','C6orf48','CPAMD8','CIRBP']
genes=['MECOM']
genes3=var_names.intersection(genes)
#sc.pl.draw_graph(adata, color=['louvain','RELN', 'ZIC1'], legend_loc='on data',save='_louvain_CMZ_RSCs_marker')
sc.pl.dotplot(adata, var_names=['SFRP2', 'SERPINF1', 'RAX', 'RELN'], groupby='louvain',save='DotPlot_louvain_marker')

zeileis_colors = np.array(sc.pl.palettes.zeileis_28)
new_colors = np.array(adata.uns['louvain_colors'])
adata.uns['louvain_colors'] = new_colors

adata.uns['iroot'] = np.flatnonzero(adata.obs['louvain']  == '18')[0]
sc.tl.dpt(adata)
sc.pl.draw_graph(adata, color=['louvain', 'dpt_pseudotime'], legend_loc='on data',save='dpt_pseudotime_louvain')
obs=adata.obs
adata.obs.to_csv('obs_dpt_pseudotime_louvain.csv')




sc.pp.normalize_total(adata, inplace=True)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=6000)
adata.write("fetal_data_adata.h5ad")
sc.tl.pca(adata, svd_solver='arpack',n_comps=50)
#sc.pp.neighbors(adata, n_neighbors=10, n_pcs=50)
sc.tl.umap(adata)
adata.obs['seurat_clusters'] = adata.obs['seurat_clusters'].astype('category')
sc.pl.umap(adata,color=['seurat_clusters'],save="python_umap_plot")



