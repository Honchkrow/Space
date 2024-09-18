import anndata as ad
import scanpy as sc
import os,csv,re
import pandas as pd
import numpy as np
import scanpy as sc
import math
import SpaGCN as spg
from scipy.sparse import issparse
import random, torch
import warnings
warnings.filterwarnings("ignore")
from sklearn.metrics import normalized_mutual_info_score,adjusted_rand_score
from sklearn import metrics
import cv2

k=20
# img=cv2.imread("/home/zw/stproject/Data/V1_Breast_Cancer_Block_A_Section_1/V1_Breast_Cancer_Block_A_Section_1_image.tif")
section_id = 'Breast_Cancer'

adata = sc.read_visium(path="../Data/V1_Breast_Cancer_Block_A_Section_1",
                       count_file='filtered_feature_bc_matrix.h5')
adata.var_names_make_unique()

spatial = adata.obsm['spatial']
adata.obs["x_pixel"]=spatial[:,0]
adata.obs["y_pixel"]=spatial[:,1]
adata.obs["x_pixel"]=spatial[:,0]
adata.obs["y_pixel"]=spatial[:,1]

x_pixel=adata.obs["x_pixel"].tolist()
y_pixel=adata.obs["y_pixel"].tolist()
s = 1
b = 49
adj = spg.calculate_adj_matrix(x=x_pixel, y=y_pixel, x_pixel=x_pixel, y_pixel=y_pixel,beta=b, alpha=s,
                               histology=False)

adata.var_names_make_unique()
spg.prefilter_genes(adata, min_cells=3)  # avoiding all genes are zeros
spg.prefilter_specialgenes(adata)
# Normalize and take log for UMI
sc.pp.normalize_per_cell(adata)
sc.pp.log1p(adata)

p = 0.5
# Find the l value given p
l = spg.search_l(p, adj, start=0.01, end=1000, tol=0.01, max_run=100)
r_seed = t_seed = n_seed = 9                                                                          # change
# Seaech for suitable resolution
res = spg.search_res(adata, adj, l, k, start=0.7, step=0.1, tol=5e-3, lr=0.05, max_epochs=20,
                     r_seed=r_seed, t_seed=t_seed, n_seed=n_seed)

clf = spg.SpaGCN()
clf.set_l(l)
# Set seed
random.seed(r_seed)
torch.manual_seed(t_seed)
np.random.seed(n_seed)
# Run
clf.train(adata, adj, init_spa=True, init="louvain", res=res, tol=5e-3, lr=0.05, max_epochs=200)
y_pred, prob = clf.predict()
adata.obs["pred"] = y_pred
adata.obs["pred"] = adata.obs["pred"].astype('category')
Ann_df =  pd.read_csv("/home/zw/stproject/Data/V1_Breast_Cancer_Block_A_Section_1/metadata.tsv", sep="	", header=0, na_filter=False,
                          index_col=0)
adata.obs['Ground Truth'] = Ann_df.loc[adata.obs_names, 'fine_annot_type']

# adata.obs["refined_pred"] = spg.spatial_domains_refinement_ez_mode(sample_id=adata.obs.index.tolist(),
#                                                                    pred=adata.obs["pred"].tolist(), x_array=x_array,
#                                                                    y_array=y_array, shape="hexagon")
obs_df = adata.obs.dropna()
ari1=metrics.adjusted_rand_score(obs_df["pred"] ,obs_df["Ground Truth"])
print(ari1)

# adata.obs["refined_pred"]=spg.spatial_domains_refinement_ez_mode(sample_id=adata.obs.index.tolist(), pred=adata.obs["pred"].tolist(), x_array=x_array, y_array=y_array, shape="hexagon")
# ari2 = metrics.adjusted_rand_score(obs_df["refined_pred"] ,obs_df["Ground Truth"] )
# print("ari2", ari2)


nmi=normalized_mutual_info_score(obs_df["pred"] ,obs_df["Ground Truth"] )
print(nmi)
sc.pl.spatial(adata, color="pred", title='SpaGCN (ARI=%.2f)' % ari1,
           save="SpaGCbreast.pdf")

obs_df = adata.obs.dropna()
#  ARI = adjusted_rand_score(obs_df['stGCL_refined'], obs_df['Ground Truth'])         # GCL?               # error!
#  print('Adjusted rand index = %.2f' %ARI)                                                                # error!

print('%.5f ---------------------------------------------------------' %ari1)
# adata.obs.to_csv("../results/breast/domains/spagcn_9_breast_%.2f.csv"%ari1)                                 # change