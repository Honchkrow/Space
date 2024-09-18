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
import matplotlib.pyplot as plt

#传参所用
# import sys
# seed = int(sys.argv[1])
seed = 666
k=8
adata = sc.read("../Data/SRARmap_pa/MS_raw_Dataset11_Dataset11_data.h5ad")
adata.var_names_make_unique()
adata = adata[adata.obs['slice_id'] == "BZ5", :]
adata.obs['Ground Truth'] = adata.obs['ct']
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
r_seed = t_seed = n_seed = seed                                                                          # change
# Seaech for suitable resolution
res = 0.5

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
plt.rcParams["figure.figsize"] = (3, 3)
sc.pl.embedding(adata, basis="spatial", color="pred", s=6, show=False, title='spagcn',save="spagcn_pa.pdf")

obs_df = adata.obs.dropna()
#  ARI = adjusted_rand_score(obs_df['stGCL_refined'], obs_df['Ground Truth'])         # GCL?               # error!
#  print('Adjusted rand index = %.2f' %ARI)                                                                # error!

print('%.5f ---------------------------------------------------------' %ari1)
# adata.obs.to_csv("../results/PA/domains/spagcn_x_PA_%.2f.csv"%ari1)                                 # change
# adata.obs.to_csv("../results/PA/domains/slice2_spagcn_{}_PA_{:.2f}.csv".format(seed,ari1))