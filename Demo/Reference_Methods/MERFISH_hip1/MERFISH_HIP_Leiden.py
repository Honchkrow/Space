import warnings
warnings.filterwarnings("ignore")

import pandas as pd
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import os
import sys
from sklearn.metrics.cluster import adjusted_rand_score
from datetime import datetime
from sklearn import metrics

#传参所用
import sys
# seed = int(sys.argv[1])
seed = 0
adata = sc.read("../Data/Mouse_hippocampus_MERFISH/hip_adata-0.04.h5ad")

adata.var_names_make_unique()
# prefilter_specialgenes(adata)
sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=3000)
sc.pp.normalize_per_cell(adata)
sc.pp.log1p(adata)

adata.obs['Ground Truth'] = adata.obs['ground_truth']
sc.pp.pca(adata, n_comps=30)            # 30
sc.pp.neighbors(adata, use_rep='X_pca')
sc.tl.leiden(adata, resolution=0.5 ,random_state= seed)     # 0.5
sc.tl.umap(adata)

plt.rcParams["figure.figsize"] = (3, 3)
sc.pl.embedding(adata, basis="spatial", color="leiden", s=6, show=False, title='leiden',save="leiden_hip.pdf")
plt.axis('off')

Ann_df = adata.obs.dropna()

ari = metrics.adjusted_rand_score(Ann_df["Ground Truth"], Ann_df["leiden"])
print(ari)

from sklearn.metrics import normalized_mutual_info_score
nmi=normalized_mutual_info_score(Ann_df["Ground Truth"], Ann_df["leiden"])
print('normalized mutual info score = %.5f' % nmi)

print(ari,"-------------------------------------")
# adata.obs.to_csv("../results/HIP/domains/Leiden_9_HIP_%.2f.csv"%ari)  # change
# adata.obs.to_csv("../results/HIP/domains/hip024_leiden_{}_hip_{:.2f}.csv".format(seed,ari))