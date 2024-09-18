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
# import sys
# seed = int(sys.argv[1])
seed = 1

adata = sc.read("../Data/SRARmap_pa/MS_raw_Dataset11_Dataset11_data.h5ad")
adata.var_names_make_unique()
adata = adata[adata.obs['slice_id'] == "BZ5", :]
adata.obs['Ground Truth'] = adata.obs['ct']
# prefilter_specialgenes(adata)
sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=3000)
sc.pp.normalize_per_cell(adata)
sc.pp.log1p(adata)

sc.pp.pca(adata, n_comps=30)     #30
sc.pp.neighbors(adata, use_rep='X_pca')
sc.tl.louvain(adata, resolution=1.5, random_state= seed) #2
sc.tl.umap(adata)



Ann_df = adata.obs.dropna()

# tt = adata.obs
ARI = adjusted_rand_score(Ann_df['louvain'], Ann_df['Ground Truth'])
print('Adjusted rand index = %.4f' % ARI)

domains='louvain'
title = 'SCANPY (ARI=%.2f)' % ARI
plt.rcParams["figure.figsize"] = (3, 3)
sc.pl.embedding(adata, basis="spatial", color="louvain", s=6, show=False, title='scanpy',save="scanpy_pa.pdf")

from sklearn.metrics import normalized_mutual_info_score
nmi=normalized_mutual_info_score(Ann_df['louvain'], Ann_df['Ground Truth'])
print('normalized mutual info score = %.5f' % nmi)

print(ARI,"=========================")
# adata.obs.to_csv("../results/PA/domains/scanpy_10_pa_%.2f.csv"%ARI)  # change
# adata.obs.to_csv("../results/PA/domains/slice2_scanpy_{}_PA_{:.2f}.csv".format(seed,ARI))
