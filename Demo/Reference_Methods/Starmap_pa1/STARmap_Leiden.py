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
adata = sc.read("../Data/SRARmap_pa/MS_raw_Dataset11_Dataset11_data.h5ad")
adata.var_names_make_unique()
adata = adata[adata.obs['slice_id'] == "BZ5", :]
adata.obs['Ground Truth'] = adata.obs['ct']
# prefilter_specialgenes(adata)
sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=3000)
sc.pp.normalize_per_cell(adata)
sc.pp.log1p(adata)
sc.pp.pca(adata, n_comps=40)    #30
sc.pp.neighbors(adata, use_rep='X_pca')
sc.tl.leiden(adata, resolution=1.2, random_state= seed) #0.5
sc.tl.umap(adata)

plt.rcParams["figure.figsize"] = (3, 3)
sc.pl.embedding(adata, basis="spatial", color="leiden", s=6, show=False, title='leiden',save="leiden_pa.pdf")
plt.axis('off')

Ann_df = adata.obs.dropna()

ari = metrics.adjusted_rand_score(Ann_df["Ground Truth"], Ann_df["leiden"])
print(ari)

from sklearn.metrics import normalized_mutual_info_score
nmi=normalized_mutual_info_score(Ann_df["Ground Truth"], Ann_df["leiden"])
print('normalized mutual info score = %.5f' % nmi)

# adata.obs.to_csv("../results/PA/domains/silce2_leiden_{}_PA_{:.2f}.csv".format(seed,ari))
# adata.obs.to_csv("../results/PA/domains/Leiden_9_PA_%.2f.csv"%ari)  # change