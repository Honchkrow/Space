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
# import sys
# seed = int(sys.argv[1])
seed = 0
adata = sc.read("../Data/BARISTASeq/BARISTASeq_Sun2021Integrating_Slice_1_data.h5ad")

adata.var_names_make_unique()
# prefilter_specialgenes(adata)
sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=3000)
# sc.pp.normalize_per_cell(adata)
# sc.pp.log1p(adata)

adata.obs['Ground Truth'] = adata.obs['layer']

sc.pp.pca(adata, n_comps=50)
sc.pp.neighbors(adata, use_rep='X_pca',n_neighbors=50)
sc.tl.leiden(adata, resolution=0.5,random_state=seed) #change(0.5)
sc.tl.umap(adata)

plt.rcParams["figure.figsize"] = (3, 3)
sc.pl.embedding(adata, basis="spatial", color="leiden", s=6, show=False, title='leiden')
plt.axis('off')

Ann_df = adata.obs.dropna()
Ann_df = Ann_df[Ann_df['Ground Truth'] != 'outside_VISp']

ari = metrics.adjusted_rand_score(Ann_df["Ground Truth"], Ann_df["leiden"])
print(ari)

from sklearn.metrics import normalized_mutual_info_score
nmi=normalized_mutual_info_score(Ann_df["Ground Truth"], Ann_df["leiden"])
print('normalized mutual info score = %.5f' % nmi)

plt.rcParams["figure.figsize"] = (3, 3)
sc.pl.scatter(adata, alpha=1, x="x", y="y", color="leiden", legend_fontsize=18, show=True,save="leiden_BARISTASeq.pdf",
                   size=100000 / adata.shape[0])
# sc.pl.scatter(adata, alpha=1, x="x", y="y", color="Ground Truth", legend_fontsize=18, show=True,
#                    size=100000 / adata.shape[0])

# adata.obs.to_csv("/home/zw/stproject/stref/results/Seq/domains_new/slice3_Leiden_{}_Seq_{:.2f}.csv".format(seed,ari))  # change