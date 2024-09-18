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

adata = sc.read("../Data/BARISTASeq/BARISTASeq_Sun2021Integrating_Slice_1_data.h5ad")

adata.var_names_make_unique()


adata.obs['Ground Truth'] = adata.obs['layer']

sc.pp.pca(adata, n_comps=30)
sc.pp.neighbors(adata, use_rep='X_pca')
sc.tl.louvain(adata, resolution=1.5,random_state=9)   # cahnge(2)
sc.tl.umap(adata)

# tt.to_csv("Data/merfish_{}.csv".format(section_id))

df=adata.obs
df.replace("outside_VISp", np.nan, inplace=True)
# df = df.dropna()
adata.obs[('Ground Truth')]=df['Ground Truth']
adata.obs['louvain']=df['louvain']

Ann_df = adata.obs.dropna()
Ann_df = Ann_df[Ann_df['Ground Truth'] != 'outside_VISp']
# tt = adata.obs
ARI = adjusted_rand_score(Ann_df['louvain'], Ann_df['Ground Truth'])
print('Adjusted rand index = %.4f' % ARI)

domains='louvain'
title = 'SCANPY (ARI=%.2f)' % ARI
plt.rcParams["figure.figsize"] = (3, 3)
ax = sc.pl.scatter(adata, alpha=1, x="x", y="y", color=domains, legend_fontsize=18, show=False,save="scanpy_BARISTASeq.pdf",
                   size=100000 / adata.shape[0])
ax.set_title(title, fontsize=23)
ax.set_aspect('equal', 'box')
ax.set_xticks([])
ax.set_yticks([])
ax.axes.invert_yaxis()
# plt.savefig("../results/Seq/figures/Scanpy_Seq_10.pdf")
plt.close()

from sklearn.metrics import normalized_mutual_info_score
nmi=normalized_mutual_info_score(Ann_df['louvain'], Ann_df['Ground Truth'])
print('normalized mutual info score = %.5f' % nmi)

print(ARI)
# adata.obs.to_csv("../results/Seq/domains/slice3xxxx_scanpy_x_Seq_%.2f.csv"%ARI)      # change
