import stlearn as st
from pathlib import Path
from sklearn import metrics
import pandas as pd
import numpy as np
import scanpy as sc
st.settings.set_figure_params(dpi=180)
import matplotlib.pyplot as plt
import os,csv,re

#传参所用
# import sys
# seed = int(sys.argv[1])
seed = 0
k=8
# adata = sc.read("/home/zw/stproject/Data/Mouse_hippocampus_MERFISH/hip_adata-0.04.h5ad")
adata = sc.read("../Data/BARISTASeq/BARISTASeq_Sun2021Integrating_Slice_1_data.h5ad")

adata.var_names_make_unique()
adata.obs['Ground Truth'] = adata.obs['layer']

adata.obs["x_pixel"] = adata.obs["x_um"]
adata.obs["y_pixel"] = adata.obs["y_um"]
x_pixel=adata.obs["x_pixel"].tolist()
y_pixel=adata.obs["x_pixel"].tolist()

# pre-processing for gene count table
adata.var_names_make_unique()
st.pp.filter_genes(adata, min_cells=3)
st.pp.normalize_total(adata)
st.pp.log1p(adata)
st.pp.scale(adata)

st.em.run_pca(adata, n_comps=15)
st.pp.neighbors(adata, n_neighbors=25)
# st.tl.clustering.louvain(adata)
st.tl.clustering.kmeans(adata, n_clusters=6, use_data="X_pca", key_added="X_pca_kmeans",random_state= seed)

obs_df = adata.obs.dropna()
obs_df = obs_df[obs_df['Ground Truth'] != 'outside_VISp']
ARI = metrics.adjusted_rand_score(obs_df['X_pca_kmeans'], obs_df['Ground Truth'])
print(ARI)

domains='X_pca_kmeans'
title = 'stLearn (ARI=%.2f)' % ARI
ax = sc.pl.scatter(adata, alpha=1, x="x_pixel", y="y_pixel", color=domains, legend_fontsize=18, show=False,
                   size=100000 / adata.shape[0])

ax.set_title(title, fontsize=23)
ax.set_aspect('equal', 'box')
ax.set_xticks([])
ax.set_yticks([])
ax.axes.invert_yaxis()
# plt.savefig("../results/Seq/figures/stLearn_Seq_9.pdf")
plt.close()

# adata.obs.to_csv("../results/Seq/domains/stLearn_9_Seq_%.2f.csv"%ARI)  # change
# adata.obs.to_csv("../results/Seq/domains_new/slice3_stlearn_{}_Seq_{:.2f}.csv".format(seed,ARI))
