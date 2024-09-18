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
import sys
# seed = int(sys.argv[1])

seed = 0
k=8
adata = sc.read("../Data/Mouse_hippocampus_MERFISH/hip_adata-0.04.h5ad")

adata.var_names_make_unique()
spatial = adata.obsm['spatial']
adata.obs["x_pixel"]=spatial[:,0]
adata.obs["y_pixel"]=spatial[:,1]
adata.obs["x_pixel"]=spatial[:,0]
adata.obs["y_pixel"]=spatial[:,1]
adata.obs['Ground Truth'] = adata.obs['ground_truth']
x_pixel=adata.obs["x_pixel"].tolist()
y_pixel=adata.obs["y_pixel"].tolist()

# pre-processing for gene count table
adata.var_names_make_unique()
st.pp.filter_genes(adata, min_cells=3)
st.pp.normalize_total(adata)
st.pp.log1p(adata)
st.pp.scale(adata)

st.em.run_pca(adata, n_comps=15)
st.pp.neighbors(adata, n_neighbors=25)      #25
# st.tl.clustering.louvain(adata)
st.tl.clustering.kmeans(adata, n_clusters=8, use_data="X_pca", key_added="X_pca_kmeans",random_state= seed)


df=adata.obs
df.replace("outside_VISp", np.nan, inplace=True)
df = df.dropna()
adata.obs['X_pca_kmeans']=df['X_pca_kmeans']
obs_df = adata.obs.dropna()
ARI = metrics.adjusted_rand_score(obs_df['X_pca_kmeans'], obs_df['Ground Truth'])
print(ARI)

domains='X_pca_kmeans'
sc.pl.embedding(adata, basis="spatial", color="X_pca_kmeans", s=6, show=False, title='leiden',save="stlearn_hip.pdf")

print(ARI,"---------------------------------")
# adata.obs.to_csv("../results/HIP/domains/stLearn_9_HIP_%.2f.csv"%ARI)  # change
# adata.obs.to_csv("../results/HIP/domains/hip024_stlearn_{}_hip_{:.2f}.csv".format(seed,ARI))
