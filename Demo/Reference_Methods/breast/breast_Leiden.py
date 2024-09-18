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


section_id = "V1_Breast_Cancer_Block_A_Section_1"
k = 20
print(section_id,k)
# input_dir = "/home/zw/stproject/Data/V1_Breast_Cancer_Block_A_Section_1"
# adata = sc.read_visium(path=input_dir, count_file='filtered_feature_bc_matrix.h5')
adata = sc.read_visium(path="../Data/V1_Breast_Cancer_Block_A_Section_1",
                       count_file='filtered_feature_bc_matrix.h5')
adata.var_names_make_unique()

sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=3000)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

sc.pp.pca(adata, n_comps=30)
sc.pp.neighbors(adata, use_rep='X_pca')
# sc.tl.louvain(adata, resolution=0.5)
# sc.tl.umap(adata)

sc.tl.leiden(adata, resolution=0.5 ,random_state=0)                                                          # change
sc.tl.umap(adata)

plt.rcParams["figure.figsize"] = (3, 3)
sc.pl.embedding(adata, basis="spatial", color="leiden", s=6, show=False, title='SCANPY')
plt.axis('off')

Ann_df = pd.read_csv("../Data/V1_Breast_Cancer_Block_A_Section_1/metadata.tsv", sep="	", header=0,
                     na_filter=False,index_col=0)
Ann_df["leiden"]=adata.obs['leiden']
Ann_df.dropna(inplace=True)

ari = metrics.adjusted_rand_score(Ann_df["fine_annot_type"], Ann_df["leiden"])
print(ari)

from sklearn.metrics import normalized_mutual_info_score
nmi=normalized_mutual_info_score(Ann_df["fine_annot_type"], Ann_df["leiden"])
print('normalized mutual info score = %.5f' % nmi)

plt.rcParams["figure.figsize"] = (3, 3)
sc.pl.spatial(adata, color="leiden", title='Leiden (ARI=%.2f)' % ari,
           save="scanpy")

# adata.obs.to_csv("../results/breast/domains/leiden_9_breast_%.2f.csv"%ari)                                  # change