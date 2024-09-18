import warnings
warnings.filterwarnings("ignore")
import STAGATE
import pandas as pd
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import os
import sys
import tensorflow as tf
tf.compat.v1.disable_eager_execution()
from sklearn.metrics.cluster import adjusted_rand_score

# # the location of R (used for the mclust clustering)
os.environ['R_HOME'] = '/home/dell/anaconda3/envs/stagate/lib/R'
# os.environ['R_USER'] = 'D:\ProgramData\Anaconda3\Lib\site-packages\rpy2'

section_id = 'Breast_Cancer'

adata = sc.read_visium(path="../Data/V1_Breast_Cancer_Block_A_Section_1",
                       count_file='filtered_feature_bc_matrix.h5')
adata.var_names_make_unique()

#Normalization
sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=3000)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

Ann_df =  pd.read_csv("../Data/V1_Breast_Cancer_Block_A_Section_1/metadata.tsv", sep="	", header=0, na_filter=False,
                          index_col=0)
adata.obs['Ground Truth'] = Ann_df.loc[adata.obs_names, 'fine_annot_type']

plt.rcParams["figure.figsize"] = (3, 3)
sc.pl.spatial(adata, img_key="hires", color=["Ground Truth"])

STAGATE.Cal_Spatial_Net(adata, rad_cutoff=400)
STAGATE.Stats_Spatial_Net(adata)

adata = STAGATE.train_STAGATE(adata, alpha=0, random_seed=5)                                         ## change



sc.pp.neighbors(adata, use_rep='STAGATE')
sc.tl.umap(adata)
adata = STAGATE.mclust_R(adata, used_obsm='STAGATE', num_cluster=20)

obs_df = adata.obs.dropna()
ARI = adjusted_rand_score(obs_df['mclust'], obs_df['Ground Truth'])
print('Adjusted rand index = %.2f' %ARI)
# adata.obs.to_csv("../results/breast/domains/stagate_5_breast_%.2f.csv"%ARI)                                    # change
plt.rcParams["figure.figsize"] = (3, 3)
sc.pl.umap(adata, color="mclust", title='STAGATE (ARI=%.2f)'%ARI,save="stagate_breast")

plt.rcParams["figure.figsize"] = (3, 3)
sc.pl.spatial(adata, color="mclust", title='STAGATE (ARI=%.2f)'%ARI,save="stagate_breast")

from sklearn.metrics import normalized_mutual_info_score
nmi=normalized_mutual_info_score(obs_df['mclust'], obs_df['Ground Truth'])
print('normalized mutual info score = %.5f' % nmi)