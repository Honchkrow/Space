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

#传参所用
import sys
# seed = int(sys.argv[1])
seed = 4
# # the location of R (used for the mclust clustering)
os.environ['R_HOME'] = '/home/dell/anaconda3/envs/stagate/lib/R'
# os.environ['R_USER'] = 'D:\ProgramData\Anaconda3\Lib\site-packages\rpy2'

adata = sc.read("../Data/SRARmap_pa/MS_raw_Dataset11_Dataset11_data.h5ad")
adata.var_names_make_unique()
adata = adata[adata.obs['slice_id'] == "BZ5", :]
adata.obs['Ground Truth'] = adata.obs['ct']

#Normalization
sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=3000)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

STAGATE.Cal_Spatial_Net(adata, rad_cutoff=400)
STAGATE.Stats_Spatial_Net(adata)

adata = STAGATE.train_STAGATE(adata, alpha=0, random_seed= seed)                                         ## change

sc.pp.neighbors(adata, use_rep='STAGATE')
sc.tl.umap(adata)
adata = STAGATE.mclust_R(adata, used_obsm='STAGATE', num_cluster=4)

obs_df = adata.obs.dropna()
ARI = adjusted_rand_score(obs_df['mclust'], obs_df['Ground Truth'])
print('Adjusted rand index = %.2f' %ARI)
# adata.obs.to_csv("../results/PA/domains/stagate_9_PA_%.2f.csv"%ARI)                                    # change
plt.rcParams["figure.figsize"] = (3, 3)
sc.pl.umap(adata, color="mclust", title='STAGATE (ARI=%.2f)'%ARI,save="stagate_PA")

plt.rcParams["figure.figsize"] = (3, 3)
sc.pl.embedding(adata, basis="spatial", color="mclust", s=6, show=False, title='stagate',save="stagate_PA.pdf")

# adata.obs.to_csv("../results/PA/domains/slice2_stagate_{}_PA_{:.2f}.csv".format(seed,ARI))