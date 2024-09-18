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
# import sys
# seed = int(sys.argv[1])
seed = 4
# # the location of R (used for the mclust clustering)
os.environ['R_HOME'] = '/home/dell/anaconda3/envs/stagate/lib/R'
# os.environ['R_USER'] = 'D:\ProgramData\Anaconda3\Lib\site-packages\rpy2'

adata = sc.read("../Data/BARISTASeq/BARISTASeq_Sun2021Integrating_Slice_1_data.h5ad")
adata.var_names_make_unique()

#Normalization
sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=3000)
sc.pp.normalize_total(adata, target_sum=1e4)
# sc.pp.log1p(adata)

adata.obs['Ground Truth'] = adata.obs['layer']


STAGATE.Cal_Spatial_Net(adata, rad_cutoff=50)
STAGATE.Stats_Spatial_Net(adata)

adata = STAGATE.train_STAGATE(adata,  alpha=0, random_seed=seed)                                         ## change



sc.pp.neighbors(adata, use_rep='STAGATE')
sc.tl.umap(adata)
adata = STAGATE.mclust_R(adata, used_obsm='STAGATE', num_cluster=6)

obs_df = adata.obs.dropna()
obs_df = obs_df[obs_df['Ground Truth'] != 'outside_VISp']
ARI = adjusted_rand_score(obs_df['mclust'], obs_df['Ground Truth'])
print('Adjusted rand index = %.2f' %ARI)
# adata.obs.to_csv("../results/Seq/domains/slice1_stagate_0_Seq_%.2f.csv"%ARI)                                    # change
obs_df = obs_df[obs_df['Ground Truth'] != 'outside_VISp']
print(ARI)
# adata.obs.to_csv("../results/Seq/domains_new/slice3_stagate_{}_Seq_{:.2f}.csv".format(seed,ARI))