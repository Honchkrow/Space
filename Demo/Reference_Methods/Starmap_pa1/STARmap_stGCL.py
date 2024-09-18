import warnings
warnings.filterwarnings("ignore")
import pandas as pd
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import os
import torch
from sklearn.metrics.cluster import adjusted_rand_score
import time
import stGCL as stGCL
from stGCL.process import prefilter_genes,prefilter_specialgenes,set_seed,refine_nearest_labels
from stGCL import utils,train_model
import ot
# # the location of R (used for the mclust clustering)
os.environ['R_HOME'] = '/home/dell/anaconda3/envs/stagate/lib/R'
# os.environ['R_USER'] = 'D:\ProgramData\Anaconda3\Lib\site-packages\rpy2'

#传参所用
import sys
# seed = int(sys.argv[1])
seed = 9
use_image=False
radius=50
set_seed(seed)
radius = 70
top_genes = 3000
epoch = 800
def refine_label(adata, radius=50, key='cluster'):
    n_neigh = radius
    new_type = []
    old_type = adata.obs[key].values

    # calculate distance
    position = adata.obsm['spatial']
    distance = ot.dist(position, position, metric='euclidean')
    n_cell = distance.shape[0]

    for i in range(n_cell):
        vec = distance[i, :]
        index = vec.argsort()
        neigh_type = []
        for j in range(1, n_neigh + 1):
            neigh_type.append(old_type[index[j]])
        max_type = max(neigh_type, key=neigh_type.count)
        new_type.append(max_type)

    new_type = [str(i) for i in list(new_type)]
    # adata.obs['label_refined'] = np.array(new_type)

    return new_type
adata = sc.read("../Data/SRARmap_pa/MS_raw_Dataset11_Dataset11_data.h5ad")
adata.var_names_make_unique()
adata = adata[adata.obs['slice_id'] == "BZ5", :]
adata.obs['Ground Truth'] = adata.obs['ct']

adata.var_names_make_unique()
# prefilter_specialgenes(adata)
sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=3000)
sc.pp.normalize_per_cell(adata)
sc.pp.log1p(adata)

prefilter_genes(adata, min_cells=3)  # avoiding all genes are zeros
prefilter_specialgenes(adata)
#Normalization
sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=top_genes)
sc.pp.normalize_total(adata, target_sum=1e4)
# sc.pp.log1p(adata)
sc.pp.scale(adata, zero_center=False, max_value=10)

utils.Cal_Spatial_Net(adata, rad_cutoff=400)    #400
adata = train_model.train(adata,4,n_epochs=epoch,use_image = use_image)
adata = utils.mclust_R(adata, used_obsm='stGCL', num_cluster=4)

new_type = refine_label(adata, radius, key='mclust')
adata.obs['stGCL_refined'] = new_type

obs_df = adata.obs.dropna()
ARI_stGCL = adjusted_rand_score(obs_df['mclust'], obs_df['Ground Truth'])
print('Adjusted rand index = %.2f' %ARI_stGCL)

ARI = adjusted_rand_score(obs_df['stGCL_refined'], obs_df['Ground Truth'])
print('Adjusted rand index = %.2f' %ARI)

print('%.5f ---------------------------------------------------------' %ARI)
# adata.obs.to_csv("../results/PA/domains/stgcl_{}_PA_{:.2f}.csv".format(seed,ARI))                                 # change
# adata.obs.to_csv("../results/PA/domains_new/slice3_stgcl_{}_PA_{:.2f}.csv".format(seed,ARI))