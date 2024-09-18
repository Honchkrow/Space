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

# # the location of R (used for the mclust clustering)
os.environ['R_HOME'] = '/home/dell/anaconda3/envs/stagate/lib/R'
# os.environ['R_USER'] = 'D:\ProgramData\Anaconda3\Lib\site-packages\rpy2'
use_image=True
radius=50
seed=0                                                                                              # change
set_seed(seed)
radius = 70
top_genes = 3000
epoch = 100

section_id = 'Breast_Cancer'

adata = sc.read_visium(path="../Data/V1_Breast_Cancer_Block_A_Section_1",
                       count_file='filtered_feature_bc_matrix.h5')
adata.var_names_make_unique()
prefilter_genes(adata, min_cells=3)  # avoiding all genes are zeros
prefilter_specialgenes(adata)
#Normalization
sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=top_genes)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.scale(adata, zero_center=False, max_value=10)
im_re = pd.read_csv("../Data/V1_Breast_Cancer_Block_A_Section_1/"
                    "image_representation/ViT_pca_representation.csv", header=0, index_col=0, sep=',')
adata.obsm["im_re"] = im_re
utils.Cal_Spatial_Net(adata, rad_cutoff=400)
adata = train_model.train(adata,20,n_epochs=epoch)
adata = utils.mclust_R(adata, used_obsm='stGCL', num_cluster=20)

new_type = refine_nearest_labels(adata, radius, key='mclust')
adata.obs['stGCL_refined'] = new_type

Ann_df =  pd.read_csv("../Data/V1_Breast_Cancer_Block_A_Section_1/metadata.tsv", sep="	", header=0, na_filter=False,
                          index_col=0)
adata.obs['Ground Truth'] = Ann_df.loc[adata.obs_names, 'fine_annot_type']
obs_df = adata.obs.dropna()
ARI_stGCL = adjusted_rand_score(obs_df['mclust'], obs_df['Ground Truth'])
print('Adjusted rand index = %.2f' %ARI_stGCL)


ARI = adjusted_rand_score(obs_df['stGCL_refined'], obs_df['Ground Truth'])
print('Adjusted rand index = %.2f' %ARI)

print('%.5f ---------------------------------------------------------' %ARI)
# adata.obs.to_csv("../results/breast/domains/stgcl_{}_breast_{:.2f}.csv".format(seed,ARI))                                 # change