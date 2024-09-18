import warnings
warnings.filterwarnings("ignore")
import MENDER
import scanpy as sc
import pandas as pd
import numpy as np
from sklearn.metrics import *
import time

#传参所用
# import sys
# seed = int(sys.argv[1])
seed = 2
adata_raw = sc.read("../Data/Mouse_hippocampus_MERFISH/hip_adata-0.04.h5ad")

adata_raw.var_names_make_unique()
adata = adata_raw.copy()
adata.obs['Ground Truth'] = adata.obs['ground_truth']
######### determine cell state using standard Leiden [start]  #########
# this step can be optionally skipped if reliable cell type annotation is available
sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=4000)
sc.pp.normalize_total(adata, inplace=True)
# sc.pp.log1p(adata)

sc.pp.pca(adata)
sc.pp.neighbors(adata)
sc.tl.leiden(adata,resolution=2,key_added='ct')
adata.obs['ct'] = adata.obs['ct'].astype('category')
######### determine cell state using standard Leiden [end]  #########
# input parameters of MENDER
scale = 2

# main body of MENDER
msm = MENDER.MENDER_single(
    adata,
    # determine which cell state to use
    # we use the cell state got by Leiden
    ct_obs='ct' ,

    random_seed= seed                                                                               # change
)
# set the MENDER parameters
msm.set_MENDER_para(
    # default of n_scales is 6
    n_scales=scale,

    # for single cell data, nn_mode is set to 'radius'
    # for spot data, nn_mode is set to 'ring', since each spot is surrounded by certain number of spots (6 for visium and 4 for ST)

    nn_mode='ring',

    # default of n_scales is 15 um (see the manuscript for why).
    # MENDER also provide a function 'estimate_radius' for estimating the radius
    # if nn_mode is set to 'ring', nn_para means the number of spots around the central spot, i.e., 6 for Visium and 4 for ST
    nn_para=6,
)
# construct the context representation
msm.run_representation(

    # the number of processings
)
# set the spatial clustering parameter
# positive values for the expected number of domains
# negative values for the clustering resolution
msm.run_clustering_normal(8)

obs_df = msm.adata_MENDER.obs.dropna()
ARI = adjusted_rand_score(obs_df['MENDER'], obs_df['Ground Truth'])
print('Adjusted rand index = %.5f' % ARI)
from sklearn.metrics import normalized_mutual_info_score
nmi=normalized_mutual_info_score(obs_df['MENDER'], obs_df['Ground Truth'])
print('normalized mutual info score = %.5f' % nmi)

sc.pl.embedding(msm.adata_MENDER, basis="spatial", color="MENDER", s=6, show=False, title='leiden',save="mender_hip.pdf")

# msm.adata_MENDER.obs.to_csv("../results/HIP/domains/mender_9_HIP_%.2f.csv"%ARI)  # change
# msm.adata_MENDER.obs.to_csv("../results/HIP/domains/hip024_mender_{}_hip_{:.2f}.csv".format(seed,ARI))
print(ARI,"-------------------------")