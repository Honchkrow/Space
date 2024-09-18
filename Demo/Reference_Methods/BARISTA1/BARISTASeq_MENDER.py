import warnings
warnings.filterwarnings("ignore")
import MENDER
import scanpy as sc
import pandas as pd
import numpy as np
from sklearn.metrics import *
import time

adata_raw= sc.read("../Data/BARISTASeq/BARISTASeq_Sun2021Integrating_Slice_1_data.h5ad")
adata_raw.var_names_make_unique()
adata = adata_raw.copy()
adata.obs['Ground Truth'] = adata.obs['layer']
######### determine cell state using standard Leiden [start]  #########
# this step can be optionally skipped if reliable cell type annotation is available
sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=4000)
sc.pp.normalize_total(adata, inplace=True)
# sc.pp.log1p(adata)

sc.pp.pca(adata)
sc.pp.neighbors(adata)
sc.tl.leiden(adata,resolution=2,key_added='ct')                                                                         #change(2)
adata.obs['ct'] = adata.obs['ct'].astype('category')
######### determine cell state using standard Leiden [end]  #########
# input parameters of MENDER
scale = 6

# main body of MENDER
msm = MENDER.MENDER_single(
    adata,
    # determine which cell state to use
    # we use the cell state got by Leiden
    ct_obs='ct' ,

    random_seed= 17                                                                                                      # change
)
# set the MENDER parameters
msm.set_MENDER_para(
    # default of n_scales is 6
    n_scales=scale,

    # for single cell data, nn_mode is set to 'radius'
    nn_mode='radius',

    # default of n_scales is 15 um (see the manuscript for why).
    # MENDER also provide a function 'estimate_radius' for estimating the radius
    nn_para=15,

)
# construct the context representation
msm.run_representation(

    # the number of processings
)
# set the spatial clustering parameter
# positive values for the expected number of domains
# negative values for the clustering resolution
msm.run_clustering_normal(6)

obs_df = msm.adata_MENDER.obs.dropna()
obs_df = obs_df[obs_df['Ground Truth'] != 'outside_VISp']
ARI = adjusted_rand_score(obs_df['MENDER'], obs_df['Ground Truth'])
print('Adjusted rand index = %.5f' % ARI)
from sklearn.metrics import normalized_mutual_info_score
nmi=normalized_mutual_info_score(obs_df['MENDER'], obs_df['Ground Truth'])
print('normalized mutual info score = %.5f' % nmi)
sc.pl.scatter(msm.adata_MENDER, alpha=1, x="x", y="y", color="MENDER", legend_fontsize=18, show=True,save="mender_BARISTASeq.pdf",
                   size=100000 / msm.adata_MENDER.shape[0])

print(ARI)
# msm.adata_MENDER.obs.to_csv("../results/Seq/domains/slice1_mender_0_Seq_%.2f.csv"%ARI)  # change