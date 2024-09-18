import warnings
warnings.filterwarnings("ignore")
import MENDER
import scanpy as sc
import pandas as pd
import numpy as np
from sklearn.metrics import *
import time

adata_raw = sc.read_visium(path="../Data/V1_Breast_Cancer_Block_A_Section_1",
                       count_file='filtered_feature_bc_matrix.h5')
adata_raw.var_names_make_unique()
adata = adata_raw.copy()



######### determine cell state using standard Leiden [start]  #########
# this step can be optionally skipped if reliable cell type annotation is available
sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=4000)
sc.pp.normalize_total(adata, inplace=True)
sc.pp.log1p(adata)

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

    random_seed= 9                                                                               # change
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
msm.run_clustering_normal(-0.09)

Ann_df = pd.read_csv("../Data/V1_Breast_Cancer_Block_A_Section_1/metadata.tsv", sep="	", header=0,
                     na_filter=False,index_col=0)
msm.adata_MENDER.obs['Ground Truth'] = Ann_df.loc[adata.obs_names, 'fine_annot_type']
obs_df = msm.adata_MENDER.obs.dropna()
ARI = adjusted_rand_score(obs_df['MENDER'], obs_df['Ground Truth'])
print('Adjusted rand index = %.5f' % ARI)
from sklearn.metrics import normalized_mutual_info_score
nmi=normalized_mutual_info_score(obs_df['MENDER'], obs_df['Ground Truth'])
print('normalized mutual info score = %.5f' % nmi)
sc.pl.spatial(msm.adata_MENDER,color='MENDER',size=1,save="mender_breast.pdf" ,title='MENDER (ARI=%.2f)' % ARI)

# adata.obs.to_csv("../results/breast/domains/mender_9_breast_%.2f.csv"%ARI)  # change