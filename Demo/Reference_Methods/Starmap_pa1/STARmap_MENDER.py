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
seed = 17

adata_raw = sc.read("../Data/SRARmap_pa/MS_raw_Dataset11_Dataset11_data.h5ad")
adata_raw.var_names_make_unique()
adata = adata_raw.copy()
adata = adata[adata.obs['slice_id'] == "BZ5", :] # 新加的
batch_obs = 'slice_id'
gt_obs = 'gt'
# input parameters of MENDER
scale = 6
# the default radius is 15 um, the unit of spatial coordinate of this dataset is 0.1 um
radius = 150

# estimate number of domains
n_cls = np.unique(adata.obs[gt_obs]).shape[0]


# record running time
time_st = time.time()


adata = adata_raw.copy()


# main body of MENDER
msm = MENDER.MENDER(
    adata,
    batch_obs = batch_obs,
    # determine which cell state to use
    # we use the cell state got by Leiden
    ct_obs='ct' ,
    random_seed= seed
)


# set the MENDER parameters


msm.prepare()
msm.set_MENDER_para(
    # default of n_scales is 6
    n_scales=scale,

    # for single cell data, nn_mode is set to 'radius'
    nn_mode='radius',

    # default of n_scales is 15 um (see the manuscript for the analysis).
    # MENDER also provide a function 'estimate_radius' for estimating the radius
    nn_para=radius,

)
# construct the context representation
msm.run_representation_mp(
    200
    # the number of processings
)

# set the spatial clustering parameter
# positive values for the expected number of domains
# negative values for the clustering resolution
msm.run_clustering_normal(n_cls)

time_ed = time.time()
time_cost = time_ed-time_st
msm.output_cluster_all(obs='MENDER',obs_gt=gt_obs)
print('MENDER prediction')
msm.output_cluster_all(obs=gt_obs,obs_gt=None)
print('ground truth')
obs_df = msm.adata_MENDER.obs.dropna()
ARI = adjusted_rand_score(obs_df['MENDER'], obs_df['gt'])
print('Adjusted rand index = %.5f' % ARI)
from sklearn.metrics import normalized_mutual_info_score
nmi=normalized_mutual_info_score(obs_df['MENDER'], obs_df['gt'])
print('normalized mutual info score = %.5f' % nmi)

sc.pl.embedding(msm.adata_MENDER, basis="spatial", color="MENDER", s=6, show=False, title='mender',save="mender_hip.pdf")
# msm.adata_MENDER.obs.to_csv("../results/PA/domains/mender_9_PA_%.2f.csv"%ARI)  # change
print(f'running time: {time_cost} s')

print(ARI,"===================")
# adata.obs.to_csv("../results/PA/domains/slice2_mender_{}_PA_{:.2f}.csv".format(seed,ARI))