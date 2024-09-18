import scanpy as sc
import pandas as pd
from sklearn import metrics
import torch

import matplotlib.pyplot as plt
import seaborn as sns
import SEDR
import os
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

#传参所用
import sys
# seed = int(sys.argv[1])
seed = 666
n_cluster = 4
os.environ['R_HOME'] = '/home/dell/anaconda3/envs/stagate/lib/R'
random_seed = 666
SEDR.fix_seed(random_seed)
device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')

adata = sc.read("../Data/SRARmap_pa/MS_raw_Dataset11_Dataset11_data.h5ad")
adata.var_names_make_unique()
adata = adata[adata.obs['slice_id'] == "BZ5", :]
adata.obs['Ground Truth'] = adata.obs['ct']


sc.pp.filter_genes(adata, min_cells=50)
sc.pp.filter_genes(adata, min_counts=10)
sc.pp.normalize_total(adata, target_sum=1e6)
sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=2000)
adata = adata[:, adata.var['highly_variable'] == True]
sc.pp.scale(adata)
adata_X = adata.X
adata.obsm['X_pca'] = adata_X
graph_dict = SEDR.graph_construction(adata, 8)
sedr_net = SEDR.Sedr(adata.obsm['X_pca'], graph_dict, mode='clustering', device=device)
using_dec = True
if using_dec:
    sedr_net.train_with_dec(N=1)
else:
    sedr_net.train_without_dec(N=1)
sedr_feat, _, _, _ = sedr_net.process()
adata.obsm['SEDR'] = sedr_feat
SEDR.mclust_R(adata, n_cluster, use_rep='SEDR', key_added='SEDR')

sub_adata = adata[~pd.isnull(adata.obs['Ground Truth'])]
ARI = metrics.adjusted_rand_score(sub_adata.obs['Ground Truth'], sub_adata.obs['SEDR'])

plt.rcParams["figure.figsize"] = (3, 3)
sc.pl.embedding(adata, basis="spatial", color="SEDR", s=6, show=False, title='SEDR',save="SEDR_pa.pdf")

print('%.5f ---------------------------------------------------------' %ARI)
# adata.obs.to_csv("../results/PA/domains/sedr_9_PA_%.2f.csv"%ARI)                                 # change
# adata.obs.to_csv("../results/PA/domains/slice2_sedr_{}_PA_{:.2f}.csv".format(seed,ARI))