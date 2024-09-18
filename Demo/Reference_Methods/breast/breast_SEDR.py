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
os.environ['R_HOME'] = '/home/dell/anaconda3/envs/stagate/lib/R'
random_seed = 9                                                                                      # change
SEDR.fix_seed(random_seed)
device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
ARIlist = []
dataset="V1_Breast_Cancer_Block_A_Section_1"
n_clusters=20
print(dataset)
input_dir = os.path.join('Data', dataset)
adata = sc.read_visium(path="../Data/V1_Breast_Cancer_Block_A_Section_1",
                       count_file='filtered_feature_bc_matrix.h5')
adata.var_names_make_unique()
df_meta  = pd.read_csv("../Data/V1_Breast_Cancer_Block_A_Section_1/metadata.tsv", sep="	", header=0,
                     na_filter=False,index_col=0)
adata.obs['Ground_Truth'] = df_meta.loc[adata.obs_names, 'fine_annot_type']
sc.pp.filter_genes(adata, min_cells=50)
sc.pp.filter_genes(adata, min_counts=10)
sc.pp.normalize_total(adata, target_sum=1e6)
sc.pp.highly_variable_genes(adata, flavor="seurat_v3",  n_top_genes=2000)
adata = adata[:, adata.var['highly_variable'] == True]
sc.pp.scale(adata)

from sklearn.decomposition import PCA  # sklearn PCA is used because PCA in scanpy is not stable.
adata_X = PCA(n_components=200, random_state=42).fit_transform(adata.X)
adata.obsm['X_pca'] = adata_X
graph_dict = SEDR.graph_construction(adata, 12)
print(graph_dict)
sedr_net = SEDR.Sedr(adata.obsm['X_pca'], graph_dict, mode='clustering', device=device)
using_dec = True
if using_dec:
    sedr_net.train_with_dec(N=1)
else:
    sedr_net.train_without_dec(N=1)
sedr_feat, _, _, _ = sedr_net.process()
adata.obsm['SEDR'] = sedr_feat
SEDR.mclust_R(adata, n_clusters, use_rep='SEDR', key_added='SEDR')

sub_adata = adata[~pd.isnull(adata.obs['Ground_Truth'])]
ARI = metrics.adjusted_rand_score(sub_adata.obs['Ground_Truth'], sub_adata.obs['SEDR'])

fig, axes = plt.subplots(1,2,figsize=(4*2, 4))
sc.pl.spatial(adata, color='Ground_Truth', ax=axes[0], show=False)
sc.pl.spatial(adata, color='SEDR', ax=axes[1], show=False)
axes[0].set_title('Manual Annotation')
axes[1].set_title('Clustering: (ARI=%.4f)' % ARI)
plt.tight_layout()
plt.show()


sc.pp.neighbors(adata, use_rep='SEDR', metric='cosine')
sc.tl.umap(adata)

fig, axes = plt.subplots(1,2,figsize=(4*2, 3))
sc.pl.umap(adata, color='Ground_Truth', ax=axes[0], show=False)
sc.pl.umap(adata, color='SEDR', ax=axes[1], show=False)
axes[0].set_title('Manual Annotation')
axes[1].set_title('Clustering')

for ax in axes:
    ax.set_aspect(1)

plt.tight_layout()
plt.show()

print('%.5f ---------------------------------------------------------' %ARI)
# adata.obs.to_csv("../results/breast/domains/sedr_9_breast_%.2f.csv"%ARI)                                 # change