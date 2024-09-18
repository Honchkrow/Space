import os
import torch
import pandas as pd
import scanpy as sc
from sklearn import metrics
import multiprocessing as mp
from GraphST import GraphST
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics.cluster import adjusted_rand_score
os.environ['R_HOME'] = '/home/dell/anaconda3/envs/stagate/lib/R'


device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
dataset="BARISTASeq"
n_clusters=6
adata = sc.read("../Data/BARISTASeq/BARISTASeq_Sun2021Integrating_Slice_1_data.h5ad")
adata.var_names_make_unique()
adata.obs['Ground Truth'] = adata.obs['layer']
# # Normalization
# sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=3000)
# sc.pp.normalize_total(adata, target_sum=1e4)
# sc.pp.log1p(adata)

model = GraphST.GraphST(adata, device=device ,random_seed=0)                               # change
adata = model.train()

radius = 20
from GraphST.utils import clustering
clustering(adata, n_clusters, radius=radius, refinement=True) #For DLPFC dataset, we use optional refinement step.

df=adata.obs
df.replace("outside_VISp", np.nan, inplace=True)
df = df.dropna()
adata.obs['Ground Truth']=df['Ground Truth']
adata.obs['domain']=df['domain']

Ann_df = adata.obs.dropna()
# Ann_df.to_csv("Data/stagate_st.csv")
ARI = adjusted_rand_score(Ann_df['domain'], Ann_df["Ground Truth"])
print(ARI)
print('=========================')


coor = pd.DataFrame(adata.obsm['spatial'])
coor.index = adata.obs.index
coor.columns = ['imagerow', 'imagecol']
adata.obs["y_pixel"]=coor['imagerow']
adata.obs["x_pixel"]=coor['imagecol']


domains='domain'
title = 'GraphST (ARI=%.2f)' % ARI
ax = sc.pl.scatter(adata, alpha=1, x="y_pixel", y="x_pixel", color=domains, legend_fontsize=18, show=False,
                   size=100000 / adata.shape[0])

ax.set_title(title, fontsize=23)
ax.set_aspect('equal', 'box')
ax.set_xticks([])
ax.set_yticks([])
ax.axes.invert_yaxis()
# plt.savefig("../results/Seq/figures/GraphST_Seq_8.pdf")
plt.close()

# adata.obs.to_csv("../results/Seq/domains/GraphST_8_Seq_%.2f.csv"%ARI)  # change


from sklearn.metrics import normalized_mutual_info_score
nmi=normalized_mutual_info_score(Ann_df['domain'], Ann_df["Ground Truth"])
# print('normalized mutual info score = %.5f' % nmi)