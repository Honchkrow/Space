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

#传参所用
import sys
# seed = int(sys.argv[1])
seed = 0
device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
dataset="PA"
n_clusters=4
adata = sc.read("../Data/SRARmap_pa/MS_raw_Dataset11_Dataset11_data.h5ad")
adata.var_names_make_unique()
adata = adata[adata.obs['slice_id'] == "BZ5", :]
# sc.pl.embedding(adata,basis="spatial",color="gt")
adata.obs['Ground Truth'] = adata.obs['gt']
# # Normalization
sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=3000)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

model = GraphST.GraphST(adata, device=device ,random_seed= seed)                               # change
adata = model.train()

radius = 20
from GraphST.utils import clustering
clustering(adata, n_clusters, radius=radius, refinement=True) #For DLPFC dataset, we use optional refinement step.

sc.pl.embedding(adata,basis="spatial",color="domain")


Ann_df = adata.obs.dropna()
# Ann_df.to_csv("Data/stagate_st.csv")
ARI = adjusted_rand_score(Ann_df['domain'], Ann_df["Ground Truth"])
print(ARI,"=============================")

# adata.obs.to_csv("../results/PA/domains/GraphST_9_PA_%.2f.csv"%ARI)  # change
# adata.obs.to_csv("../results/PA/domains/slice2_graphst_{}_PA_{:.2f}.csv".format(seed,ARI))