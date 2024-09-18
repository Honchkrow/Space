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
# import sys
seed = 4

device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
dataset="HIP"
n_clusters=8
adata = sc.read("../Data/Mouse_hippocampus_MERFISH/hip_adata-0.04.h5ad")
adata.var_names_make_unique()
# # Normalization
# sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=3000)
# sc.pp.normalize_total(adata, target_sum=1e4)
# sc.pp.log1p(adata)

model = GraphST.GraphST(adata, device=device ,random_seed= seed)                               # change
adata = model.train()

radius = 20
from GraphST.utils import clustering
clustering(adata, n_clusters, radius=radius, refinement=True) #For DLPFC dataset, we use optional refinement step.

adata.obs['Ground Truth'] = adata.obs['ground_truth']

Ann_df = adata.obs.dropna()
# Ann_df.to_csv("Data/stagate_st.csv")
ARI = adjusted_rand_score(Ann_df['domain'], Ann_df["Ground Truth"])
print(ARI)

# adata.obs.to_csv("../results/HIP/domains/GraphST_9_hip_%.2f.csv"%ARI)  # change
# adata.obs.to_csv("../results/HIP/domains/hip024_graphst_{}_hip_{:.2f}.csv".format(seed,ARI))