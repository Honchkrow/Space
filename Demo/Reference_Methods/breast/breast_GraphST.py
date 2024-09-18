import os
import torch
import pandas as pd
import scanpy as sc
from sklearn import metrics
import multiprocessing as mp
from GraphST import GraphST
import numpy as np
os.environ['R_HOME'] = '/home/dell/anaconda3/envs/stagate/lib/R'

device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
ARIlist = []
dataset="V1_Breast_Cancer_Block_A_Section_1"
n_clusters=20
print(dataset)
input_dir = os.path.join('Data', dataset)
adata = sc.read_visium(path="../Data/V1_Breast_Cancer_Block_A_Section_1",
                       count_file='filtered_feature_bc_matrix.h5')
adata.var_names_make_unique()


import time
import tracemalloc
start_time = time.time()
tracemalloc.start()


model = GraphST.GraphST(adata, device=device ,random_seed= 4)                               # change
adata = model.train()

radius = 50
from GraphST.utils import clustering
clustering(adata, n_clusters, radius=radius, refinement=True) #For DLPFC dataset, we use optional refinement step.


current, peak = tracemalloc.get_traced_memory()
end_time = time.time()
tracemalloc.stop()
elapsed_time = end_time - start_time
print('train_model.train() 函数的运行时间 = %.2f 秒' % elapsed_time)
print(f"当前内存使用: {current / 10**6:.2f} MB")
print(f"内存峰值使用: {peak / 10**6:.2f} MB")

df_meta  = pd.read_csv("../Data/V1_Breast_Cancer_Block_A_Section_1/metadata.tsv", sep="	", header=0,
                     na_filter=False,index_col=0)
adata.obs['fine_annot_type'] = df_meta.loc[adata.obs_names, 'fine_annot_type']


# filter out NA nodes
adata = adata[~pd.isnull(adata.obs['fine_annot_type'])]

# calculate metric ARI
ARI = metrics.adjusted_rand_score(adata.obs['domain'], adata.obs['fine_annot_type'])
adata.uns['ARI'] = ARI
ARIlist.append(ARI)
print('Dataset:', dataset)
print('ARI:', ARI)

# adata.obs.to_csv("../results/breast/domains/graphst_9_breast_%.2f.csv"%ARI)              # change