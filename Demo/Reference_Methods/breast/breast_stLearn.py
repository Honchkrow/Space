#-*- coding : utf-8 -*-
import os,csv,re
import pandas as pd
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
from sklearn import metrics
from pathlib import Path
import stlearn as st


section_id = "V1_Breast_Cancer_Block_A_Section_1"
k=20
BASE_PATH = Path("../Data/V1_Breast_Cancer_Block_A_Section_1/")
Ann_df =  pd.read_csv("../Data/V1_Breast_Cancer_Block_A_Section_1/metadata.tsv", sep="	", header=0, na_filter=False,
                          index_col=0)
TILE_PATH = Path("../Data/V1_Breast_Cancer_Block_A_Section_1//tmp/tiles")
TILE_PATH.mkdir(parents=True, exist_ok=True)
OUT_PATH = Path("out")
OUT_PATH.mkdir(parents=True, exist_ok=True)
data = st.Read10X(BASE_PATH)
# pre-processing for gene count table
# pre-processing for gene count table
st.pp.filter_genes(data, min_cells=1)
st.pp.normalize_total(data)
st.pp.log1p(data)

# run PCA for gene expression data
st.em.run_pca(data, n_comps=15)

# pre-processing for spot image
st.pp.tiling(data, TILE_PATH)

# this step uses deep learning model to extract high-level features from tile images
# may need few minutes to be completed
st.pp.extract_feature(data)

# stSME
st.spatial.SME.SME_normalize(data, use_data="raw", weights="physical_distance")
data_ = data.copy()
data_.X = data_.obsm['raw_SME_normalized']

st.pp.scale(data_)
st.em.run_pca(data_, n_comps=15)

st.tl.clustering.kmeans(data_, n_clusters=k, use_data="X_pca", key_added="X_pca_kmeans",algorithm= "auto",random_state=5)      # change
st.pl.cluster_plot(data_, use_label="X_pca_kmeans")

# louvain clustering on stSME normalised data
st.pp.neighbors(data_,n_neighbors=k,use_rep='X_pca')
st.tl.clustering.louvain(data_, resolution=1.19)
st.pl.cluster_plot(data_,use_label="louvain")

Ann_df["kmeans"]=data_.obs["X_pca_kmeans"]
Ann_df['louvain']=data_.obs['louvain']
Ann_df.dropna(inplace=True)
ari = metrics.adjusted_rand_score(Ann_df["fine_annot_type"], Ann_df["kmeans"])
ari2 = metrics.adjusted_rand_score(Ann_df["fine_annot_type"], Ann_df['louvain'])
print(ari,ari2)
from sklearn.metrics import normalized_mutual_info_score
nmi=normalized_mutual_info_score(Ann_df["fine_annot_type"], Ann_df['louvain'])
print('normalized mutual info score = %.5f' % nmi)
nmi=normalized_mutual_info_score(Ann_df["fine_annot_type"], Ann_df['kmeans'])
print('normalized mutual info score = %.5f' % nmi)
plt.rcParams["figure.figsize"] = (3, 3)
sc.pl.spatial(data_, color="louvain", title='stLearn (ARI=%.2f)' % ari,
           save="stLearnloubreast")

plt.rcParams["figure.figsize"] = (3, 3)
sc.pl.spatial(data_, color="X_pca_kmeans", title='stLearn (ARI=%.2f)' % ari,
           save="stLearnkmeansbreast")

print('%.5f ---------------------------------------------------------' %ari)
# data_.obs.to_csv("../results/breast/domains/stlearn_6_breast_%.2f.csv"%ari)                                 # change