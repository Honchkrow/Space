import os
import scanpy as sc
import pandas as pd

os.environ["R_HOME"] = "/home/zw/software/miniforge-pypy3/envs/space/lib/R"
import Space
from Space.cons_func import (
    run_GraphST,
    run_Leiden,
    run_MENDER,
    run_SCANPY,
    run_SEDR,
    run_SpaceFlow,
    run_SpaGCN,
    run_STAGATE,
    run_stGCL,
    run_stLearn,
)
from sklearn.metrics import adjusted_rand_score
from sklearn.cluster import SpectralClustering
from Space.cons_func import get_results, get_domains
from Space.utils import calculate_location_adj, plot_results_ari, get_bool_martix

adata = sc.read_visium(
    path="./Data/V1_Breast_Cancer_Block_A_Section_1", count_file="filtered_feature_bc_matrix.h5"
)
Ann_df = pd.read_csv(
    "./Data/V1_Breast_Cancer_Block_A_Section_1/metadata.tsv",
    sep="	",
    header=0,
    na_filter=False,
    index_col=0,
)
adata.var_names_make_unique()
im_re = pd.read_csv(
    "./Data/V1_Breast_Cancer_Block_A_Section_1/image_representation/ViT_pca_representation.csv",
    header=0,
    index_col=0,
    sep=",",
)
adata.obsm["im_re"] = im_re
adata.obs["gt"] = Ann_df["fine_annot_type"]
gt = adata.obs["gt"]
k = 20  # n_clusters
epochs = 120 # 100是0.63 120 是0.65
seed = 666
alpha = 1
learning_rate = 0.0001

# Define a dictionary where the key is a subfunction and the value is a dictionary of parameters
func_dict = {
    # run_GraphST: {"input_adata": adata, "n_cluster": k, "radius": 50, "random_seed": 5},
    # run_Leiden: {"input_adata": adata, "random_state": 0},
    # run_MENDER: {"input_adata": adata, "n_cluster": k, "scale": 2, "radius": 6, "random_seed": 1},
    # run_SCANPY: {"input_adata": adata, "resolution": 0.5, "random_state": 0},
    # run_SEDR: {"input_adata": adata, "n_cluster": k, "random_seed": 4},
    # run_SpaceFlow: {"input_adata": adata, "resolution": 1.5, "n_neighbors": 50, "random_seed": 8},
    # run_STAGATE: {"input_adata": adata, "n_cluster": k, "random_seed": 3},
    # run_stGCL: {
    #     "input_adata": adata,
    #     "n_cluster": k,
    #     "radius": 70,
    #     "epoch": 100,
    #     "use_image": True,
    #     "seed": 0,
    # },
    # run_stLearn: {"input_adata": adata, "n_cluster": k, "random_state": 5},
    # run_SpaGCN: {"input_adata": adata, "max_epochs": 20, "seed": 4},
}

# To save time and ensure the stability of the results, we can use the results we obtained in advance:
# mul_reults : "../Data/V1_Breast_Cancer_Block_A_Section_1/result.csv"
mul_reults = pd.read_csv(
    "./Data/V1_Breast_Cancer_Block_A_Section_1/result.csv", header=0, index_col=0
)
mul_reults = mul_reults.iloc[:, 2:]
# or
"""collected_results = get_results(func_dict, use_multithreading=False,monitor_performance=True)
mul_reults = get_domains(adata,collected_results,gt)
mul_reults = mul_reults.drop('ground_truth', axis=1)"""


plot_results_ari(mul_reults)
mul_reults = mul_reults.drop("SpaceFlow", axis=1)
mul_reults = mul_reults.drop("MENDER", axis=1)

pos_similarity = calculate_location_adj(adata.obsm["spatial"], l=123)

model = Space.Space(
    get_bool_martix(mul_reults),
    pos_similarity,
    epochs=epochs,
    gt=gt.values,
    k=k,
    seed=seed,
    alpha=alpha,
    beta=1,
    learning_rate=learning_rate,
)
con_martix = model.train()
sc = SpectralClustering(n_clusters=k, affinity="precomputed", random_state=666)
labels = sc.fit_predict(con_martix)
adata.obs["consensus"] = labels
ari = adjusted_rand_score(labels, gt.values)
print(ari)
