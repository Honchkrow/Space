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

# change slice
slide_id = "19"

adata = sc.read(f"./Data/Mouse_hippocampus_MERFISH/hip_adata-0.{slide_id}.h5ad")
adata.var_names_make_unique()

gt = adata.obs["ground_truth"]
k = 8  # n_clusters
epochs = 300
seed = 666
alpha = 4
learning_rate = 0.0001

# Define a dictionary where the key is a subfunction and the value is a dictionary of parameters
func_dict = {
    run_GraphST: {"input_adata": adata, "n_cluster": k, "random_seed": 4},
    run_Leiden: {"input_adata": adata, "random_state": 0},
    run_MENDER: {"input_adata": adata, "n_cluster": k, "scale": 2, "radius": 6, "random_seed": 2},
    run_SCANPY: {"input_adata": adata, "resolution": 2, "random_state": 3},
    run_SEDR: {"input_adata": adata, "n_cluster": k, "random_seed": 0},
    run_SpaceFlow: {"input_adata": adata, "n_neighbors": 50, "random_seed": 7},
    run_STAGATE: {"input_adata": adata, "n_cluster": k, "rad_cutoff": 60, "random_seed": 1},
    run_stGCL: {
        "input_adata": adata,
        "n_cluster": k,
        "radius": 70,
        "epoch": 800,
        "rad_cutoff": 60,
        "seed": 0,
    },  ##?
    run_stLearn: {"input_adata": adata, "n_cluster": k, "random_state": 0},
    run_SpaGCN: {"input_adata": adata, "seed": 2},
}

# To save time and ensure the stability of the results, we can use the results we obtained in advance:
# mul_reults : "../Data/Mouse_hippocampus_MERFISH/result04.csv"
mul_reults = pd.read_csv(f"./Data/Mouse_hippocampus_MERFISH/result{slide_id}.csv", header=0, index_col=0)
mul_reults = mul_reults.iloc[:, 2:]
# or
"""collected_results = get_results(func_dict, use_multithreading=True, monitor_performance=True)
mul_reults = get_domains(adata,collected_results,gt)
mul_reults = mul_reults.drop('ground_truth', axis=1)"""


plot_results_ari(mul_reults)
mul_reults = mul_reults.drop("SCANPY", axis=1)
mul_reults = mul_reults.drop("SEDR", axis=1)
mul_reults = mul_reults.drop("Graphst", axis=1)
mul_reults = mul_reults.drop("Leiden", axis=1)
mul_reults = mul_reults.drop("MENDER", axis=1)
mul_reults = mul_reults.drop("SpaGCN", axis=1)
mul_reults = mul_reults.drop("stLearn", axis=1)
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
sc = SpectralClustering(n_clusters=k, affinity="precomputed", random_state=epochs)
labels = sc.fit_predict(con_martix)
adata.obs["consensus"] = labels
ari = adjusted_rand_score(labels, adata.obs["ground_truth"].values)
print(ari)
