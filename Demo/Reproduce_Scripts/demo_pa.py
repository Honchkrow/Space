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
from Space.utils import calculate_location_adj, plot_results_ari, get_bool_martix,plot_ari_with_removal


slide_id = "14"
# slide_id = "5"
adata = sc.read("./Data/SRARmap_pa/MS_raw_Dataset11_Dataset11_data.h5ad")
adata = adata[adata.obs["slice_id"] == f"BZ{slide_id}", :]
adata.var_names_make_unique()
gt = adata.obs["gt"]
k = 4  # n_clusters
epochs = 300
seed = 666
alpha = 1
learning_rate = 0.001

# Define a dictionary where the key is a subfunction and the value is a dictionary of parameters
func_dict = {
    run_GraphST: {"input_adata": adata, "n_cluster": k},
    run_Leiden: {
        "input_adata": adata,
    },
    run_MENDER: {"input_adata": adata, "n_cluster": k},
    run_SCANPY: {
        "input_adata": adata,
    },
    run_SEDR: {"input_adata": adata, "n_cluster": k},
    run_SpaceFlow: {
        "input_adata": adata,
    },
    run_STAGATE: {"input_adata": adata, "n_cluster": k},
    run_stGCL: {"input_adata": adata, "n_cluster": k},
    run_stLearn: {"input_adata": adata, "n_cluster": k},
    run_SpaGCN: {
        "input_adata": adata,
    },
}

# To save time and ensure the stability of the results, we can use the results we obtained in advance:
# mul_reults : "../Data/SRARmap_pa/result1.csv"
mul_reults = pd.read_csv(f"Data/SRARmap_pa/result{slide_id}.csv", header=0, index_col=0)
mul_reults = mul_reults.iloc[:, 2:]
# or
"""collected_results = get_results(func_dict, use_multithreading=True, monitor_performance=True)
mul_reults = get_domains(adata,collected_results,gt)
mul_reults = mul_reults.drop('ground_truth', axis=1)"""
# plot_results_ari(mul_reults)
# mul_reults = mul_reults.drop('Leiden', axis=1)
# mul_reults = mul_reults.drop('SpaGCN', axis=1)
# mul_reults = mul_reults.drop('stLearn', axis=1)
# mul_reults = mul_reults.drop('SCANPY', axis=1)


mul_reults = plot_ari_with_removal(mul_reults,4)


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
ari = adjusted_rand_score(labels, adata.obs["gt"].values)
print(ari)
