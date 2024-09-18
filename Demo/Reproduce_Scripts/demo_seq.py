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
slide_id = "2"
adata = sc.read(f"./Data/BARISTASeq/BARISTASeq_Sun2021Integrating_Slice_{slide_id}_data.h5ad")
adata.var_names_make_unique()

gt = adata.obs["layer"]
k = 6
epochs = 300
seed = 666
alpha = 1
learning_rate = 0.001

# Define a dictionary where the key is a subfunction and the value is a dictionary of parameters
func_dict = {
    # run_GraphST: {"input_adata": adata, "n_cluster": k, "random_seed": 0},
    # run_Leiden: {"input_adata": adata, "random_state": 0},
    # run_MENDER: {"input_adata": adata, "n_cluster": k, "radius": 15, "random_seed": 17},
    # run_SCANPY: {"input_adata": adata, "resolution": 2, "random_state": 0},
    run_SEDR:{'input_adata':adata, 'n_cluster': k,'random_seed':0},
    # run_SpaceFlow: {"input_adata": adata, "n_neighbors": 50, "random_seed": 3},
    # run_STAGATE: {"input_adata": adata, "n_cluster": k, "rad_cutoff": 50, "random_seed": 4},
    # run_stGCL: {"input_adata": adata, "n_cluster": k, "radius": 70, "rad_cutoff": 50, "seed": 9},
    # run_stLearn: {"input_adata": adata, "n_cluster": k, "random_state": 0},
    # run_SpaGCN:{'input_adata':adata, 'seed':2},
}

# To save time and ensure the stability of the results, we can use the results we obtained in advance:
# mul_reults : "../Data/BARISTASeq/result1.csv"

# mul_reults = pd.read_csv(f"./Data/BARISTASeq/result{slide_id}.csv", header=0, index_col=0)
# mul_reults = mul_reults.loc[mul_reults["Ground Truth"] != "outside_VISp"]
# mul_reults = mul_reults.iloc[:, 2:].dropna()


collected_results = get_results(func_dict, use_multithreading=False, monitor_performance=True)
mul_reults = get_domains(adata,collected_results,gt)
mul_reults= mul_reults.loc[mul_reults["ground_truth"] != 'outside_VISp']
mul_reults = mul_reults.dropna()
mul_reults = mul_reults.drop('ground_truth', axis=1)



plot_results_ari(mul_reults)
outside_VISp_indices = adata.obs["layer"] == "outside_VISp"
adata = adata[~outside_VISp_indices]
gt = adata.obs["layer"]

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
ari = adjusted_rand_score(labels, adata.obs["layer"].values)
print(ari)
