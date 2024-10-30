import os
import scanpy as sc
import pandas as pd

os.environ["R_HOME"] = "/home/zw/software/miniforge3/envs/space/lib/R"
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
from Space.utils import (
    calculate_location_adj,
    plot_results_ari,
    get_bool_martix,
    plot_ari_with_removal,
    refine_label,
)

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
k = 20
epochs = 150
seed = 666
alpha = 1
learning_rate = 0.0001

# Define a dictionary where the key is a subfunction and the value is a dictionary of parameters
# users should change the parameters according to the data and requirements of each method
func_dict = {
    run_GraphST: {"input_adata": adata, "n_cluster": k, "radius": 50, "random_seed": 200},
    run_Leiden: {"input_adata": adata, "random_state": 0},
    run_MENDER: {"input_adata": adata, "n_cluster": k, "scale": 2, "radius": 6, "random_seed": 1},
    run_SCANPY: {"input_adata": adata, "resolution": 0.5, "random_state": 0},
    run_SEDR: {"input_adata": adata, "n_cluster": k, "random_seed": 0},
    run_SpaceFlow: {"input_adata": adata, "resolution": 1.5, "n_neighbors": 50, "random_seed": 8},
    run_STAGATE: {"input_adata": adata, "n_cluster": k, "random_seed": 6},
    run_stGCL: {
        "input_adata": adata,
        "n_cluster": k,
        "radius": 70,
        "epoch": 100,
        "use_image": True,
        "seed": 0,
    },
    run_stLearn: {
        "input_adata": adata,
        "n_cluster": k,
    },
    run_SpaGCN: {"input_adata": adata, "max_epochs": 20, "seed": 4},
}


collected_results = get_results(func_dict, use_multithreading=False, monitor_performance=True)
mul_reults = get_domains(adata, collected_results)

mul_reults = plot_ari_with_removal(mul_reults, 8)


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
sClustering = SpectralClustering(n_clusters=k, affinity="precomputed", random_state=666)
labels = sClustering.fit_predict(con_martix)
adata.obs["Space"] = labels
adata.obs["Space"] = adata.obs["Space"].astype("str")
ari = adjusted_rand_score(labels, gt.values)
print(ari)

new_type = refine_label(adata, 70, key='Space')
adata.obs['Space_refined'] = new_type

ari = adjusted_rand_score(new_type, gt.values)
print("refined",ari)

