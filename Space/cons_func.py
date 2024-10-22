import scanpy as sc
import torch
import warnings
warnings.filterwarnings('ignore')
import pandas as pd
import matplotlib.pyplot as plt
import tensorflow as tf
tf.compat.v1.disable_eager_execution()
import ot
import random
import numpy as np
import concurrent.futures


def run_GraphST(input_adata, n_cluster ,radius = 20,random_seed=0):
    print("run_GraphST")
    from GraphST import GraphST
    device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
    adata = input_adata.copy()
    model = GraphST.GraphST(adata, device=device, random_seed=random_seed)
    adata = model.train()
    from GraphST.utils import clustering
    clustering(adata, n_cluster, radius=radius, refinement=True)

    return adata.obs['domain'], adata.obsm['emb_pca']

def run_Leiden(input_adata, resolution=0.5,n_top_genes=3000,random_state=0,n_comps=30):
    print("run_Leiden")
    adata = input_adata.copy()
    sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=n_top_genes)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    sc.pp.pca(adata, n_comps=n_comps)
    sc.pp.neighbors(adata, use_rep='X_pca')
    sc.tl.leiden(adata, resolution=resolution)
    # sc.tl.umap(adata)

    return adata.obs['leiden'], adata.obsm['X_pca']

def run_MENDER(input_adata, n_cluster, resolution=2,n_top_genes=4000,scale = 6,random_seed=17,radius = 150):
    print("run_MENDER")
    import Space.MENDER as MENDER
    adata = input_adata.copy()
    sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=n_top_genes)
    sc.pp.normalize_total(adata, inplace=True)
    sc.pp.pca(adata)
    sc.pp.neighbors(adata)
    sc.tl.leiden(adata, resolution=resolution, key_added='ct')
    adata.obs['ct'] = adata.obs['ct'].astype('category')
    msm = MENDER.MENDER_single(adata, ct_obs='ct', random_seed=random_seed)
    msm.set_MENDER_para(n_scales=scale, nn_mode='radius', nn_para=radius )
    msm.run_representation()
    msm.run_clustering_normal(n_cluster)
    mender_emb = msm.adata_MENDER.obsm["X_pca"]
    return msm.adata_MENDER.obs['MENDER'], mender_emb

def run_SCANPY(input_adata, n_comps=30, resolution=0.5,random_state=1):
    print("run_SCANPY")
    adata = input_adata.copy()
    sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=3000)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    sc.pp.pca(adata, n_comps=n_comps)
    sc.pp.neighbors(adata, use_rep='X_pca')
    sc.tl.louvain(adata, resolution=resolution)

    return adata.obs['louvain'],adata.obsm['X_pca']

def run_SEDR(input_adata, n_cluster, random_seed = 0,n_top_genes=2000,n=12,n_components=200):
    print("run_SEDR")
    import Space.SEDR as SEDR
    random_seed = random_seed
    SEDR.fix_seed(random_seed)
    device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
    adata = input_adata.copy()
    sc.pp.filter_genes(adata, min_cells=50)
    sc.pp.filter_genes(adata, min_counts=10)
    sc.pp.normalize_total(adata, target_sum=1e6)
    sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=n_top_genes)
    adata = adata[:, adata.var['highly_variable'] == True]
    sc.pp.scale(adata)
    from sklearn.decomposition import PCA  # sklearn PCA is used because PCA in scanpy is not stable.
    adata_X = PCA(n_components=n_components, random_state=42).fit_transform(adata.X)
    adata.obsm['X_pca'] = adata_X
    graph_dict = SEDR.graph_construction(adata, n)
    sedr_net = SEDR.Sedr(adata.obsm['X_pca'], graph_dict, mode='clustering', device=device)
    using_dec = True
    if using_dec:
        sedr_net.train_with_dec(N=1)
    else:
        sedr_net.train_without_dec(N=1)
    sedr_feat, _, _, _ = sedr_net.process()
    adata.obsm['SEDR'] = sedr_feat
    SEDR.mclust_R(adata, n_cluster, use_rep='SEDR', key_added='SEDR')

    return adata.obs['SEDR'], adata.obsm['SEDR']

def run_SpaceFlow(input_adata, n_top_genes=3000, resolution=0.5,epochs=1000,random_seed=3, n_neighbors=10):
    print("run_SpaceFlow")
    import Space.SpaceFlow as SpaceFlow
    adata = input_adata.copy()
    sc.pp.filter_genes(adata, min_cells=3)
    sf = SpaceFlow.SpaceFlow(adata=adata)
    sf.preprocessing_data(n_top_genes=n_top_genes, n_neighbors=n_neighbors)
    emb = sf.train(spatial_regularization_strength=0.1, z_dim=50, lr=1e-3, epochs=epochs, max_patience=50,
             embedding_save_filepath="./embedding.csv",
             min_stop=100, random_seed=random_seed, gpu=0, regularization_acceleration=True, edge_subset_sz=1000000)

    domain = sf.segmentation(n_neighbors=50, resolution=resolution)

    adata.obs['SpaceFlow'] = domain.values
    adata.obs['SpaceFlow'] = adata.obs['SpaceFlow'].astype(str)

    return adata.obs['SpaceFlow'], emb

def run_STAGATE(input_adata, n_cluster ,rad_cutoff=400, n_top_genes=3000,random_seed=0,  alpha=0):
    print("run_STAGATE")
    import Space.STAGATE as STAGATE

    adata = input_adata.copy()
    sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=n_top_genes)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    STAGATE.Cal_Spatial_Net(adata, rad_cutoff=rad_cutoff)
    # STAGATE.Stats_Spatial_Net(adata)
    adata = STAGATE.train_STAGATE(adata, alpha=alpha, random_seed=random_seed)#
    adata = STAGATE.mclust_R(adata, used_obsm='STAGATE', num_cluster=n_cluster)

    return adata.obs['mclust'], adata.obsm['STAGATE']

def run_stGCL(input_adata, n_cluster ,seed = 9,radius = 50, top_genes = 3000, epoch = 1000,rad_cutoff=400, use_image= False):
    print("run_stGCL")
    import Space.stGCL as stGCL
    from Space.stGCL.process import prefilter_genes, prefilter_specialgenes, set_seed, refine_nearest_labels#
    from Space.stGCL import utils, train_model#
    use_image = use_image
    set_seed(seed)
    top_genes = top_genes
    epoch = epoch
    def refine_label(adata, radius=50, key='cluster'):
        n_neigh = radius
        new_type = []
        old_type = adata.obs[key].values
        position = adata.obsm['spatial']
        distance = ot.dist(position, position, metric='euclidean')
        n_cell = distance.shape[0]
        for i in range(n_cell):
            vec = distance[i, :]
            index = vec.argsort()
            neigh_type = []
            for j in range(1, n_neigh + 1):
                neigh_type.append(old_type[index[j]])
            max_type = max(neigh_type, key=neigh_type.count)
            new_type.append(max_type)
        new_type = [str(i) for i in list(new_type)]
        return new_type
    adata = input_adata.copy()
    prefilter_genes(adata, min_cells=3)
    prefilter_specialgenes(adata)
    sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=top_genes)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.scale(adata, zero_center=False, max_value=10)
    utils.Cal_Spatial_Net(adata, rad_cutoff=rad_cutoff)
    adata = train_model.train(adata, knn=n_cluster, n_epochs=epoch, use_image=use_image)
    adata = utils.mclust_R(adata, used_obsm='stGCL', num_cluster=n_cluster)
    new_type = refine_label(adata, radius, key='mclust')
    adata.obs['stGCL_refined'] = new_type

    return adata.obs['stGCL_refined'], adata.obsm['stGCL']

def run_stLearn(input_adata, n_cluster,TILE_PATH="./stlearn_tiles/",
                n_comps=15,min_cells=1,
                random_state=9,used_image= None):
    print("run_stLearn")
    import stlearn as st
    data = input_adata.copy()
    st.pp.filter_genes(data, min_cells=min_cells)
    st.pp.normalize_total(data)
    st.pp.log1p(data)
    if used_image is not None:
        st.em.run_pca(data, n_comps=n_comps)
        library_id = list(data.uns["spatial"].keys())[0]
        data.uns["spatial"][library_id]["use_quality"]=used_image
        scale = data.uns["spatial"][library_id]["scalefactors"][
            "tissue_" + used_image + "_scalef"]
        image_coor = data.obsm["spatial"] * scale
        data.obs["imagerow"] = image_coor[:, 0]
        data.obs["imagecol"] = image_coor[:, 1]
        st.pp.tiling(data, TILE_PATH)
        st.pp.extract_feature(data)
        st.spatial.SME.SME_normalize(data, use_data="raw", weights="physical_distance")
        data.X = data.obsm['raw_SME_normalized']

    st.pp.scale(data)
    st.em.run_pca(data, n_comps=n_comps)
    st.tl.clustering.kmeans(data, n_clusters=n_cluster, use_data="X_pca", key_added="X_pca_kmeans", algorithm="auto",
                            random_state=random_state)
    return data.obs['X_pca_kmeans'], data.obsm['X_pca']

def run_SpaGCN (input_adata, p=0.5,res=0.5,seed=666,max_epochs=200,s = 1,b = 49):
    print("run_SpaGCN")
    import  Space.SpaGCN as spg#

    adata = input_adata.copy()

    adata.var_names_make_unique()
    spg.prefilter_genes(adata, min_cells=3)
    # spg.prefilter_specialgenes(adata)
    sc.pp.normalize_per_cell(adata)
    sc.pp.log1p(adata)

    spatial = adata.obsm['spatial']
    adata.obs["x_pixel"] = spatial[:, 0]
    adata.obs["y_pixel"] = spatial[:, 1]
    adata.obs["x_pixel"] = spatial[:, 0]
    adata.obs["y_pixel"] = spatial[:, 1]
    x_pixel = adata.obs["x_pixel"].tolist()
    y_pixel = adata.obs["y_pixel"].tolist()

    adj = spg.calculate_adj_matrix(x=x_pixel, y=y_pixel, x_pixel=x_pixel, y_pixel=y_pixel, beta=b, alpha=s,
                                   histology=False)
    l = spg.search_l(p, adj, start=0.01, end=1000, tol=0.01, max_run=100)
    r_seed = t_seed = n_seed = seed
    clf = spg.SpaGCN()
    clf.set_l(l)
    random.seed(r_seed)
    torch.manual_seed(t_seed)
    np.random.seed(n_seed)
    clf.train(adata, adj, init_spa=True, init="louvain", res=res, tol=5e-3, lr=0.05, max_epochs=max_epochs)
    y_pred, prob = clf.predict()
    adata.obs["pred"] = y_pred
    adata.obs["pred"] = adata.obs["pred"].astype('category')

    return(adata.obs['pred'],None)


def call_func_with_params(func, params, monitor_performance):
    import time
    import memory_profiler
    import GPUtil

    start_time = time.time() if monitor_performance else None
    mem_usage_start = memory_profiler.memory_usage()[0] if monitor_performance else None
    gpu_usage_start = GPUtil.getGPUs()[0].memoryUsed if (monitor_performance and GPUtil.getGPUs()) else None
    result = func(**params)
    mem_usage_end = memory_profiler.memory_usage()[0] if monitor_performance else None
    gpu_usage_end = GPUtil.getGPUs()[0].memoryUsed if (monitor_performance and GPUtil.getGPUs()) else None
    end_time = time.time() if monitor_performance else None
    performance_data = {} #记录运行时间和内存使用
    if monitor_performance:
        execution_time = end_time - start_time if start_time and end_time else None
        memory_used = mem_usage_end - mem_usage_start if mem_usage_start and mem_usage_end else None
        gpu_memory_used = gpu_usage_end - gpu_usage_start if gpu_usage_start and gpu_usage_end else None
        performance_data = {
            'execution_time': execution_time,
            'memory_used': memory_used,
            'gpu_memory_used': gpu_memory_used
        }

    return func.__name__, result, performance_data


def get_results(func_dict, use_multithreading=False, monitor_performance=False):
    results = {}
    if use_multithreading:
        with concurrent.futures.ThreadPoolExecutor() as executor:
            future_to_func = {
                executor.submit(call_func_with_params, func, params, monitor_performance): func
                for func, params in func_dict.items()
            }
            for future in concurrent.futures.as_completed(future_to_func):
                func_name, result, performance_data = future.result()
                func_name = func_name.replace('run_', '')
                results[func_name] = {
                    'performance_data': performance_data,
                    'cluster_result': result[0],
                    'emb': result[1]
                }
    else:
        for func, params in func_dict.items():
            func_name, result, performance_data = call_func_with_params(func, params, monitor_performance)
            func_name = func_name.replace('run_', '')
            results[func_name] = {
                'performance_data': performance_data,
                'cluster_result': result[0],
                'emb':result[1]
            }
    return results

def get_domains(adata,results,gt=None):
    combined_df = pd.DataFrame(index=adata.obs_names)

    for func_name, data in results.items():
        result = data['cluster_result']
        df = result.to_frame()
        df.columns.values[0] = func_name
        combined_df = pd.concat([combined_df, df], axis=1)
    if gt is not None:
        combined_df["ground_truth"] = gt
    return combined_df