# 定义子函数，每个子函数返回两个参数，并添加延迟
def run_GraphST(input_adata, n_cluster ,radius = 20,random_seed=0):
    import os
    import torch
    from consensus.GraphST import GraphST#
    from sklearn.metrics.cluster import adjusted_rand_score

    os.environ['R_HOME'] = '/home/zw/software/miniconda3/envs/stpython/lib/R'
    device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
    adata = input_adata.copy()
    model = GraphST.GraphST(adata, device=device, random_seed=random_seed)  # change
    adata = model.train()
    from GraphST.utils import clustering
    clustering(adata, n_cluster, radius=radius, refinement=True)

    return adata.obs['domain'], adata.obsm['emb_pca']

def run_Leiden(input_adata, resolution=0.5,n_top_genes=3000,random_state=0):
    import warnings
    warnings.filterwarnings("ignore")
    import scanpy as sc

    adata = input_adata.copy()
    sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=n_top_genes)
    sc.pp.pca(adata, n_comps=50)
    sc.pp.neighbors(adata, use_rep='X_pca', n_neighbors=50)
    sc.tl.leiden(adata, resolution=resolution, random_state=random_state)
    sc.tl.umap(adata)

    return adata.obs['leiden'], adata.obsm['X_pca']

def run_MENDER(input_adata, n_cluster, resolution=2,n_top_genes=4000,scale = 6,random_seed=17):
    import warnings
    warnings.filterwarnings("ignore")
    import MENDER as MENDER
    import scanpy as sc

    adata = input_adata.copy()
    sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=n_top_genes)
    sc.pp.normalize_total(adata, inplace=True)
    sc.pp.pca(adata)
    sc.pp.neighbors(adata)
    sc.tl.leiden(adata, resolution=resolution, key_added='ct')
    adata.obs['ct'] = adata.obs['ct'].astype('category')
    msm = MENDER.MENDER_single(adata, ct_obs='ct', random_seed=random_seed)
    msm.set_MENDER_para(n_scales=scale, nn_mode='radius', nn_para=15, )
    msm.run_representation()
    msm.run_clustering_normal(n_cluster)
    mender_emb = msm.adata_MENDER.obsm["X_pca"]

    return msm.adata_MENDER.obs['MENDER'], mender_emb

def run_SCANPY(input_adata, n_comps=30, resolution=1.5,random_state=1):
    import warnings
    warnings.filterwarnings("ignore")
    import scanpy as sc

    adata = input_adata.copy()
    sc.pp.pca(adata, n_comps=n_comps)
    sc.pp.neighbors(adata, use_rep='X_pca')
    sc.tl.louvain(adata, resolution=resolution, random_state=random_state)
    sc.tl.umap(adata)

    return adata.obs['louvain'],adata.obsm['X_pca']

def run_SEDR(input_adata, n_cluster, random_seed = 666,n_top_genes=2000):
    import scanpy as sc
    import SEDR
    import torch
    import os
    import warnings
    warnings.filterwarnings('ignore')

    os.environ['R_HOME'] = '/home/zw/software/miniconda3/envs/stpython/lib/R'
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
    adata_X = adata.X
    adata.obsm['X_pca'] = adata_X
    graph_dict = SEDR.graph_construction(adata, 8)
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

def run_SpaceFlow(input_adata, n_top_genes=3000, resolution=0.5,epochs=1000,random_seed=3):
    import scanpy as sc
    from consensus.SpaceFlow import SpaceFlow#
    import pandas as pd
    import matplotlib.pyplot as plt

    adata = input_adata.copy()
    sc.pp.filter_genes(adata, min_cells=3)
    sf = SpaceFlow.SpaceFlow(adata=adata)
    sf.preprocessing_data(n_top_genes=n_top_genes)
    sf.train(spatial_regularization_strength=0.1, z_dim=50, lr=1e-3, epochs=epochs, max_patience=50,
             embedding_save_filepath="./embedding.csv",
             min_stop=100, random_seed=random_seed, gpu=0, regularization_acceleration=True, edge_subset_sz=1000000)

    sf.segmentation(domain_label_save_filepath="/home/zw/stproject/Data/V1_Breast_Cancer_Block_A_Section_1/domains.csv",
                    n_neighbors=50, resolution=resolution)  # 0.5
    domian = pd.read_csv("/home/zw/stproject/Data/V1_Breast_Cancer_Block_A_Section_1/domains.csv", index_col=None,
                         header=None)
    emb = pd.read_csv("./embedding.csv", index_col=None, header=None, sep="	", )
    sf.plot_segmentation(segmentation_figure_save_filepath="./domain_segmentation.pdf", colormap="tab20",
                         scatter_sz=1., rsz=4., csz=4., wspace=.4, hspace=.5, left=0.125, right=0.9, bottom=0.1,
                         top=0.9)
    adata.obs['SpaceFlow'] = domian.values
    adata.obs['SpaceFlow'] = adata.obs['SpaceFlow'].astype(str)
    adata.obsm["emb"] = emb.values
    sc.pp.neighbors(adata, use_rep='emb')
    sc.tl.umap(adata)
    plt.rcParams["figure.figsize"] = (3, 3)
    emb = sf.train(spatial_regularization_strength=0.1, z_dim=50, lr=1e-3, epochs=epochs, max_patience=50,
                   embedding_save_filepath="./embedding.csv",
                   min_stop=100, random_seed=random_seed, gpu=0, regularization_acceleration=True, edge_subset_sz=1000000)

    return adata.obs['SpaceFlow'], emb

def run_STAGATE(input_adata, n_cluster ,rad_cutoff=50, n_top_genes=3000,random_seed=4):
    import warnings
    warnings.filterwarnings("ignore")
    import consensus.STAGATE as STAGATE#
    import scanpy as sc
    import tensorflow as tf

    tf.compat.v1.disable_eager_execution()
    adata = input_adata.copy()
    sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=n_top_genes)
    sc.pp.normalize_total(adata, target_sum=1e4)
    STAGATE.Cal_Spatial_Net(adata, rad_cutoff=rad_cutoff)
    STAGATE.Stats_Spatial_Net(adata)
    adata = STAGATE.train_STAGATE(adata,  alpha=0, random_seed=random_seed)#
    sc.pp.neighbors(adata, use_rep='STAGATE')
    sc.tl.umap(adata)
    adata = STAGATE.mclust_R(adata, used_obsm='STAGATE', num_cluster=n_cluster)

    return adata.obs['mclust'], adata.obsm['STAGATE']

def run_stGCL(input_adata, n_cluster ,seed = 9,radius = 50, top_genes = 3000, epoch = 1000,rad_cutoff=50):
    import warnings
    warnings.filterwarnings("ignore")
    import scanpy as sc
    import os

    from consensus.stGCL.process import prefilter_genes, prefilter_specialgenes, set_seed, refine_nearest_labels#
    from consensus.stGCL import utils, train_model#
    import ot
    os.environ['R_HOME'] = '/home/zw/software/miniconda3/envs/stpython/lib/R'
    use_image = False
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

def run_stLearn(input_adata, n_cluster,random_state=0):
    import stlearn as st

    st.settings.set_figure_params(dpi=180)
    adata = input_adata.copy()
    spatial = adata.obsm['spatial']
    adata.obs["x_pixel"] = spatial[:, 0]
    adata.obs["y_pixel"] = spatial[:, 1]
    # adata.obs["x_pixel"] = adata.obs["x_um"]
    # adata.obs["y_pixel"] = adata.obs["y_um"]
    adata.var_names_make_unique()
    st.pp.filter_genes(adata, min_cells=3)
    st.pp.normalize_total(adata)
    st.pp.log1p(adata)
    st.pp.scale(adata)
    st.em.run_pca(adata, n_comps=15)
    st.pp.neighbors(adata, n_neighbors=25)
    st.tl.clustering.kmeans(adata, n_clusters=n_cluster, use_data="X_pca", key_added="X_pca_kmeans", random_state=random_state)

    return adata.obs['X_pca_kmeans'], adata.obsm['X_pca']

def run_SpaGCN (input_adata, p=0.5,res=0.5,seed=666,max_epochs=200):
    import scanpy as sc
    import numpy as np
    import consensus.SpaGCN as spg#
    import random, torch
    import warnings

    warnings.filterwarnings("ignore")
    adata = input_adata.copy()
    spatial = adata.obsm['spatial']
    adata.obs["x_pixel"] = spatial[:, 0]
    adata.obs["y_pixel"] = spatial[:, 1]
    adata.obs["x_pixel"] = spatial[:, 0]
    adata.obs["y_pixel"] = spatial[:, 1]
    x_pixel = adata.obs["x_pixel"].tolist()
    y_pixel = adata.obs["y_pixel"].tolist()
    s = 1
    b = 49
    adj = spg.calculate_adj_matrix(x=x_pixel, y=y_pixel, x_pixel=x_pixel, y_pixel=y_pixel, beta=b, alpha=s,
                                   histology=False)
    adata.var_names_make_unique()
    spg.prefilter_genes(adata, min_cells=3)  # avoiding all genes are zeros
    spg.prefilter_specialgenes(adata)
    sc.pp.normalize_per_cell(adata)
    sc.pp.log1p(adata)
    p = p
    l = spg.search_l(p, adj, start=0.01, end=1000, tol=0.01, max_run=100)
    r_seed = t_seed = n_seed = seed
    res = res
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



# 定义一个调用子函数并返回结果的包装函数，同时记录运行时间和内存使用
def call_func_with_params(func, params, monitor_performance): #func是调用子函数
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


# 定义总函数，接受一个包含子函数及其参数的字典，并汇总所有返回值
def call_functions_and_collect_results(func_dict, use_multithreading=False, monitor_performance=False):
    results = {}
    if use_multithreading:
        with concurrent.futures.ThreadPoolExecutor() as executor:
            future_to_func = {
                executor.submit(call_func_with_params, func, params, monitor_performance): func
                for func, params in func_dict.items()
            }
            for future in concurrent.futures.as_completed(future_to_func):
                func_name, result, performance_data = future.result() #返回执行结果
                results[func_name] = {
                    'result': result,
                    'performance_data': performance_data,
                    'cluster_result': result[0],
                    'emb': result[1]
                }
    else:
        for func, params in func_dict.items():
            func_name, result, performance_data = call_func_with_params(func, params, monitor_performance)
            results[func_name] = {
                'result': result,
                'performance_data': performance_data,
                'cluster_result': result[0],
                'emb':result[1]
            }
    return results

# 输出对比方法的潜在特征
def umap(func_dict):
    for func, params in func_dict.items():
        func_name, result, performance_data = call_func_with_params(func, params, False)
        print(f"{func_name} embbeding : ")
        print(result[1])
