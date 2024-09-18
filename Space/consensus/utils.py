import os
import torch
import pandas as pd
from sklearn.decomposition import PCA
import numba
from sklearn.metrics import adjusted_rand_score
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import ot
import random
from torch.backends import cudnn
import torch

@numba.njit("f4(f4[:], f4[:])")
def euclid_dist(t1,t2):
    sum=0
    for i in range(t1.shape[0]):
        sum+=(t1[i]-t2[i])**2
    return np.sqrt(sum)
@numba.njit("f4[:,:](f4[:,:])", parallel=True, nogil=True)
def pairwise_distance(X):
    n=X.shape[0]
    adj=np.empty((n, n), dtype=np.float32)
    for i in numba.prange(n):
        for j in numba.prange(n):
            adj[i][j]=euclid_dist(X[i], X[j])
    return adj
def calculate_location_adj(x, y,l):
    #x,y,x_pixel, y_pixel are lists
    print("Calculateing location adj matrix")
    X=np.array([x, y]).T.astype(np.float32)
    adj=pairwise_distance(X)
    adj_exp=np.exp(-1*(adj**2)/(2*(l**2)))
    return adj_exp

def plot_results_ari(mul_reults):
    methods = mul_reults.columns
    ari_matrix = pd.DataFrame(np.zeros((len(methods), len(methods))), index=methods, columns=methods)
    for i in range(len(methods)):
        for j in range(len(methods)):
            if i == j:
                # 同一方法的ARI自然为1
                ari_matrix.iloc[i, j] = 1
            elif i < j:
                # 计算不同方法的ARI
                ari = adjusted_rand_score(mul_reults[methods[i]], mul_reults[methods[j]])
                ari_matrix.iloc[i, j] = ari
                ari_matrix.iloc[j, i] = ari  # ARI是对称的

    # 绘制热图
    plt.figure(figsize=(10, 8))
    sns.heatmap(ari_matrix, annot=True, cmap='coolwarm', square=True)
    plt.title('ARI')
    plt.show()
    print(ari_matrix.mean())

def get_bool_martix(mul_reults):
    m = len(mul_reults)  # 文件的数量
    methods = mul_reults.columns  # 聚类方法的名称
    matrices_dict = {}
    for method in methods:
        clustering_results = mul_reults[method].values
        matrix = (clustering_results[:, None] == clustering_results).astype(int)
        # 保存当前方法的矩阵
        matrices_dict[method] = matrix
    return matrices_dict


def dfs(matrix, visited, vertex, cluster_label, cluster_matrix):
    visited[vertex] = True
    cluster_matrix[vertex] = cluster_label

    for i in range(len(matrix)):
        if matrix[vertex][i] == 1 and not visited[i]:
            dfs(matrix, visited, i, cluster_label, cluster_matrix)


def bool_matrix_to_clusters(bool_matrix):
    n = bool_matrix.shape[0]

    # 创建一个n*1的矩阵用于存储聚类结果
    cluster_matrix = np.zeros((n, 1), dtype=int)
    visited = np.zeros(n, dtype=bool)

    cluster_label = 1
    for i in range(n):
        if not visited[i]:
            dfs(bool_matrix, visited, i, cluster_label, cluster_matrix)
            cluster_label += 1

    return cluster_matrix.flatten()

def refine_label(adata, radius=50, key='cluster'):
    n_neigh = radius
    new_type = []
    old_type = adata.obs[key].values

    # calculate distance
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
    # adata.obs['label_refined'] = np.array(new_type)

    return new_type
def set_seed(seed=0):
    os.environ['PYTHONHASHSEED'] = str(seed)
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    cudnn.deterministic = True
    torch.backends.cudnn.deterministic = True
    cudnn.benchmark = False
    os.environ['PYTHONHASHSEED'] = str(seed)
    os.environ['CUBLAS_WORKSPACE_CONFIG'] = ':4096:8'