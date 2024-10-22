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
def calculate_location_adj(spatial,l):
    print("Calculateing location adj matrix")
    posistion = pd.DataFrame(spatial)
    x = posistion[0].tolist()
    y = posistion[1].tolist()
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
                ari_matrix.iloc[j, i] = ari
    plt.figure(figsize=(10, 8))
    sns.heatmap(ari_matrix, annot=True, cmap='coolwarm', square=True)
    plt.title('ARI')
    plt.show()
    print(ari_matrix.mean())


def plot_ari_with_removal(mul_reults, remove_lowest_mean=0):
    methods = mul_reults.columns
    ari_matrix = pd.DataFrame(np.zeros((len(methods), len(methods))), index=methods, columns=methods)

    # 计算 ARI 矩阵
    for i in range(len(methods)):
        for j in range(len(methods)):
            if i == j:
                # 同一方法的ARI自然为1
                ari_matrix.iloc[i, j] = 1
            elif i < j:
                # 计算不同方法的ARI
                ari = adjusted_rand_score(mul_reults[methods[i]], mul_reults[methods[j]])
                ari_matrix.iloc[i, j] = ari
                ari_matrix.iloc[j, i] = ari

    # 绘制热图
    plt.figure(figsize=(10, 8))
    sns.heatmap(ari_matrix, annot=True, cmap='coolwarm', square=True)
    plt.title('ARI')
    plt.show()
    print(ari_matrix)
    # 输出 ARI 均值
    print(ari_matrix.mean())

    # 如果 remove_lowest_mean > 0，则删除ARI均值最小的行和列
    if remove_lowest_mean > 0:
        mean_ari = ari_matrix.mean()
        # 找到 mean 值最小的前 remove_lowest_mean 个方法
        min_mean_methods = mean_ari.nsmallest(remove_lowest_mean).index.tolist()
        print(f"Removing methods with lowest ARI means: {min_mean_methods}")
        # 从 mul_reults 和 ari_matrix 中删除对应的行和列
        mul_reults = mul_reults.drop(columns=min_mean_methods)
        # ari_matrix = ari_matrix.drop(index=min_mean_methods, columns=min_mean_methods)
    return mul_reults



def get_bool_martix(mul_reults):
    m = len(mul_reults)
    methods = mul_reults.columns
    matrices_dict = {}
    for method in methods:
        clustering_results = mul_reults[method].values
        matrix = (clustering_results[:, None] == clustering_results).astype(int)
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

