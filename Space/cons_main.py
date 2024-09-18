
def Space(input_adata,gt,k,func_dict,epochs = 10,seed = 666,alpha = 1,learning_rate = 0.01):

    from consensus.cons_func import call_functions_and_collect_results


    collected_results_mt = call_functions_and_collect_results(func_dict, use_multithreading=False,
                                                              monitor_performance=True)

    import pandas as pd
    from sklearn.metrics import adjusted_rand_score
    from consensus import consensus
    from sklearn.cluster import SpectralClustering
    from consensus.utils import calculate_location_adj,plot_results_ari,get_bool_martix

    combined_df = pd.DataFrame()
    adata = input_adata.copy()
    adata.obs['Ground Truth'] = gt
    GT = adata.obs['Ground Truth']
    GT = GT.to_frame()
    combined_df = pd.concat([combined_df, GT], axis=1)
    for func_name, data in collected_results_mt.items():
        result = data['cluster_result']
        df = result.to_frame()
        df.columns.values[0] = func_name
        combined_df = pd.concat([combined_df, df], axis=1)
    mul_reults = combined_df

    mul_reults = mul_reults.loc[mul_reults["Ground Truth"] != 'outside_VISp']
    mul_reults = mul_reults.iloc[:, 1:].dropna()
    plot_results_ari(mul_reults)
    matrices_dict = get_bool_martix(mul_reults)
    outside_VISp_indices = adata.obs['Ground Truth'] == 'outside_VISp'
    adata = adata[~outside_VISp_indices]
    posistion = pd.DataFrame(adata.obsm['spatial'])
    x_pixel = posistion[0].tolist()
    y_pixel = posistion[1].tolist()
    pos_similarity = calculate_location_adj(x=x_pixel, y=y_pixel, l=123)
    model = consensus.consensus(matrices_dict, pos_similarity, epochs=epochs, gt=gt.values, k=k,
                                seed=seed, alpha=alpha, beta=1, learning_rate=learning_rate)
    con_martix = model.train()
    sc = SpectralClustering(n_clusters=k, affinity='precomputed', random_state=666)
    labels = sc.fit_predict(con_martix)
    adata.obs["consensus"] = labels
    ari = adjusted_rand_score(labels, gt.values)

    return adata.obs["consensus"],ari
