
import scanpy as sc
from SpaceFlow import SpaceFlow
import pandas as pd
from sklearn.metrics.cluster import adjusted_rand_score
from sklearn.metrics import normalized_mutual_info_score
import matplotlib.pyplot as plt
# sc.pp.filter_genes(adata, min_cells=3)
#
# sf = SpaceFlow.SpaceFlow(adata=adata)
#传参所用
import sys
# seed = int(sys.argv[1])
seed = 3
adata = sc.read("../Data/SRARmap_pa/MS_raw_Dataset11_Dataset11_data.h5ad")
adata.var_names_make_unique()
adata = adata[adata.obs['slice_id'] == "BZ5", :]
adata.obs['Ground Truth'] = adata.obs['ct']
sc.pp.filter_genes(adata, min_cells=3)

sf = SpaceFlow.SpaceFlow(adata=adata)

sf.preprocessing_data(n_top_genes=3000)

sf.train(spatial_regularization_strength=0.1, z_dim=50, lr=1e-3, epochs=1000, max_patience=50, embedding_save_filepath="./embedding.csv",
         min_stop=100, random_seed= seed, gpu=0, regularization_acceleration=True, edge_subset_sz=1000000)
                                                                                                         # change
sf.segmentation(domain_label_save_filepath="../Data/V1_Breast_Cancer_Block_A_Section_1/domains.csv", n_neighbors=50, resolution=0.5)    #1.5
domian=pd.read_csv("../Data/V1_Breast_Cancer_Block_A_Section_1/domains.csv",index_col=None,header=None)
emb=pd.read_csv("./embedding.csv",index_col=None,header=None,sep="	",)


sf.plot_segmentation(segmentation_figure_save_filepath="./domain_segmentation.pdf", colormap="tab20",
                     scatter_sz=1., rsz=4., csz=4., wspace=.4, hspace=.5, left=0.125, right=0.9, bottom=0.1, top=0.9)

adata.obs['SpaceFlow'] = domian.values
adata.obs['SpaceFlow']=adata.obs['SpaceFlow'].astype(str)

obs_df = adata.obs.dropna()
# obs_df.to_csv("result/{}_type_stMMR.csv".format(section_id))
ARI = adjusted_rand_score(obs_df['SpaceFlow'], obs_df['Ground Truth'])
print('Adjusted rand index = %.5f' % ARI)
from sklearn.metrics import normalized_mutual_info_score
nmi=normalized_mutual_info_score(obs_df['SpaceFlow'], obs_df['Ground Truth'])
print('normalized mutual info score = %.5f' % nmi)
plt.rcParams["figure.figsize"] = (3, 3)
sc.pl.embedding(adata, basis="spatial", color="SpaceFlow", s=6, show=False, title='SpaceFlow',save="SpaceFlow_pa.pdf")

adata.obsm["emb"]=emb.values
sc.pp.neighbors(adata, use_rep='emb')
sc.tl.umap(adata)
plt.rcParams["figure.figsize"] = (3, 3)
sc.pl.umap(adata, color="SpaceFlow", title='SpaceFlow (ARI=%.2f)'%ARI,save="SpaceFlow_hip")

# sf.pseudo_Spatiotemporal_Map(pSM_values_save_filepath="./pSM_values.tsv", n_neighbors=20, resolution=1.0)
#
# sf.plot_pSM(pSM_figure_save_filepath="./pseudo-Spatiotemporal-Map.pdf", colormap="roma", scatter_sz=1.,
#             rsz=4., csz=4., wspace=.4, hspace=.5, left=0.125, right=0.9, bottom=0.1, top=0.9)

print('%.5f ---------------------------------------------------------' %ARI)
# adata.obs.to_csv("../results/PA/domains/spaceflow_9_PA_%.2f.csv"%ARI)                                 # change
# adata.obs.to_csv("../results/PA/domains/slice2_spaceflow_{}_PA_{:.2f}.csv".format(seed,ARI))