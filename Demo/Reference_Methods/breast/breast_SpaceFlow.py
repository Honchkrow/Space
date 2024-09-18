
import scanpy as sc
from SpaceFlow import SpaceFlow
import pandas as pd
from sklearn.metrics.cluster import adjusted_rand_score
from sklearn.metrics import normalized_mutual_info_score
import matplotlib.pyplot as plt
# sc.pp.filter_genes(adata, min_cells=3)
#
# sf = SpaceFlow.SpaceFlow(adata=adata)
section_id = "V1_Breast_Cancer_Block_A_Section_1"
k=20
print(section_id, k)
adata = sc.read_visium(path="../Data/V1_Breast_Cancer_Block_A_Section_1",
                       count_file='filtered_feature_bc_matrix.h5')
adata.var_names_make_unique()
sc.pp.filter_genes(adata, min_cells=3)
sf = SpaceFlow.SpaceFlow(adata=adata)

sf.preprocessing_data(n_top_genes=3000)

sf.train(spatial_regularization_strength=0.1, z_dim=50, lr=1e-3, epochs=1000, max_patience=50, embedding_save_filepath="./embedding.csv",
         min_stop=100, random_seed=9, gpu=0, regularization_acceleration=True, edge_subset_sz=1000000)
                                                                                                         # change
sf.segmentation(domain_label_save_filepath="../Data/V1_Breast_Cancer_Block_A_Section_1/domains.csv", n_neighbors=50, resolution=1.5)
domian=pd.read_csv("../Data/V1_Breast_Cancer_Block_A_Section_1/domains.csv",index_col=None,header=None)
emb=pd.read_csv("./embedding.csv",index_col=None,header=None,sep="	",)


sf.plot_segmentation(segmentation_figure_save_filepath="./domain_segmentation.pdf", colormap="tab20",
                     scatter_sz=1., rsz=4., csz=4., wspace=.4, hspace=.5, left=0.125, right=0.9, bottom=0.1, top=0.9)
Ann_df = pd.read_csv("../Data/V1_Breast_Cancer_Block_A_Section_1/metadata.tsv", sep="	", header=0, na_filter=False,
                     index_col=0)

# domian.columns=['SpaceFlow']
adata.obs['Ground Truth'] = Ann_df.loc[adata.obs_names, 'fine_annot_type']
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
sc.pl.spatial(adata, color="SpaceFlow",save="spaceflow_breast",title='SpaceFlow (ARI=%.2f)' % ARI)

adata.obsm["emb"]=emb.values
sc.pp.neighbors(adata, use_rep='emb')
sc.tl.umap(adata)
plt.rcParams["figure.figsize"] = (3, 3)
sc.pl.umap(adata, color="SpaceFlow", title='SpaceFlow (ARI=%.2f)'%ARI,save="SpaceFlow_breast")

# sf.pseudo_Spatiotemporal_Map(pSM_values_save_filepath="./pSM_values.tsv", n_neighbors=20, resolution=1.0)
#
# sf.plot_pSM(pSM_figure_save_filepath="./pseudo-Spatiotemporal-Map.pdf", colormap="roma", scatter_sz=1.,
#             rsz=4., csz=4., wspace=.4, hspace=.5, left=0.125, right=0.9, bottom=0.1, top=0.9)

print('%.5f ---------------------------------------------------------' %ARI)
# adata.obs.to_csv("../results/breast/domains/spaceflow_9_breast_%.2f.csv"%ARI)                                 # change