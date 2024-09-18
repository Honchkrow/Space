# Space

Reconciling Multiple Spatial Domain Identification Algorithms via Consensus Clustering

## 1. Introduction

**Space** is a spatial domain identification method from <u>spa</u>tially resolved transcriptomics (SRT) data using <u>c</u>onsensus clust<u>e</u>ring. It integrates **10 SOTA algorithms**. Space selects reliable algorithms by measuring their consistency. Then, it constructs a consensus matrix to integrate the outputs from multiple algorithms. We introduce **similarity loss**, **spatial loss**, and **low-rank loss** in Space to enhance accuracy and optimize computational efficiency.

<center>
    <img style="border-radius: 0.3125em;
    box-shadow: 0 2px 4px 0 rgba(34,36,38,.12),0 2px 10px 0 rgba(34,36,38,.08);" 
    src="./Images/main_figure.svg">
    <br>
    <div style="color:orange; border-bottom: 1px solid #d9d9d9;
    display: inline-block;
    color: #999;
    padding: 2px;">Space Workflow</div>
</center>

The integrated methods:
- [x] **GraphST**
- [x] **Leiden**
- [x] **MENDER**
- [x] **Louvain**
- [x] **SEDR**
- [x] **SpaceFlow**
- [x] **SpaGCN**
- [x] **STAGATE**
- [x] **stGCL**
- [x] **stLearn**


The **organization** of this repository file is as follows:

```shell
Space/  
├── Data/                    # data for reproducibility
├── Demo/  
│   ├── Reference_Methods/   # the scripts of 10 SOTA algorithms
│   └── Reproduce_Scripts/   # the scripts for reproducing results in manuscript
├── Images/                  # images
├── Space/                   # Space source code
└── environment.yml          # conda/mamba environment file
```


## 2. Installation Tutorial

The deployment of Space requires a Linux/Unix machine. We recommend using [conda](https://anaconda.org/anaconda/conda)/[mamba](https://github.com/conda-forge/miniforge) and create a virtual environment to manage all the dependencies. If you did not install conda before, please install [conda](https://anaconda.org/anaconda/conda)/[mamba](https://github.com/conda-forge/miniforge) first.

We provide the environment file, allowing users to quickly deploy Space using the following command.

```shell
# clone or download this repository
git clone https://github.com/Honchkrow/Space

# enter the folder
cd Space

# install environment using environment.yml
conda env create -n Space -f environment.yml
```

*<font color=red>Note:</font> The environment name can be changed by replacing "-n Space" to "-n environment_name".*

*<font color=red>Note:</font> If errors about **<font color=blue>unavailable or invalid channel</font>** occur, please check that whether the **.condarc** file in your ~ directory had been modified. Modifing .condarc file may cause wrong channel error. In this case, just rename/backup your .condarc file. Once the installation finished, this file can be recoveried. Of course, you can delete .condarc file if necessary.*

Once the environment is created, the users can enter the environment using the following command.

```shell
conda activate Space
```

## 3. How to use Space

In this section, we will use a SRT dataset to provide a detailed introduction to the functionalities of Space.

### 3.1 Datasets

In the manuscript for Space, we present the results of Space on four different datasets. These datasets are:

- [Human breast cancer](https://support.10xgenomics.com/spatial-gene-expression/datasets/1.1.0/V1_Breast_Cancer_Block_A_Section_1): This dataset contains 3,798 spots and 36,601 genes, along with 20 manually annotated regions.
- [Mouse hypothalamus](https://datadryad.org/stash/dataset/doi:10.5061/dryad.8t8s248): This dataset includes 5 tissue sections, with 8 regions manually annotated, and the number of spots varies from 5,488 to 5,926.
- [Mouse primary visual area](https://spacetx.github.io/data.html): This dataset comprises three slices, containing 3390, 4491, and 3545 spots, respectively, for a total of 79 genes. Manual annotation was performed on six visual cortex layers, ranging from VISP_I to VISP_VI, as well as the white matter region (VISpwm).
- [Mouse visual cortex](https://www.starmapresources.com/data): This dataset includes three tissue sections_BZ5, BZ9, and BZ14_which have spot counts of 1049, 1053, and 1088, respectively, amounting to a total of 166 genes. Manual annotation was conducted for four regions.

We have prepared two types of data. The first type consists of results on four datasets, obtained using ten SOTA algorithms. This data is in CSV format, which allows users to quickly reproduce the results of the Space article. **This data is already integrated into this repository, so users do not need to download it separately.** The second type is the processed SRT data, which includes gene expression matrices, spatial location information, and H&E images. Due to the large size of this data, it cannot be uploaded to GitHub. Therefore, users will need to download it.

#### Download the processed SRT datasets (not mandatory)

To facilitate user access, we have uploaded the processed SRT datasets to **[Google Drive](https://drive.google.com/drive/folders/1rXn5_HYpFo514hQXepZnaNJAD5xq9o4Z?usp=drive_link)** and **[BaiduYun](https://pan.baidu.com/s/1qxoq0ttp0BzsvLzBtDUgIg?pwd=3vvv)**. Users can directly download and use them.

To facilitate users in quickly reproducing our results, they can merge the extracted 'Data' folder with the 'Data' folder in the Space project. This can be done immediately after downloading and unzipping the files.The **organization** of this project will become:

```Shell
# Only shows the BARISTASeq dataset.
# Mouse_hippocampus_MERFISH, SRARmap_pa and V1_Breast_Cancer_Block_A_Section_1 are the same.
Space/  
├── Data/
│   ├── BARISTASeq/
│   │   ├── BARISTASeq_Sun2021Integrating_Slice_1_data.h5ad
│   │   ├── BARISTASeq_Sun2021Integrating_Slice_2_data.h5ad
│   │   ├── BARISTASeq_Sun2021Integrating_Slice_3_data.h5ad
│   │   ├── result1.csv
│   │   ├── result2.csv
│   │   └── result3.csv
│   ├── Mouse_hippocampus_MERFISH/  # not show
│   ├── SRARmap_pa/  # not show
│   └── V1_Breast_Cancer_Block_A_Section_1/  # not show
├── Demo/  
│   ├── Reference_Methods/
│   └── Reproduce_Scripts/
├── Images/
├── Space/
└── environment.yml
```


### 3.2 Reproducing the results of Space (simple tutorial)

To reproduce the results of the Space article, users can run the scripts in the **Demo** folder. The scripts are organized into two folders: **Reference_Methods** and **Reproduce_Scripts**. The **Reference_Methods** folder contains scripts for reproducing the results of the ten SOTA algorithms. The **Reproduce_Scripts** folder contains scripts for reproducing the results of the Space.

#### How to integrate the results of different algorithms using Space

Here, for quick illustration, we directly apply Space to the results obtained from 10 SOTA methods. These methods have already been executed. The scripts are asved in **Reference_Methods** folder. The results of these methods are saved in the **Data** folder.

First, load the necessary packages and set R environment. 

*<font color=red>Please note that in the code below, the R environment must be the one installed within Space. Users need to replace it according to the installation directory of Space.</font>*

```python
import os
import scanpy as sc
import pandas as pd
import Space
from sklearn.metrics import adjusted_rand_score
from sklearn.cluster import SpectralClustering
from Space.cons_func import get_results, get_domains
from Space.utils import calculate_location_adj, plot_results_ari, get_bool_martix

# The mclust is used.
# Please modify this path!
os.environ["R_HOME"] = "/home/zw/software/miniforge-pypy3/envs/space/lib/R"
```

Next, load the dataset.

```python
# read the expression data
adata = sc.read_visium(
    path="./Data/V1_Breast_Cancer_Block_A_Section_1", 
    count_file="filtered_feature_bc_matrix.h5"
)

# read the metadata
Ann_df = pd.read_csv(
    "./Data/V1_Breast_Cancer_Block_A_Section_1/metadata.tsv",
    sep="\t",
    header=0,
    na_filter=False,
    index_col=0,
)
adata.var_names_make_unique()

# read the image representation
im_re = pd.read_csv(
    "./Data/V1_Breast_Cancer_Block_A_Section_1/image_representation/ViT_pca_representation.csv",
    header=0,
    index_col=0,
    sep=",",
)
adata.obsm["im_re"] = im_re
adata.obs["gt"] = Ann_df["fine_annot_type"]
gt = adata.obs["gt"]
```

Then, set the parameters.

```python
k = 20                   # number of clusters
epochs = 120             # epoch in training
seed = 666               # random seed
alpha = 1                # recommended value
learning_rate = 0.0001   # learning rate in training
```

Now, read the results from 10 SOTA methods.

```python
mul_reults = pd.read_csv(
    "./Data/V1_Breast_Cancer_Block_A_Section_1/result.csv", 
    header=0, 
    index_col=0
)
mul_reults = mul_reults.iloc[:, 2:]

```


Now, we can observe the consistency between the results of different methods.

```python
plot_results_ari(mul_reults)
```

<center>
    <img style="border-radius: 0.3125em;
    box-shadow: 0 2px 4px 0 rgba(34,36,38,.12),0 2px 10px 0 rgba(34,36,38,.08);" 
    src="./Images/consistency.png">
    <br>
    <div style="color:orange; border-bottom: 1px solid #d9d9d9;
    display: inline-block;
    color: #999;
    padding: 2px;">Consistency between different methods</div>
</center>


```python
# Discard methods that show poor consistency
mul_reults = mul_reults.drop("SpaceFlow", axis=1)
mul_reults = mul_reults.drop("MENDER", axis=1)

# compute the positional similarity matrix
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

# tarining model
con_martix = model.train()

# set spectral cluster model
sc = SpectralClustering(n_clusters=k, affinity="precomputed", random_state=666)

# clustering
labels = sc.fit_predict(con_martix)

adata.obs["consensus"] = labels

ari = adjusted_rand_score(labels, gt.values)

print(ari)
```

you will obtain a result from Space with an ARI of 0.648.

**<font color=red>In most cases, Space does not yield a fixed result. This is not due to an issue with Space, but because some methods exhibit randomness even when the random seed is fixed. Please refer to [this](https://github.com/QIFEIDKN/STAGATE/issues/10) for more information. However, the variations in the results we obtain are minimal. The outcomes are stable across multiple runs.</font>**


#### Visualization


#### Domain-specific gene analysis


#### Trajectory inference


#### Use with Scanpy or Seurat


### 3.3 How to choose and use different baseline algorithms


## 4 Citation

