# Space

Reconciling Multiple Spatial Domain Identification Algorithms via Consensus Clustering

## 1. Introduction

**Space** is a spatial domain identification method from <u>**spa**</u>tially resolved transcriptomics (SRT) data using <u>**c**</u>onsensus clust<u>**e**</u>ring. It integrates **10 SOTA algorithms**. Space selects reliable algorithms by measuring their consistency. Then, it constructs a consensus matrix to integrate the outputs from multiple algorithms. We introduce **similarity loss**, **spatial loss**, and **low-rank loss** in Space to enhance accuracy and optimize computational efficiency.

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

```Shell
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

```Shell
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

```Shell
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

To facilitate user access, we have uploaded these four datasets to **[Google Drive](https://drive.google.com/drive/folders/1rXn5_HYpFo514hQXepZnaNJAD5xq9o4Z?usp=drive_link)**. Users can directly download and use them.
