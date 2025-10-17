# BiCLUM

**BiCLUM** is a method for integrating multi-omics data, including scRNA-seq, scATAC-seq, and CITE-seq datasets. This repository contains code and examples for running BiCLUM on real multi-omics datasets.

![image](https://github.com/LiminLi-xjtu/BiCLUM/blob/master/img/BiCLUM_arch.jpg)

## System Requirements

BiCLUM is implemented in **Python 3.8.19**. All results presented in the paper were obtained using an **NVIDIA GeForce RTX 3090** GPU.


## Requirements

To run BiCLUM, you need to install the following dependencies:

- `anndata==0.11.1`
- `annoy==1.17.3`
- `fbpca==1.0`
- `hnswlib==0.8.0`
- `joblib==1.4.2`
- `leidenalg==0.10.2`
- `matplotlib==3.7.5`
- `munkres==1.1.4`
- `networkx==3.1`
- `numpy==1.22.3`
- `pandas==2.0.3`
- `PyYAML==6.0.2`
- `scanpy==1.10.4`
- `scikit-learn==1.3.2`
- `scipy==1.14.1`
- `seaborn==0.13.2`
- `torch==2.3.0`
- `tqdm==4.66.5`


## Example Usage

The input data for BiCLUM should be in the .h5ad format. If your data is in a different format, refer to the [Scanpy](https://scanpy.readthedocs.io/en/stable/) or [anndata](https://anndata.readthedocs.io/en/stable/) tutorials for instructions on how to convert your data into the .h5ad format.

For scATAC data, the gene activity score matrix can be generated using methods like like [ArchR](https://www.archrproject.com/) or [Signac](https://stuartlab.org/signac/). The transformation code for the scATAC modality in the **PBMC (paired)** dataset is available in the `transform` folder.


### Datasets

- **BMMC (paired)** and **PBMC (paired)**: Paired multi-omics data consisting of **scRNA-seq** and **scATAC-seq**.
- **BMMC (unpaired)** and **PBMC (unpaired)**: Unpaired multi-omics data consisting of **scRNA-seq** and **scATAC-seq**.
- **Kidney**: Unpaired dataset combining **snRNA-seq** and **snATAC-seq** data.
- **BMCITE**: CITE-seq dataset with **scRNA-seq** and **protein modalities**, processed into two batches: `BMCITE_s1d1_s1d2` and `BMCITE_s1d2_s3d7`.

#### Data Sources

- The **BMMC** and **BMCITE** datasets (including `BMCITE_s1d1_s1d2` and `BMCITE_s1d2_s3d7`) were obtained from the work of [Luecken, Malte D., et al.](https://datasets-benchmarks-proceedings.neurips.cc/paper/2021/hash/158f3069a435b314a80bdcb024f8e422-Abstract-round2.html) and downloaded from the **Gene Expression Omnibus (GEO)** repository (accession number: [GSE194122](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE194122)).

- The **PBMC (paired)** dataset was obtained from **10x Genomics** ([PBMC (paired) dataset link](https://www.10xgenomics.com/resources/datasets/pbmc-from-a-healthy-donor-granulocytes-removed-through-cell-sorting-3-k-1-standard-2-0-0)).
- The **PBMC (unpaired)** dataset was obtained from the work of [Wang C, et al.](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02116-x) and downloaded from ([PBMC (unpaied) dataset link](https://github.com/liulab-dfci/MAESTRO/tree/master/data)).

- The **Kidney** dataset was obtained from the work of [Muto, Yoshiharu, et al.](https://www.nature.com/articles/s41467-021-22368-w) and retrieved from the **GEO** repository (accession number: [GSE151302](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE151302)).

You can access these multi-omics datasets and their corresponding gene activity score matrices on [Zenodo](https://zenodo.org/uploads/14506611).

### Configuration

The configuration file contains the parameter settings and data information for different datasets. You can find it in the `config/` folder. The settings include important hyperparameters such as learning rates, batch sizes, and specific model parameters tailored for each dataset, as well as the paths to the input datasets. You can modify these settings according to your requirements.

Example parameter settings in the config file:

```bash
# Dataset and transformation settings
dataset_name: BMMC_paired  # Name of the dataset (e.g., BMMC (paired)). Specify which dataset you are working with.
dataset_dir: ../data/BMMC_paired  # Path to the BMMC (paired) data directory. This directory should include files like rna.h5ad, atac.h5ad, and scGAM_ArchR.h5ad.
                          # These datasets can be downloaded from [Zenodo](https://zenodo.org/uploads/14506611).
GAM_name: ArchR  # Transformation method used to convert the ATAC-seq data (atac.h5ad) into a gene activity score matrix.
                 # Options could include 'ArchR' or 'Signac', depending on your preprocessing method.
dataset_type: RNA_ATAC  # Type of data integration. This defines which modalities you are integrating. 
                       # Options could be:
                       # - 'RNA_ATAC' for RNA and ATAC-seq integration
                       # - 'RNA_Protein' for RNA and Protein integration
batch: None  # Batch information for BMCITE (CITE-seq) data. If no batch information is available, set this parameter to None.
paired: True  # Whether the data is paired (True/False). 'True' if the data are from paired modalities (e.g., scRNA-seq and scATAC-seq), 
              # 'False' otherwise (e.g., unpaired datasets).

# Model hyperparameters
n_high_var: 2000  # The number of highly variable genes (HVGs) selected for preprocessing. 
dim: 100  # The dimensionality of PCA embeddings in the preprocessing step. It reduces the number of features while retaining variance.
neighbors_mnn: 500  # The number of nearest neighbors used for Mutual Nearest Neighbor (MNN) construction. 
metric: cosine  # The distance metric to use for MNN construction. 
                # You can also choose other distance metrics like 'euclidean' depending on the nature of your data.
use_rep: hvg_count  # The representation used for gene features. Options include:
                   # - 'hvg_count' for using raw counts of high-variance genes
                   # - 'hvg_norm' for using normalized data
                   # - 'low_emb' for using low-dimensional embeddings from previous steps as feature representations.

latent_dim: 50  # The dimensionality of the latent space for the model. 

# Regularization parameters
tau_cell: 0.5  # Hyperparameter for the cell-level contrastive loss function, default is 0.5. 
tau_feature: 0.5  # Hyperparameter for the feature-level contrastive loss function, default is 0.5.
gamma_a: 1 # The weight for modality 1 reconstruction regularization, default is 1.
gamma_b: 1 # The weight for modality 2 reconstruction regularization, default is 1.
alpha: 10000  # The weight for cell-level regularization, default is 10000. 
beta: 10000  # The weight for feature-level regularization, default is 10000. 

# Training settings
seed: 123  # The random seed for reproducibility
batch_size: 256  # The batch size for training
learning_rate: 0.0001  # The learning rate for optimization, default is 0.0001
weight_decay: 0.00005  # The weight decay (L2 regularization) applied to model parameters to prevent overfitting and encourage simpler models, deafult is 0.00005.
epoch: 1000  # The number of training epochs

```

### Running BiCLUM

To run **BiCLUM** for multi-omics data integration, follow these steps depending on the type of datasets you are working with.

#### 1. **For scRNA and scATAC Integration:**

To integrate **scRNA-seq** and **scATAC-seq** data, use the following command:

```bash
python Run_RNA_ATAC.py dataset_name
```

Where `dataset_name` is the name of your dataset (e.g., **PBMC_paired**, **Kidney**, etc.).

**Example**: To run the **PBMC_paired** dataset (scRNA-seq and scATAC-seq):

```bash
python Run_RNA_ATAC.py PBMC_paired
```

#### 2. **For scRNA and Protein Integration (CITE-seq):**

To integrate **scRNA-seq** and **protein** data (CITE-seq), use the following command:

```bash
python Run_BMCITE.py dataset_name
```

**Example**: To run the **BMCITE_s1d1_s1d2** dataset (scRNA-seq and Protein):

```bash
python Run_BMCITE.py BMCITE_s1d1_s1d2
```

### Outputs

1. **Model Outputs**: The parameter file after model training is saved in the `output/dataset_name/` folder.
2. **Embeddings and Losses**: The cell embeddings (`inte_cell`) and feature embeddings (`inte_fea`) for multimodal data, along with the model's loss list (`loss_list`), are saved in the `./results/dataset_name/GAM_name.npy` file.
3. **Model Evaluation**: The indicators for evaluating model performance, including **omics mixing**, **cell type preservation**, **transfer accuracy**, and **FOSCTTM**, are saved in the `eva/dataset_name/eva_GAM_name.csv` file. 

### Visualizing the Results

Once the model is trained and the embeddings are saved, you can visualize the results with the following commands:

#### 1. **For scRNA and scATAC Integration:**

```bash
python Vis_RNA_ATAC.py dataset_name
```

#### 2. **For scRNA and Protein Integration (CITE-seq):**

```bash
python Vis_BMCITE.py dataset_name
```

Replace `dataset_name` with the name of your dataset (e.g., **PBMC_paired**, **BMCITE_s1d1_s1d2**).

The visualization results, including plots such as UMAP or PAGA graph, are saved in the `vis/dataset_name/` folder.

---

### Notes:
- **output folder**: Contains the model parameters after training.
- **results folder**: Stores the embeddings of cells and features, as well as the loss list.
- **eva folder**: Includes model performance indicators for evaluating the success of the integration.
- **vis folder**: Holds the generated visualizations from the trained model.
