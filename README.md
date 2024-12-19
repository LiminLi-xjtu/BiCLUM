# BiCLUM

**BiCLUM** is a method for integrating multi-omics data, including scRNA-seq, scATAC-seq, and CITE-seq datasets. This repository contains code and examples for running BiCLUM on real multi-omics datasets.

![image](https://github.com/LiminLi-xjtu/BiCLUM/blob/master/img/BiCLUM_arch.jpg)


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

## Installation

You can install the development version of **BiCLUM** via pip directly from GitHub:

```bash
pip install git+https://github.com/LiminLi-xjtu/BiCLUM.git
```

## Example Usage

The input data for BiCLUM should be in the .h5ad format. If your data is in a different format, refer to the [Scanpy](https://scanpy.readthedocs.io/en/stable/) or [anndata](https://anndata.readthedocs.io/en/stable/) tutorials for instructions on how to convert your data into the .h5ad format.

For scATAC data, the gene activity score matrix can be generated using methods like like [ArchR](https://www.archrproject.com/) or [Signac](https://stuartlab.org/signac/). The transformation code for the scATAC modality in the **PBMC** dataset is available in the `transform` folder.


### Datasets

- **BMMC_s1d1** and **PBMC**: Paired multi-omics data (scRNA-seq and scATAC-seq).
- **Kidney**: Unpaired dataset (snRNA-seq and snATAC-seq).
- **BMCITE**: CITE-seq dataset (scRNA-seq and protein modalities), processed into two batches: `BMCITE_s1d1_s1d2` and `BMCITE_s1d2_s3d7`.

You can access these datasets on [Zenodo](https://zenodo.org/uploads/14506611).


### Running BiCLUM

To run **BiCLUM** for multi-omics data integration, follow these steps depending on the type of datasets you are working with.

#### 1. For scRNA and scATAC Integration, use the following command:

```bash
python Run_RNA_ATAC.py dataset_name
```

Where `dataset_name` is the name of your dataset (e.g., **PBMC**, **Kidney**, etc.).

**Example**: To run the **PBMC** dataset (scRNA-seq and scATAC-seq):

```bash
python Run_RNA_ATAC.py PBMC
```

#### 2. **For scRNA and Protein Integration (CITE-seq)**:

```bash
python Run_BMCITE.py dataset_name
```

**Example**: To run the **BMCITE_s1d1_s1d2** dataset (scRNA-seq and Protein):

```bash
python Run_BMCITE.py BMCITE_s1d1_s1d2
```


### Visualizing the Results

Use the following command to visualize the results:

#### 1. **For scRNA and scATAC Integration**:

```bash
python Vis_RNA_ATAC.py dataset_name
```

#### 2. **For scRNA and Protein Integration (CITE-seq)**:

```bash
python Vis_BMCITE.py dataset_name
```

Replace `dataset_name` with the name of the datasets (e.g., **PBMC**, **BMCITE_s1d1_s1d2**).
