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
pip install git+https://github.com/LiminLi-XJTU/BiCLUM.git
```

## Example Usage

The input data for **BiCLUM** should be in the `.h5ad` format. If your data is in a different format, refer to the [Scanpy](https://scanpy.readthedocs.io/en/stable/) or [anndata](https://anndata.readthedocs.io/en/latest/) tutorials for instructions on how to convert your data into the `.h5ad` format.

## Example Usage

The input data for BiCLUM should be in the .h5ad format. If your data is in a different format, refer to the [Scanpy](https://scanpy.readthedocs.io/en/stable/) or [anndata](https://anndata.readthedocs.io/en/stable/) tutorials for instructions on how to convert your data into the .h5ad format.

For scATAC data, the gene activity score matrix can be generated using methods like like [ArchR](https://www.archrproject.com/) or [Signac](https://stuartlab.org/signac/). The transformation code for the scATAC modality in the **PBMC** dataset is available in the `transform` folder.


### Datasets

- **BMMC_s1d1** and **PBMC**: Paired multi-omics data (scRNA-seq and scATAC-seq).
- **Kidney**: Unpaired dataset (snRNA-seq and snATAC-seq).
- **BMCITE**: CITE-seq dataset (scRNA-seq and protein modalities), processed into two batches: `BMCITE_s1d1_s1d2` and `BMCITE_s1d2_s3d7`.

You can access these datasets on [Zenodo](https://zenodo.org/uploads/14506611).


### Running BiCLUM

To run **BiCLUM** on different types of multi-omics data, use the following commands:

#### For scRNA and scATAC integration:

```bash
python Run_RNA_ATAC.py dataset_name
```

For example, to run the **Kidney** dataset:

```bash
python Run_RNA_ATAC.py Kidney
```

#### For scRNA and Protein integration (CITE-seq):

```bash
python Run_BMCITE.py dataset_name
```

For example, to run the **BMCITE_s1d1_s1d2** dataset:

```bash
python Run_BMCITE.py BMCITE_s1d1_s1d2
```




