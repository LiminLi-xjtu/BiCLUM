
import os
import numpy as np
import torch
import pandas as pd
import anndata
import scanpy

from Positive_pairs import construct_pairs_mnn_cell, construct_pairs_mnn_fea
from data_loader import train_load, test_cell_load, test_fea_load


def load_data(dataset_type="RNA_ATAC", dataset_dir=None, GAM="ArchR", batch=None):
    """Load data """

    if dataset_type == "RNA_ATAC":
        Other_data = anndata.read(os.path.join(dataset_dir, 'GASM/scGAM_' + GAM + '.h5ad'))
        RNA_data = anndata.read(os.path.join(dataset_dir, 'GASM/rna.h5ad'))

        Other_data.obs['omic_id'] = 'Other'
        RNA_data.obs['omic_id'] = 'RNA'

        adata_fm = Other_data.concatenate(RNA_data, join='inner', batch_key='domain_id')
        RNA_data = adata_fm[adata_fm.obs["omic_id"] == 'RNA']
        Other_data = adata_fm[adata_fm.obs["omic_id"] == 'Other']

        label = []
        label.append(np.array(Other_data.obs["cell_type"]))
        label.append(np.array(RNA_data.obs["cell_type"]))
        meta_batch = adata_fm.obs

        print("number of genes=%d"%(adata_fm.shape[1]))
        return Other_data, RNA_data, label, meta_batch


    elif dataset_type == 'RNA_Protein':

        Other_data = anndata.read(os.path.join(dataset_dir, 'bm.adt.' + batch + '.h5ad'))
        RNA_data = anndata.read(os.path.join(dataset_dir, 'bm.rna.' + batch + '.h5ad'))

        Other_data.obs['omic_id'] = 'Other'
        RNA_data.obs['omic_id'] = 'RNA'

        label = []
        label.append(np.array(Other_data.obs["cell_type"]))
        label.append(np.array(RNA_data.obs["cell_type"]))
        meta_batch = pd.concat([Other_data.obs, RNA_data.obs])

        print("number of genes=%d"%(Other_data.shape[1]))
        return Other_data, RNA_data, label, meta_batch


def processing(gam, rna, data_id, gam_type):

    rna.X = rna.X.astype(int)
    rna.layers["raw"] = rna.X.copy()
    scanpy.pp.highly_variable_genes(rna, n_top_genes=2000, flavor="seurat_v3")

    if gam_type == 'snapATAC2' or gam_type == 'Signac':
        gam.X = gam.X.astype(int)
        gam.layers["raw"] = gam.X.copy()
        scanpy.pp.highly_variable_genes(gam, n_top_genes=2000, flavor="seurat_v3")
    else:
        gam.layers["raw"] = gam.X.copy()
        try:
            scanpy.pp.highly_variable_genes(gam, n_top_genes=2000)
        except:
            gam.var['highly_variable'] = rna.var['highly_variable']

    hvg_rna = rna.var['highly_variable']
    hvg_gam = gam.var['highly_variable']
    hvg = hvg_rna | hvg_gam
    gam.var['highly_variable'] = hvg
    rna.var['highly_variable'] = hvg

    scanpy.pp.normalize_total(gam)
    scanpy.pp.log1p(gam)
    scanpy.pp.scale(gam)
    scanpy.tl.pca(gam, n_comps=100, svd_solver="auto")

    scanpy.pp.normalize_total(rna)
    scanpy.pp.log1p(rna)
    scanpy.pp.scale(rna)
    scanpy.tl.pca(rna, n_comps=100, svd_solver="auto")

    path = '../data/processed_data/' + data_id + '/'
    if not os.path.exists(path):
        os.makedirs(path)

    if 'orig.ident' in rna.obs:
        rna.obs['orig.ident'] = rna.obs['orig.ident'].astype(str)
    if 'author_cell_type' in rna.obs:
        rna.obs['author_cell_type'] = rna.obs['author_cell_type'].astype(str)
    if 'artif_dupl' in rna.var:
        rna.var['artif_dupl'] = rna.var['artif_dupl'].astype(str)
    if rna.raw is not None and '_index' in rna.raw.var.columns:
        rna.raw.var.drop(columns=['_index'], inplace=True)

    if 'orig.ident' in gam.obs:
        gam.obs['orig.ident'] = gam.obs['orig.ident'].astype(str)
    if 'author_cell_type' in gam.obs:
        gam.obs['author_cell_type'] = gam.obs['author_cell_type'].astype(str)
    if 'artif_dupl' in gam.var:
        gam.var['artif_dupl'] = gam.var['artif_dupl'].astype(str)
    if gam.raw is not None and '_index' in gam.raw.var.columns:
        gam.raw.var.drop(columns=['_index'], inplace=True)

    rna.write(os.path.join(path, gam_type + "-rna-pp.h5ad"), compression="gzip")
    gam.write(os.path.join(path, gam_type + "-gam-pp.h5ad"), compression="gzip")

    return gam, rna


def data_input(args):

    print(f"load data")
    other, rna, label, meta_batch = load_data(args.dataset_type, args.dataset_dir, args.GAM_name, args.batch)

    print(f"preprocess RNA and Other data")
    if args.dataset_type == 'RNA_ATAC':
        processing(other, rna, args.dataset_name, args.GAM_name)
        # RNA_data = ad.read_h5ad(
        #     os.path.join('../data/processed_data/', args.dataset_name + '/' + args.GAM_name + "-rna-pp.h5ad"))
        # Other_data = ad.read_h5ad(
        #     os.path.join('../data/processed_data/', args.dataset_name + '/' + args.GAM_name + "-gam-pp.h5ad"))
    else:
        RNA_data = rna
        Other_data = other

    # input for encoder
    dataset = []
    if args.dataset_type == 'RNA_ATAC':
        dataset.append(Other_data.obsm['X_pca'])
        dataset.append(RNA_data.obsm['X_pca'])
    else:
        dataset.append(Other_data.obsm['X_apca'])
        dataset.append(RNA_data.obsm['X_pca'])

    # output for decoder
    if args.dataset_type == 'RNA_ATAC':
        other_hvg = Other_data.var.query("highly_variable").index
        rna_hvg = RNA_data.var.query("highly_variable").index
    else:
        fea_pairs = pd.read_csv(os.path.join(args.dataset_dir, 'graph.csv'))
        other_hvg = Other_data.var_names
        rna_hvg_1 = set(RNA_data.var_names[np.where(RNA_data.var['vst.variable']==1)])
        rna_hvg_2 = set(fea_pairs.values[:, 1])
        rna_hvg = list(rna_hvg_1.union(rna_hvg_2))

    datahvg = []
    if args.use_rep == 'hvg_count':
        Other_data = Other_data[:, other_hvg]
        RNA_data = RNA_data[:, rna_hvg]
        datahvg.append(Other_data.layers['raw'])
        datahvg.append(RNA_data.layers['raw'])
    elif args.use_rep == 'hvg_norm':
        Other_data = Other_data[:, other_hvg]
        RNA_data = RNA_data[:, rna_hvg]
        datahvg.append(Other_data.X)
        datahvg.append(RNA_data.X)
    elif args.use_rep == 'low_emb':
        datahvg = dataset

    for i in range(len(dataset)):
        dataset[i] = torch.from_numpy(dataset[i].astype('float32')).to(args.device)
        datahvg[i] = torch.from_numpy(datahvg[i].toarray().astype('float32')).to(args.device)


    m1, m2 = np.shape(dataset[0])[1], np.shape(dataset[1])[1]
    m1_hvg, m2_hvg = np.shape(datahvg[0])[1], np.shape(datahvg[1])[1]
    n1, n2 = np.shape(dataset[0])[0], np.shape(dataset[1])[0]

    args.in_features = [m1, m2]
    args.out_features = [m1_hvg, m2_hvg]
    args.in_cells = [n1, n2]

    print(f"construct positive pairs")
    if args.dataset_type == 'RNA_ATAC':
        mnn_cell = construct_pairs_mnn_cell([Other_data.X, RNA_data.X], args, meta_batch)
        mnn_fea = construct_pairs_mnn_fea(other_hvg, rna_hvg)
    else:
        other_names = fea_pairs.values[:, 2]
        rna_names = fea_pairs.values[:, 1]
        mnn_cell = construct_pairs_mnn_cell([Other_data[:,other_names].X.toarray(), RNA_data[:,rna_names].X.toarray()], args, meta_batch)
        mnn_fea = construct_pairs_mnn_fea(other_names, rna_names)

    train_loader = train_load(mnn_cell, mnn_fea, args)
    test_cell_loader = test_cell_load(args)
    test_fea_loader = test_fea_load(args)

    return dataset, datahvg, train_loader, test_cell_loader, test_fea_loader, args

