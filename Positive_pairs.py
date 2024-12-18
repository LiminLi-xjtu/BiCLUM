import numpy as np
from data_loader import *
from scDML.calculate_NN import mnn
import pandas as pd
import networkx as nx
from itertools import chain
import os

def get_mnn(data_matrix, omic_index, k=5, approx=True, approx_method="hnswlib",
                 metric="cosine", flag="in"):
    """
    This function is modified from the original code found at:
    https://github.com/eleozzr/scDML
    Original author: xiaokangyu
    """

    cell_names = np.array(range(len(data_matrix)))

    batch_unique = np.unique(omic_index)
    cells_batch = []
    for i in batch_unique:
        cells_batch.append(cell_names[omic_index == i])

    ref = list(cells_batch[0])
    target = list(cells_batch[1])
    ds1 = data_matrix[target]
    ds2 = data_matrix[ref]
    names1 = target
    names2 = ref
    match = mnn(ds1, ds2, names1, names2, knn=k, approx=approx, approx_method=approx_method,
                metric=metric, flag=flag)
    mnns = match

    return mnns


def mnn_pair(emb_matrix, omic_index, K=5,metric="cosine",
            num1=None, num2=None, flag="in"):

    """
    obtain MNN pairs
    """

    mnn_inter_batch_approx = get_mnn(data_matrix=emb_matrix, omic_index=omic_index, k=K,
                                          approx=True, metric=metric, flag=flag)

    mnn_inter_batch = np.array([list(i) for i in mnn_inter_batch_approx])

    mnn_sort=np.lexsort(mnn_inter_batch[:,::-1].T)
    mnns=mnn_inter_batch[mnn_sort,:]
    # anchors = mnns[0:int(len(mnns)/2),:]
    # anchors[:,1] = anchors[:,1] - num1

    return mnns

def construct_pairs_mnn_cell(dataset, args, meta_batch):

    emb_matrix1 = np.vstack((dataset[0], dataset[1]))
    omic_index = meta_batch['omic_id']
    n1, n2 = np.shape(dataset[0])[0], np.shape(dataset[1])[0]
    mnn_cell = mnn_pair(emb_matrix1, omic_index, K=args.neighbors_mnn, metric=args.metric, flag="out", # mnn_method = args.mnn_method, type="mnn", flag="out",
                         num1=n1, num2=n2)
    return mnn_cell


def construct_pairs_mnn_fea(other_name, rna_name):

    mnn_fea = np.vstack([np.arange(len(other_name)), np.arange(len(rna_name))]).astype(np.int64).T

    return mnn_fea

