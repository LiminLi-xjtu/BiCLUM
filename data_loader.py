import numpy as np
from torch.utils.data import Dataset, DataLoader
import torch
import scipy.sparse as sps
from collections import defaultdict





class GetDatasetMNN(Dataset):
    def __init__(self, anchor_mnn, args):

        self.anchor_mnn = anchor_mnn
        self.sample = args.in_cells
        self.num_samples = sum(self.sample)
        self.pos_mnn_dict, self.mnn_idx = self.mnn_pairs_dict()


    def mnn_pairs_dict(self):

        mnn_graph = sps.csr_matrix((np.ones(self.anchor_mnn.shape[0]), (self.anchor_mnn[:, 0], self.anchor_mnn[:, 1])),
                                   dtype=np.int8, shape=(self.num_samples, self.num_samples))

        mnn_dict = defaultdict(list)
        mnn_idx = []

        for i in range(self.num_samples):
            tmp_idx = mnn_graph[i, :]
            row, col = np.nonzero(tmp_idx)
            mnn_idx.append(i)
            mnn_dict[i].append(col)

        return mnn_dict, mnn_idx

    def __getitem__(self, idx):

        flag = len(self.pos_mnn_dict[idx][0]) > 0

        if idx < self.sample[0]:
            if flag:
                pos_mnn_anchors = self.pos_mnn_dict[idx][0]
                p1 = np.random.choice(pos_mnn_anchors) - self.sample[0]
            else:
                p1 = np.random.choice(self.sample[1])
            idx_other, idx_rna = idx, p1

        else:
            if flag:
                pos_mnn_anchors = self.pos_mnn_dict[idx][0]
                p0 = np.random.choice(pos_mnn_anchors)
            else:
                p0 = np.random.choice(self.sample[0])
            p1 = idx-self.sample[0]
            idx_other, idx_rna = p0, p1

        return idx_other, idx_rna, flag

    def __len__(self):
        return self.num_samples


class GetDatasetFeaturePair(Dataset):
    def __init__(self, anchor_feature, args):

        self.anchor_feature = anchor_feature
        # self.device = args.device
        self.features = args.out_features
        self.num_features = self.features[0]

        self.pos_feature_dict, self.feature_idx = self.feature_pairs_dict()

    def feature_pairs_dict(self):

        n_feature1, n_feature2 = self.features
        feature_graph = sps.csr_matrix((np.ones(self.anchor_feature.shape[0]), (self.anchor_feature[:, 0], self.anchor_feature[:, 1])),
                                   dtype=np.int8, shape=self.features)

        feature_dict = defaultdict(list)

        feature_idx = []

        for i in range(n_feature1):
            tmp_idx= feature_graph[i, :]
            row, col = np.nonzero(tmp_idx)
            feature_idx.append(i)
            feature_dict[i].append(col)

        return feature_dict, feature_idx

    def __getitem__(self, idx):

        idx_other_feature, idx_rna_feature = idx, idx

        return idx_other_feature, idx_rna_feature

    def __len__(self):
        return self.num_features



class GetDatasetModalityFeature(Dataset):
    def __init__(self, num_features):

        self.num_features = num_features

    def __getitem__(self, idx):

        return idx

    def __len__(self):
        return self.num_features


class GetDatasetModalityCell(Dataset):
    def __init__(self, num_cells):

        self.num_cells = num_cells

    def __getitem__(self, idx):

        return idx

    def __len__(self):

        return self.num_cells


def train_load(cell_mnn, feature_pair, args):

    train_dataset_MNN = GetDatasetMNN(cell_mnn, args)
    train_loader_MNN = DataLoader(
        dataset=train_dataset_MNN,
        batch_size=args.batch_size,
        shuffle=True,
        drop_last=True
    )

    train_dataset_Feature = GetDatasetFeaturePair(feature_pair, args)
    train_loader_Feature = DataLoader(
        dataset=train_dataset_Feature,
        batch_size=args.batch_size,
        shuffle=True,
        drop_last=False
    )


    return train_loader_MNN, train_loader_Feature


def test_fea_load(args):

    train_dataset_Other_Feature = GetDatasetModalityFeature(args.out_features[0])
    train_loader_Other_Feature = DataLoader(
        dataset=train_dataset_Other_Feature,
        batch_size=args.batch_size,
        shuffle=True,
        drop_last=False
    )

    train_dataset_RNA_Feature = GetDatasetModalityFeature(args.out_features[1])
    train_loader_RNA_Feature = DataLoader(
        dataset=train_dataset_RNA_Feature,
        batch_size=args.batch_size,
        shuffle=True,
        drop_last=False
    )


    return train_loader_Other_Feature, train_loader_RNA_Feature


def test_cell_load(args):
    test_dataset_Other = GetDatasetModalityCell(args.in_cells[0])
    test_loader_Other = DataLoader(
        dataset=test_dataset_Other,
        batch_size=args.batch_size,
        shuffle=True,
        drop_last=False
    )

    test_dataset_RNA = GetDatasetModalityCell(args.in_cells[1])
    test_loader_RNA = DataLoader(
        dataset=test_dataset_RNA,
        batch_size=args.batch_size,
        shuffle=True,
        drop_last=False
    )


    return test_loader_Other, test_loader_RNA


