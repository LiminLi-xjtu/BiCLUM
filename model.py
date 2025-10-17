import torch
import torch.nn.functional as F
from tqdm import tqdm
import numpy as np

from loss import contrastive_loss

class Model(torch.nn.Module):
    def __init__(self, in_features, out_features, in_cells,
                 latent_dim=50, h_depth=2, dropout=0.2, h_dim=256):
        super(Model, self).__init__()

        self.h_depth = h_depth
        self.scale = torch.nn.ParameterList([torch.nn.Parameter(torch.zeros(out_features[k])) for k in range(len(out_features))])
        self.bias = torch.nn.ParameterList([torch.nn.Parameter(torch.zeros(out_features[k])) for k in range(len(out_features))])

        for k in range(len(in_features)):
            setattr(self, f"encoder_c{k}", self.create_encoder(in_features[k], latent_dim, h_dim, dropout))
            setattr(self, f"encoder_f{k}", self.create_encoder(in_cells[k], latent_dim, h_dim, dropout))


    def create_encoder(self, in_dim, latent_dim, h_dim, dropout):
        layers = []
        ptr_dim = in_dim
        for layer in range(self.h_depth):
            layers.append(torch.nn.Linear(ptr_dim, h_dim))
            layers.append(torch.nn.BatchNorm1d(h_dim))
            layers.append(torch.nn.LeakyReLU(negative_slope=0.2))
            layers.append(torch.nn.Dropout(p=dropout))
            ptr_dim = h_dim
        layers.append(torch.nn.Linear(ptr_dim, latent_dim))
        return torch.nn.Sequential(*layers)

    def dec_c(self, u, v, k):
        z = F.softplus(self.scale[k]) * (u @ v.t()) + self.bias[k]
        return z

    def compute_cemb(self, x, k):
        h = getattr(self, f"encoder_c{k}")(x.view(x.size(0), -1))
        return h

    def compute_femb(self, x, k):
        h = getattr(self, f"encoder_f{k}")(x.view(x.size(0), -1))
        return h


    def encodeBatch(self, dataset, datahvg, test_loader, test_fea_loader, args):
        self.eval()

        with torch.no_grad():

            test_loaders = [(test_loader[0], test_fea_loader[0], dataset[0], datahvg[0]),  # other cell
                            (test_loader[1], test_fea_loader[1], dataset[1], datahvg[1])]  # RNA

            n_other, n_rna, = args.in_cells
            m_other, m_rna = args.out_features
            Z_cell = [np.zeros((n_other, args.latent_dim)), np.zeros((n_rna, args.latent_dim))]
            Z_feature = [np.zeros((m_other, args.latent_dim)), np.zeros((m_rna, args.latent_dim))]

            for k, (loader_cell, loader_feature, data, data_hvg) in enumerate(test_loaders):
                for idx_cell in loader_cell:
                    z_tmp = self.compute_cemb(data[idx_cell, :], k)
                    Z_cell[k][idx_cell] = z_tmp.detach().cpu().numpy()
                for idx_feature in loader_feature:
                    feature_data = data_hvg.T[idx_feature, :]
                    z_tmp = self.compute_femb(feature_data, k)
                    Z_feature[k][idx_feature] = z_tmp.detach().cpu().numpy()

            return Z_cell, Z_feature


    def compute_feature_level(self, loader_Feature, datahvg, args, device):

        feature_z = [torch.zeros((args.out_features[0], args.latent_dim)).to(device),
                     torch.zeros((args.out_features[1], args.latent_dim)).to(device)]
        loss_feature_pairs = 0

        for _, (feature_other_idx, feature_rna_idx) in enumerate(loader_Feature):

            fea_idx = [feature_other_idx, feature_rna_idx]
            fea_other_hvg = datahvg[0].T[fea_idx[0], :].to(torch.float32)
            fea_rna_hvg = datahvg[1].T[fea_idx[1], :].to(torch.float32)
            fea = [fea_other_hvg, fea_rna_hvg]

            for k in range(len(fea)):
                feature_z[k][fea_idx[k]] = self.compute_femb(fea[k], k)

            loss_feature_pairs += contrastive_loss(fea[0].size()[0], args.tau_feature, feature_z[0][fea_idx[0]], feature_z[1][fea_idx[1]])


        loss_feature_pairs /= len(loader_Feature)

        return feature_z, loss_feature_pairs


    def compute_cell_level(self, cell_enc, cell_dec, cell_flag, args, feature_z):

        cell_z = []
        cell_r = []

        for k in range(len(cell_enc)):
            cell_z.append(self.compute_cemb(cell_enc[k], k))
            cell_r.append(self.dec_c(cell_z[k], feature_z[k], k))

        loss_rec_cell = 0

        loss_rec_cell += args.gamma_a * F.mse_loss(cell_r[0], cell_dec[0])
        loss_rec_cell += args.gamma_b * F.mse_loss(cell_r[1], cell_dec[1])

        loss_cell_pairs = contrastive_loss(int(torch.sum(cell_flag)), args.tau_cell, cell_z[0][cell_flag, :],
                                           cell_z[1][cell_flag, :])

        return loss_rec_cell, loss_cell_pairs


    def fit(self, train_loader, dataset, datahvg, early_stopping, device, args):

        loader_MNN, loader_Feature = train_loader[0], train_loader[1]
        optimizer = torch.optim.Adam(self.parameters(), lr=args.learning_rate, weight_decay=args.weight_decay)
        loss_list = []
        t_progress = tqdm(range(args.epoch), desc='Training')

        for epoch in t_progress:
            self.train()
            loss_epoch = 0

            for step, (cell_others, cell_rnas, cell_flag) in enumerate(loader_MNN):

                cell_enc = [dataset[0][cell_others,:], dataset[1][cell_rnas, :]]
                cell_dec = [datahvg[0][cell_others, :], datahvg[1][cell_rnas, :]]

                feature_z, loss_feature_pairs = self.compute_feature_level(loader_Feature, datahvg, args, device)
                loss_rec_cell, loss_cell_pairs = self.compute_cell_level(cell_enc, cell_dec, cell_flag, args, feature_z)
                loss = loss_rec_cell + args.alpha * loss_cell_pairs + args.beta * loss_feature_pairs

                optimizer.zero_grad()
                loss.backward()
                optimizer.step()
                loss_epoch += loss.item()

            loss_list.append(loss_epoch / len(loader_MNN))

            early_stopping(loss_epoch, self)

            if early_stopping.early_stop:
                print("Early stopping")
                self.loss_list = loss_list
                break

        self.loss_list = loss_list



