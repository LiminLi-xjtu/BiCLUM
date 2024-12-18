
import argparse
import random
import os
import torch
import numpy as np
import sys

from util import load_config, EarlyStopping
from data_preprocess import data_input
from model import Model


batch = sys.argv[1]
config_args = load_config('./config/BMCITE_' + batch)
args = argparse.Namespace(**config_args)



outdir = 'output/' + args.dataset_name + '/'
if not os.path.exists(outdir):
    os.makedirs(outdir)


if __name__ == "__main__":

    use_cuda = torch.cuda.is_available()
    device = torch.device('cuda:0' if use_cuda else 'cpu')
    args.device = device

    seed = args.seed
    np.random.seed(seed)
    random.seed(seed + 1)
    torch.manual_seed(seed + 2)
    torch.cuda.manual_seed(seed + 3)
    torch.backends.cudnn.deterministic = True

    dataset, datahvg, train_loader, test_cell_loader, test_fea_loader, args = data_input(args)

    early_stopping = EarlyStopping(patience=30,
                                   checkpoint_file=os.path.join(outdir, args.GAM_name + '.pt'))
    # Build the model
    model = Model(args.in_features, args.out_features, args.in_cells, # args.device,
                  args.latent_dim, h_depth=2, dropout=0.2, h_dim=256).to(device)

    model.fit(train_loader, dataset, datahvg, early_stopping, device, args)
    torch.save(model.state_dict(), outdir + args.batch + '.pt')

    encoded_cell_list, encoded_fea_list = model.encodeBatch(dataset, datahvg, test_cell_loader, test_fea_loader, args) # , device=device

    eva = dict({"inte_cell": encoded_cell_list, "inte_fea": encoded_fea_list, 'loss_list': model.loss_list})

    path = './result/' + args.dataset_name + '/'

    if not os.path.exists(path):
        os.makedirs(path)

    np.save(os.path.join(path, args.GAM_name + '.npy'), eva)

