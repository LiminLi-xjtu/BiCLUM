
import numpy as np
import torch
import yaml
import json


def load_config(config_name, verbose=True):
    if '.yaml' not in config_name:
        config_name += '.yaml'
    with open(config_name, 'r') as f:
        f_str = f.read()
        dic = yaml.safe_load(f_str)
        if verbose:
            js = json.dumps(dic, sort_keys=True, indent=4, separators=(',', ':'))
            print(js)
        return dic



class EarlyStopping:

    def __init__(self, patience=10, verbose=False, checkpoint_file=''):

        self.patience = patience
        self.verbose = verbose
        self.counter = 0
        self.best_score = None
        self.early_stop = False
        self.loss_min = np.Inf
        self.checkpoint_file = checkpoint_file

    def __call__(self, loss, model):
        if np.isnan(loss):
            self.early_stop = True
        score = -loss

        if self.best_score is None:
            self.best_score = score
            self.save_checkpoint(loss, model)
        elif score <= self.best_score:
            self.counter += 1
            if self.verbose:
                print(f'EarlyStopping counter: {self.counter} out of {self.patience}')
            if self.counter >= self.patience:
                self.early_stop = True
        else:
            self.best_score = score
            self.save_checkpoint(loss, model)
            self.counter = 0

    def save_checkpoint(self, loss, model):
        '''
        Saves model when loss decrease.
        '''
        if self.verbose:
            print(f'Loss decreased ({self.loss_min:.6f} --> {loss:.6f}).  Saving model ...')
        torch.save(model.state_dict(), self.checkpoint_file)
        self.loss_min = loss


