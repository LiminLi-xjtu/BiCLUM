
from sklearn.cluster import KMeans
import sys
import argparse

from util import load_config
from vis.plot_func import *
from eva.load_result import *

dataset_name = sys.argv[1]
config_args = load_config('./config/' + dataset_name)
args = argparse.Namespace(**config_args)


methods = ['BiCLUM']

dataset_dir = "../data/"
result_dir = "./results/" + args.dataset_name
eva_dir = "./eva/" + args.dataset_name
vis_dir = "./vis/" + args.dataset_name

if not os.path.exists(vis_dir):
    os.makedirs(vis_dir)
if not os.path.exists(eva_dir):
    os.makedirs(eva_dir)


adata, anno_rna, anno_other = load_cite(args.dataset_name, dataset_dir, result_dir, methods, args.batch)
obsm_names = adata.obsm_keys()

n_clusters = len(np.unique(adata.obs['cluster']))
cell_types = np.unique(adata.obs['cell_type'])



###################################Plot umap visualization#############################################

omic_colors = ['#92a5d1', '#feb29b']
cell_type_colors_all = ['#f7b3ac', '#e9212c', '#f5c4db', '#d998b6', '#a61e4d', '#bdb5e1', '#a879b8',
                        '#f7e16f','#ff8831', '#958431', '#b0d992', '#82b7a3', '#008d00', '#27483e',
                        '#b4d7f0', '#5591dc', '#919191', '#2017a5', '#ab888d']

cell_type_colors = cell_type_colors_all[:n_clusters]

umap_plot_colors(adata, args.dataset_name, omic_colors, cell_type_colors, vis_dir)


###################################Plot PAGA graph##############################################

paga_plot_colors(adata, obsm_names, cell_type_colors, vis_dir, args.dataset_name)

###########################################################################################
