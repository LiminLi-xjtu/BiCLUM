
from sklearn.cluster import KMeans

from vis.a_plot_func import *
from eva.evaluation import *
from eva.load_result import *


dataset_name = "PBMC"
dataset_type = 'RNA_ATAC'
GAM_name = 'ArchR'
batch = None
paired = True


methods_gasm = ['BiCLUM']

dataset_dir = "../data/"
result_dir = "../1-comp_methods/results/" + dataset_name
eva_dir = "../1-comp_methods/eva/" + dataset_name
vis_dir = "../1-comp_methods/vis/" + dataset_name

if not os.path.exists(vis_dir):
    os.makedirs(vis_dir)
if not os.path.exists(eva_dir):
    os.makedirs(eva_dir)


adata_gasm, anno_rna_gasm, anno_other_gasm = load_result_gasm(dataset_name, dataset_dir, result_dir, GAM_name, methods_gasm)
obsm_names_gasm = adata_gasm.obsm_keys()

n_clusters = len(np.unique(adata_gasm.obs['cluster']))
cell_types = np.unique(adata_gasm.obs['cell_type'])


###############################Print quanttative metrics#################################################

eva_metrics = evaluate(adata_gasm, anno_rna_gasm, anno_atac_gasm, eva_dir, GAM_name)
# eva_metrics = pd.read_csv(eva_dir + '/eva_' + GAM_name + '.csv', index_col=0, header=0).T

print(eva_metrics)

###################################Plot umap visualization#############################################

omic_colors = ['#92a5d1', '#feb29b']
cell_type_colors_all = ['#f7b3ac', '#e9212c', '#f5c4db', '#d998b6', '#a61e4d', '#bdb5e1', '#a879b8',
                        '#f7e16f','#ff8831', '#958431', '#b0d992', '#82b7a3', '#008d00', '#27483e',
                        '#b4d7f0', '#5591dc', '#919191', '#2017a5', '#ab888d']

cell_type_colors = cell_type_colors_all[:n_clusters]

adata_gasm_ = anndata.AnnData(adata_gasm.obsm['raw_data'])
adata_gasm_.obsm['raw_legend'] = adata_gasm_.X
adata_gasm_.obs = adata_gasm.obs
umap_plot_colors(adata_gasm_, dataset_name, omic_colors, cell_type_colors, vis_dir)

umap_plot_colors(adata_gasm, dataset_name, omic_colors, cell_type_colors, vis_dir)


###################################Plot PAGA graph##############################################

paga_plot_colors(adata_gasm, obsm_names_gasm, cell_type_colors, vis_dir, dataset_name)

###################################Plot confusion matrix#########################################

accuracy = {}
CM = {}

for i, method in enumerate(obsm_names_gasm):
    kmeans = KMeans(n_clusters=len(np.unique(np.hstack((anno_rna_gasm, anno_other_gasm)))), random_state=0).fit(
        adata_gasm.obsm[method])
    clustering = kmeans.labels_
    Y = np.array(adata_gasm.obs['cluster'].array)
    y_preds = get_y_preds(Y, clustering, len(np.unique(Y)))
    scores = clustering_metric(Y, clustering, len(np.unique(Y)))
    accuracy[method] = scores[0]['accuracy']
    CM[method] = scores[1]

confusion_plot(adata_gasm, obsm_names_gasm, CM, accuracy, dataset_name, vis_dir, keys='confusion_plot')

###########################################################################################
