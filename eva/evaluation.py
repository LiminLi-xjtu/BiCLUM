
from eva.load_result import *
from eva.metrics import *

def evaluate(adata, anno_rna, anno_atac, eva_dir, name):

    eva_metrics = {}
    obsm_names = adata.obsm_keys()
    for i in range(len(obsm_names)):
        method = obsm_names[i]
        print(method)
        tmp = all_metrics(adata.obsm[method], np.array(adata.obs['cluster'].array), anno_rna, anno_atac, np.array(adata.obs['omic_id'].array), paired=paired)
        eva_metrics[method] = tmp

    eva_metrics = pd.DataFrame(eva_metrics)
    eva_metrics.index = ['omics_mixing', 'biology_conservation', 'transfer_accuracy', 'foscttm']

    file_path = os.path.join(eva_dir, f'eva_{name}.csv')

    eva_metrics.to_csv(file_path, index=True)

    return eva_metrics

