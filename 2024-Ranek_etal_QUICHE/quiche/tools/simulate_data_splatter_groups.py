import numpy as np
import scanpy as sc 
import anndata
import scprep
import pandas as pd

def splatter_sim(batch_cells = 100000,
                n_clusters = 3,
                n_genes = 500,
                bcv_common = 0.1,
                lib_loc = 12,
                group_prob = None,
                random_state = 0):
    """Simulates a single-cell RNA sequencing trajectory using Splatter: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-017-1305-0. 
    ~~~ Uses the scprep wrapper function: https://scprep.readthedocs.io/en/stable/_modules/scprep/run/splatter.html ~~~  
    Parameters
    For more details on the parameters, see: https://scprep.readthedocs.io/en/stable/_modules/scprep/run/splatter.html#SplatSimulate
    ----------
    Returns
    adata: anndata.AnnData
        annotated data object containing simulated single-cell RNA sequecing data (dimensions = cells x features)
    ----------
    """ 
    #set simulation parameters from real single-cell RNA sequencing dataset: https://pubmed.ncbi.nlm.nih.gov/27419872/
    params = {}
    params['group_prob'] = group_prob
    params['bcv_common'] = bcv_common
    params['mean_rate'] = 0.0173
    params['mean_shape'] = 0.54
    if lib_loc is None:
        params['lib_loc'] = 12.6
    else: 
        params['lib_loc'] = lib_loc
    params['lib_scale'] = 0.423
    params['out_prob'] = 0.000342
    params['out_fac_loc'] = 0.1
    params['out_fac_scale'] = 0.4
    params['bcv_df'] = 90.2
    results = scprep.run.SplatSimulate(method = 'groups', 
                                        batch_cells = batch_cells, 
                                        group_prob = params['group_prob'], 
                                        n_genes = n_genes,
                                        de_prob = 0.1,
                                        de_down_prob = 0.5,
                                        de_fac_loc = 0.1,
                                        de_fac_scale = 0.4, 
                                        bcv_common = params['bcv_common'],
                                        mean_rate = params['mean_rate'],
                                        mean_shape = params['mean_shape'],
                                        lib_loc = params['lib_loc'], 
                                        lib_scale = params['lib_scale'], 
                                        out_prob = params['out_prob'], 
                                        out_fac_loc = params['out_fac_loc'], 
                                        out_fac_scale = params['out_fac_scale'], 
                                        bcv_df = params['bcv_df'],
                                        seed = random_state)
    data = pd.DataFrame(results['counts'])
    group = results['group'].copy()
    metadata = pd.DataFrame({'group':group.astype('str')})
    de_genes = pd.concat([pd.DataFrame(results['de_fac_1'], columns = ['group1']),
                            pd.DataFrame(results['de_fac_2'], columns = ['group2']),
                            pd.DataFrame(results['de_fac_3'], columns = ['group3']),
                            pd.DataFrame(results['de_fac_4'], columns = ['group4']),
                            pd.DataFrame(results['de_fac_5'], columns = ['group5'])], axis = 1)            
    gene_index = []
    for i in range(0, len(de_genes.index)):
        if de_genes.loc[i].sum() != n_clusters:
            id = 'DE_group_' + '_'.join(map(str, (np.where(de_genes.loc[i] !=1)[0]))) + '.{}'.format(i)
            gene_index.append(id)
        else:
            gene_index.append(str(i))
    cell_index = pd.Index(['cell_{}'.format(i) for i in range(metadata.shape[0])])
    data.index = cell_index
    data.columns = gene_index
    metadata.index = cell_index
    adata = anndata.AnnData(data)
    adata.obs = metadata
    adata.layers['raw'] = adata.X.copy()
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    return adata

adata = splatter_sim(batch_cells = 525000, n_clusters = 5, n_genes = 500, bcv_common = 0.1, lib_loc = 12,
                    group_prob = [0.2, 0.2, 0.2, 0.2, 0.2], random_state = 0)

adata.write_h5ad(f'adata_simulated_expression_groups_large.h5ad')
