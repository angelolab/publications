import pandas as pd
import numpy as np
import os
import anndata
from sketchKH import *
from rpy2.robjects import pandas2ri
from ark.utils.plot_utils import cohort_cluster_plot
import ark.settings as settings
from scipy.spatial import cKDTree

pandas2ri.activate()

def select_random_point(annotations_by_mask = None,
                        fov_key = 'fov',
                        fov = None,
                        labels_key = 'mask_name',
                        mask = 'cancer_core',
                        radius = 200,
                        p = 2,
                        num_to_change = 80,
                        condition_id = 'immune1',
                        n_jobs = 8):
    
    fov_subset = annotations_by_mask[annotations_by_mask[fov_key] == fov].copy()
    if mask is None: #if enrichment location is not specified, then randomly place
        annotation_subset = fov_subset.copy()
        random_label = np.random.choice(list(annotation_subset.label), num_to_change, replace = False)
        bool_mask = annotation_subset[np.isin(annotation_subset.label, random_label)][labels_key]
        annotations_by_mask.loc[bool_mask.index, labels_key] = condition_id
    else:
        annotation_subset = fov_subset[annotations_by_mask[labels_key] == mask].copy()
        random_label = np.random.choice(list(annotation_subset.label))
        spatial_kdTree = cKDTree(annotation_subset.loc[:, ['centroid-1', 'centroid-0']])
        nn = spatial_kdTree.query_ball_point(annotation_subset[np.isin(annotation_subset.label, random_label)].loc[:, ['centroid-1', 'centroid-0']], r=radius, p = p, workers = n_jobs) #no self interactions
        nn = np.array(nn[0])
        indices_to_change = nn[np.random.choice(len(nn), size=num_to_change, replace=False)]
        updated_mask = annotations_by_mask[annotations_by_mask.isin(annotation_subset.iloc[indices_to_change, :])].dropna()[labels_key]
        annotations_by_mask.loc[updated_mask.index, labels_key] = condition_id

        #get all cells
        spatial_kdTree_truth = cKDTree(fov_subset.loc[:, ['centroid-1', 'centroid-0']])
        nn_truth = spatial_kdTree_truth.query_ball_point(fov_subset[np.isin(fov_subset.label, random_label)].loc[:, ['centroid-1', 'centroid-0']], r=radius, p = p, workers = n_jobs) #no self interactions
        ground_truth_mask = fov_subset[fov_subset.isin(fov_subset.iloc[nn_truth[0], :])].dropna()['ground_truth_DA']
        annotations_by_mask.loc[ground_truth_mask.index, 'ground_truth_DA'] = 1
    
    return annotations_by_mask


def stratify_sampling(annotations_by_mask_merged, count_df, ratio_df, fov_list, labels_key = 'mask_name'):
    # downsample cells per cell type by randomly removing cells, they get nan for mask label
    annotations_by_mask_subset = pd.DataFrame()
    for fov in fov_list:
        annotations_by_mask_fov = annotations_by_mask_merged.iloc[np.isin(annotations_by_mask_merged.fov, fov), :]
        cell_type_counts_fov = count_df.loc[fov].copy()
        cells2remove = cell_type_counts_fov - ratio_df
        for mask in cells2remove.index:
            annotations_by_mask_fov_mask = annotations_by_mask_fov[annotations_by_mask_fov.mask_name == mask].copy()
            cell_ids = np.random.choice(annotations_by_mask_fov_mask.label, size = cells2remove.loc[mask], replace = False)
            bool_idx = np.isin(annotations_by_mask_fov_mask.label, cell_ids)
            annotations_by_mask_fov_mask.loc[bool_idx, labels_key] = np.nan
            annotations_by_mask_subset = pd.concat([annotations_by_mask_subset, annotations_by_mask_fov_mask], axis = 0)
    return annotations_by_mask_subset

def simulate_condition(annotations_by_mask = None,
                       cell_table = None,
                       fov_key = 'fov',
                       fov_list = None,
                       labels_key = 'mask_name',
                       mask = 'cancer_core',
                       radius = 200,
                       p = 2,
                       condition_id = 'immune1',
                       num_to_change = 80,
                       seg_dir = r'/Volumes/Shared/Noah Greenwald/TNBC_Cohorts/SPAIN/segmentation/samples/deepcell_output',
                       compartment_colormap = None,
                       prevalence = 0.2,
                       ratio_df = None,
                       save_dir = None,
                       n_jobs = 8):
    
    condition_df = pd.DataFrame()
    selected_fovs = np.random.choice(fov_list, size = int(len(fov_list)*prevalence), replace=False)
    
    print(selected_fovs)
   
    for fov in fov_list:
        annotations_by_mask_merged = pd.merge(annotations_by_mask, cell_table.loc[:, [fov_key, 'label', 'centroid-1', 'centroid-0']], on = [fov_key, 'label'])
        annotations_by_mask_merged[labels_key].replace(['stroma_core', 'stroma_border'], ['stroma', 'stroma'], inplace = True)
        annotations_by_mask_merged['ground_truth_DA'] = 0

        if fov in selected_fovs:
            annotations_by_mask_merged = select_random_point(annotations_by_mask = annotations_by_mask_merged, fov_key = fov_key, fov = fov, labels_key = labels_key, mask = mask,
                                                            radius = radius, p = p, num_to_change = num_to_change, condition_id = condition_id, n_jobs = n_jobs)
        else:
            annotations_by_mask_merged = select_random_point(annotations_by_mask = annotations_by_mask_merged, fov_key = fov_key, fov = fov, labels_key = labels_key, mask = None,
                                                            radius = radius, p = p, num_to_change = num_to_change, condition_id = condition_id, n_jobs = n_jobs)

        count_df = annotations_by_mask_merged[np.isin(annotations_by_mask_merged.fov, fov)].groupby([fov_key, labels_key]).count()['label'].unstack()
        annotations_by_mask_merged = stratify_sampling(annotations_by_mask_merged, count_df, ratio_df, [fov], labels_key = labels_key)
                    
        if annotations_by_mask_merged is not None:    
            cluster_mask = cohort_cluster_plot(fovs=[fov],
                                                seg_dir=seg_dir,
                                                save_dir=save_dir,
                                                cell_data=annotations_by_mask_merged,
                                                erode=True,
                                                fov_col=settings.FOV_ID,
                                                label_col=settings.CELL_LABEL,
                                                cluster_col='mask_name',
                                                seg_suffix="_whole_cell.tiff",
                                                cmap=compartment_colormap,
                                                display_fig=False)
            
            condition_df_ = annotations_by_mask_merged[np.isin(annotations_by_mask_merged.fov, fov)].copy()
            condition_df = pd.concat([condition_df, condition_df_], axis = 0)

    condition_df_merged = pd.merge(cell_table, condition_df, on = ['fov', 'label','centroid-1', 'centroid-0'])
    adata = anndata.AnnData(condition_df_merged)
    adata.obs = condition_df_merged.loc[:, ['fov', 'label', labels_key, 'ground_truth_DA']].copy()
    adata.obsm['spatial'] = condition_df_merged.loc[:, ['centroid-1', 'centroid-0']].values
    adata.obs['condition'] = mask
    adata.obs['Patient_ID'] = adata.obs['fov'] + adata.obs['condition'] 
    adata.obs['Patient_ID'] = pd.Categorical(adata.obs['Patient_ID'])
    adata.obs_names = [f'c_{i}_{mask}_{adata.obs.Patient_ID[i]}' for i in range(0, len(adata.obs_names))]

    return adata

def simulate_data(annotations_by_mask,
                cell_table,
                adata_expression,
                fov_key = 'fov',
                fov_list = None,
                labels_key = 'mask_name',
                cond1 = 'cancer_core',
                cond2 = 'cancer_border',
                radius = 200,
                p = 2,
                condition_id = 'immune1',
                num_to_change = 80,
                compartment_colormap = None,
                prevalence = 0.2,
                ratio_df = None, 
                cell_types = ['cancer_core', 'cancer_border', 'stroma', 'immune1'],
                sim_cell_types = ['Group1', 'Group1', 'Group2', 'Group3'],
                group_key = 'group',
                spatial_key = 'spatial',
                n_jobs = 8):
    
    save_dir_cond1 = os.path.join(cond1, 'prev'+str(prevalence), 'r'+str(radius), 'n'+str(num_to_change)) 
    save_dir_cond2 = os.path.join(cond2, 'prev'+str(prevalence), 'r'+str(radius), 'n'+str(num_to_change))

    condition1 = simulate_condition(annotations_by_mask = annotations_by_mask,
                                cell_table = cell_table,
                                fov_key = fov_key,
                                fov_list = np.random.choice(fov_list, 10, replace = False),
                                labels_key = labels_key,
                                mask = cond1,
                                radius = radius,
                                p = p,
                                num_to_change=num_to_change,
                                prevalence = prevalence, 
                                condition_id = condition_id,
                                compartment_colormap = compartment_colormap,
                                ratio_df = ratio_df,
                                save_dir = save_dir_cond1,
                                n_jobs = n_jobs)

    condition2 = simulate_condition(annotations_by_mask = annotations_by_mask,
                                cell_table = cell_table,
                                fov_key = fov_key,
                                fov_list = np.random.choice(fov_list, 10, replace = False),
                                labels_key = labels_key,
                                mask = cond2,
                                radius = radius,
                                p = p,
                                num_to_change=num_to_change,
                                prevalence = prevalence, 
                                condition_id = condition_id,
                                compartment_colormap = compartment_colormap,
                                ratio_df = ratio_df,
                                save_dir = save_dir_cond2,
                                n_jobs = n_jobs)
    
    adata_spatial = anndata.concat([condition1, condition2])
    adata_spatial = adata_spatial[~adata_spatial.obs[labels_key].isna()]

    adata = assign_simulated_expression(adata_spatial, adata_expression, cell_types = cell_types, sim_cell_types = sim_cell_types,
                                              labels_key = labels_key, group_key = group_key, spatial_key = spatial_key)
    
    return adata

def assign_simulated_expression(adata_spatial,
                                adata_expression, 
                                cell_types = ['cancer_core', 'cancer_border', 'stroma', 'immune1'],
                                sim_cell_types = ['Group1', 'Group1', 'Group2', 'Group3'],
                                labels_key = 'mask_name',
                                group_key = 'group',
                                spatial_key = 'spatial'):

    relabeled_expression = np.full(shape = (adata_spatial.X.shape[0], adata_expression.X.shape[1]), fill_value = np.nan)
    for i in range(0, len(cell_types)):
        ncells = adata_spatial[adata_spatial.obs[labels_key] == cell_types[i]].shape[0]
        idx_sim_cells = np.random.choice(adata_expression[adata_expression.obs[group_key] == sim_cell_types[i]].obs_names, size = ncells, replace = False)
        relabeled_expression[adata_spatial.obs[labels_key] == cell_types[i]] = adata_expression[np.isin(adata_expression.obs_names, idx_sim_cells)].X.copy()
    adata_relabeled = anndata.AnnData(pd.DataFrame(relabeled_expression, index = adata_spatial.obs_names, columns = adata_expression.var_names))
    adata_relabeled.obs = adata_spatial.obs
    adata_relabeled.obsm[spatial_key] = adata_spatial.obsm[spatial_key]
    return adata_relabeled

# import numpy as np
# import scanpy as sc 
# import anndata
# import delve_benchmark
# import gc
# import logging
# import rpy2.robjects as robjects
# import gc 
# from sklearn.model_selection import train_test_split
# import scprep
# import pandas as pd

# def linear_mask_(metadata):
#     # makes step variable monotonically increeasing for the linear trajectory
#     metadata_ = metadata.copy()
#     mask_root = metadata_['group'] == 'Path1'
#     metadata_.loc[mask_root, 'step'] = 100 - metadata_.loc[mask_root, 'step']
#     for i in [2,3,4,5]:
#         mask = metadata_['group'] == 'Path'+str(i)
#         metadata_.loc[mask, 'step'] = 100*(i-1) + metadata_.loc[mask, 'step']
#     return metadata_

# def _sum_to_one(x):
#     x = x / np.sum(x)
#     x = x.round(3)
#     if np.sum(x) != 1:
#         x[0] += 1 - np.sum(x)
#     x = x.round(3)
#     return x

# def splatter_sim(cells_per_path = 200,
#                 n_paths = 5,
#                 n_genes = 500,
#                 bcv_common = 0.1,
#                 lib_loc = 12,
#                 path_from = None,
#                 path_skew = None,
#                 path_type = None,
#                 group_prob = None,
#                 random_state = 0):
#     """Simulates a single-cell RNA sequencing trajectory using Splatter: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-017-1305-0. 
#     ~~~ Uses the scprep wrapper function: https://scprep.readthedocs.io/en/stable/_modules/scprep/run/splatter.html ~~~  
#     Parameters
#     For more details on the parameters, see: https://scprep.readthedocs.io/en/stable/_modules/scprep/run/splatter.html#SplatSimulate
#     ----------
#     Returns
#     adata: anndata.AnnData
#         annotated data object containing simulated single-cell RNA sequecing data (dimensions = cells x features)
#     ----------
#     """ 
#     #set simulation parameters from real single-cell RNA sequencing dataset: https://pubmed.ncbi.nlm.nih.gov/27419872/
#     params = {}
#     params['group_prob'] = group_prob
#     params['bcv_common'] = bcv_common
#     params['path_from'] = path_from
#     params['path_skew'] = path_skew
#     params['mean_rate'] = 0.0173
#     params['mean_shape'] = 0.54
#     if lib_loc is None:
#         params['lib_loc'] = 12.6
#     else: 
#         params['lib_loc'] = lib_loc
#     params['lib_scale'] = 0.423
#     params['out_prob'] = 0.000342
#     params['out_fac_loc'] = 0.1
#     params['out_fac_scale'] = 0.4
#     params['bcv_df'] = 90.2
#     results = scprep.run.SplatSimulate(method = 'paths', 
#                                         batch_cells = [cells_per_path * n_paths], 
#                                         group_prob = params['group_prob'], 
#                                         n_genes = n_genes,
#                                         de_prob = 0.1,
#                                         de_down_prob = 0.5,
#                                         de_fac_loc = 0.1,
#                                         de_fac_scale = 0.4, 
#                                         bcv_common = params['bcv_common'],
#                                         dropout_type = 'none',
#                                         path_from = params['path_from'],
#                                         path_skew = params['path_skew'],
#                                         mean_rate = params['mean_rate'],
#                                         mean_shape = params['mean_shape'],
#                                         lib_loc = params['lib_loc'], 
#                                         lib_scale = params['lib_scale'], 
#                                         out_prob = params['out_prob'], 
#                                         out_fac_loc = params['out_fac_loc'], 
#                                         out_fac_scale = params['out_fac_scale'], 
#                                         bcv_df = params['bcv_df'],
#                                         seed = random_state)
#     data = pd.DataFrame(results['counts'])
#     group = results['group'].copy()
#     metadata = pd.DataFrame({'group':group.astype('str'), 'step':results['step'].astype(int)})
#     if path_type == 'linear':
#         metadata = linear_mask_(metadata)
#     elif path_type == 'branch':
#         metadata = branch_mask_(metadata)
#     de_genes = pd.concat([pd.DataFrame(results['de_fac_1'], columns = ['path1']),
#                             pd.DataFrame(results['de_fac_2'], columns = ['path2']),
#                             pd.DataFrame(results['de_fac_3'], columns = ['path3']),
#                             pd.DataFrame(results['de_fac_4'], columns = ['path4']),
#                             pd.DataFrame(results['de_fac_5'], columns = ['path5'])], axis = 1)            
#     gene_index = []
#     for i in range(0, len(de_genes.index)):
#         if de_genes.loc[i].sum() != n_paths:
#             id = 'DE_group_' + '_'.join(map(str, (np.where(de_genes.loc[i] !=1)[0]))) + '.{}'.format(i)
#             gene_index.append(id)
#         else:
#             gene_index.append(str(i))
#     cell_index = pd.Index(['cell_{}'.format(i) for i in range(metadata.shape[0])])
#     data.index = cell_index
#     data.columns = gene_index
#     metadata.index = cell_index
#     adata = anndata.AnnData(data)
#     adata.obs = metadata
#     adata.layers['raw'] = adata.X.copy()
#     sc.pp.normalize_total(adata)
#     sc.pp.log1p(adata)
#     return adata

# n_paths = 5
# group_prob = np.random.dirichlet(np.ones(n_paths) * 1.).round(3)
# group_prob = _sum_to_one(group_prob)
# path_skew = np.random.beta(10., 10., n_paths)
# adata = splatter_sim(cells_per_path = 3000, n_paths = n_paths, n_genes = 500, bcv_common = 0.1, lib_loc = 12,
#                     path_from = [0, 1, 2, 3, 4], path_skew = path_skew,
#                     group_prob = [1,1,1,1,1], path_type = 'linear', random_state = 0)
# adata.write_h5ad(f'adata_simulated_expression.h5ad')
