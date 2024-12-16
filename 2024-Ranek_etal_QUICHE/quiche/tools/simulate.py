import pandas as pd
import numpy as np
import os
import anndata
from sketchKH import *
import ark.settings as settings
from scipy.spatial import cKDTree
import quiche as qu
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Rectangle
import random
from scipy.spatial.distance import cdist
import scanpy as sc
import shutil

def select_random_point(annotations_by_mask = None,
                        fov_key = 'fov',
                        fov = None,
                        labels_key = 'mask_name',
                        mask = 'cancer_core',
                        radius = 200,
                        p = 2,
                        num_to_change = 80,
                        condition_id = 'immune1',
                        group = 'core_immune1',
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

        spatial_kdTree_truth = cKDTree(fov_subset.loc[:, ['centroid-1', 'centroid-0']])
        nn_truth = spatial_kdTree_truth.query_ball_point(fov_subset[np.isin(fov_subset.label, random_label)].loc[:, ['centroid-1', 'centroid-0']], r=radius, p = p, workers = n_jobs) #no self interactions
        ground_truth_mask = fov_subset[fov_subset.isin(fov_subset.iloc[nn_truth[0], :])].dropna()['ground_labels']

        annotations_by_mask.loc[ground_truth_mask.index, 'ground_labels'] = 1
        annotations_by_mask.loc[ground_truth_mask.index, 'DA_group'] = group    
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

    for fov in fov_list:
        annotations_by_mask_merged = pd.merge(annotations_by_mask, cell_table.loc[:, [fov_key, 'label', 'centroid-1', 'centroid-0']], on = [fov_key, 'label'])
        annotations_by_mask_merged[labels_key].replace(['stroma_core', 'stroma_border'], ['stroma', 'stroma'], inplace = True)
        annotations_by_mask_merged['ground_labels'] = 0
        annotations_by_mask_merged['DA_group'] = 'random'

        try:
            mask_label = mask.split('_')[1]
        except:
            mask_label = 'none'

        if fov in selected_fovs:
            annotations_by_mask_merged = select_random_point(annotations_by_mask = annotations_by_mask_merged, fov_key = fov_key, fov = fov, labels_key = labels_key, mask = mask,
                                                            radius = radius, p = p, num_to_change = num_to_change, condition_id = condition_id, n_jobs = n_jobs, group = f'{mask_label}_immune1')
        else:
            annotations_by_mask_merged = select_random_point(annotations_by_mask = annotations_by_mask_merged, fov_key = fov_key, fov = fov, labels_key = labels_key, mask = None,
                                                            radius = radius, p = p, num_to_change = num_to_change, condition_id = condition_id, n_jobs = n_jobs, group = f'{mask_label}_immune1')

        count_df = annotations_by_mask_merged[np.isin(annotations_by_mask_merged.fov, fov)].groupby([fov_key, labels_key]).count()['label'].unstack()
        annotations_by_mask_merged = stratify_sampling(annotations_by_mask_merged, count_df, ratio_df, [fov], labels_key = labels_key)
                    
        if annotations_by_mask_merged is not None:    
            
            condition_df_ = annotations_by_mask_merged[np.isin(annotations_by_mask_merged.fov, fov)].copy()
            condition_df = pd.concat([condition_df, condition_df_], axis = 0)

    condition_df_merged = pd.merge(cell_table, condition_df, on = ['fov', 'label','centroid-1', 'centroid-0'])
    adata = anndata.AnnData(condition_df_merged)
    adata.obs = condition_df_merged.loc[:, ['fov', 'label', labels_key, 'ground_labels', 'DA_group']].copy()
    adata.obsm['spatial'] = condition_df_merged.loc[:, ['centroid-1', 'centroid-0']].values
    adata.obs['condition'] = mask
    adata.obs['Patient_ID'] = adata.obs['fov'] + adata.obs['condition'] 
    adata.obs['Patient_ID'] = pd.Series(adata.obs['Patient_ID']).astype('category')
    adata.obs[labels_key] = pd.Series(adata.obs[labels_key]).astype('category')
    adata.obs_names = [f'c_{i}_{mask}_{adata.obs.Patient_ID[i]}' for i in range(0, len(adata.obs_names))]

    return adata

def simulate_structured_data(annotations_by_mask,
                cell_table,
                adata_expression,
                n_cond = 10,
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
                                fov_list = np.random.choice(fov_list, n_cond, replace = False),
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
                                fov_list = np.random.choice(fov_list, n_cond, replace = False),
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

def spatial_regions_unstructured(num_grids_x = 5, num_grids_y = 5, n_regions = 3, scale = 1000, n_niches = 1000, da_vec = ['A', 'C', 'E'], seed = 0,
                     colors_dict = {'A': '#B46CDA','B': '#78CE8B', 'C': '#FF8595', 'D': '#1885F2', 'E': '#D78F09'}, hex = '#e41a1c', sample_size = None, 
                     fc_dict = {'A': 1,'B': 3, 'C': 5, 'D': 1.5, 'E': 2}, show_grid = False, save_directory = 'figures', filename_save = None):

    np.random.seed(seed)
    data = {'x': np.random.rand(n_niches)*scale,
        'y': np.random.rand(n_niches)*scale}
    
    df = pd.DataFrame(data)

    grid_size_x = scale / num_grids_x
    grid_size_y = scale / num_grids_y

    if show_grid == True:
        for i in range(1, num_grids_x):
            plt.axvline(i * grid_size_x, color='gray', linestyle='--', linewidth=0.5)
        for i in range(1, num_grids_y):
            plt.axhline(i * grid_size_y, color='gray', linestyle='--', linewidth=0.5)

    df['group'] = np.nan
    df['group'] = df['group'].astype('object')
    
    df['DA_group'] = 'random'
    df['DA_group'] = df['DA_group'].astype('object')

    df['DA_group_center'] = 'random'
    df['DA_group_center'] = df['DA_group_center'].astype('object')

    available_grids = set([(x, y) for x in range(1, num_grids_x + 1) for y in range(1, num_grids_y + 1)])
    grids = []
    while len(grids) < n_regions:
        grid_ = random.sample(available_grids, 1)[0]
        if all((abs(grid_[0] - x) > 1 or abs(grid_[1] - y) > 1) for (x, y) in grids):
            grids.append(grid_)

    selected_grid_list = []
    for grid in grids:
        selected_grid_x = grid[0]
        selected_grid_y = grid[1]
        selected_grid_list.append((grid[0], grid[1]))
        selected_locations = df[
            (df['x'] >= (selected_grid_x - 1) * grid_size_x) & (df['x'] < selected_grid_x * grid_size_x) &
            (df['y'] >= (selected_grid_y - 1) * grid_size_y) & (df['y'] < selected_grid_y * grid_size_y)
        ].index

        df.loc[selected_locations, 'group'] = np.random.choice(da_vec, size=len(selected_locations))
        df.loc[selected_locations, 'DA_group'] = '_'.join(da_vec)

        centroid = ((selected_grid_x - 1 + selected_grid_x) * 0.5 * grid_size_x,
            (selected_grid_y - 1 + selected_grid_y) * 0.5 * grid_size_y)
        distances = cdist(df.loc[selected_locations, ['x', 'y']], [centroid])
        closest_indices = np.argsort(distances.flatten())[0]
        df.loc[selected_locations[closest_indices], 'DA_group_center'] = '_'.join(da_vec)

    labels = df.loc[:, 'group'].value_counts()
    remaining_sample_size = sample_size.copy()
    for label, count in labels.items():
        remaining_sample_size[label] -= count

    relabel_arr = np.concatenate([[label] * count for label, count in remaining_sample_size.items() if count > 0])
    df.loc[df['group'].isnull(), 'group'] = np.random.choice(relabel_arr, size=df['group'].isnull().sum(), replace=False)

    df['foldchange'] = pd.Series(df['group']).map(fc_dict)
    df.index = [f'Loc{i}' for i in range(1, df.shape[0]+1)]

    param_dict = {'seed': seed,
                  'num_grids_x': num_grids_x,
                  'num_grids_y': num_grids_y,
                  'n_regions': n_regions,
                  'da_vec': da_vec,
                  'n_niches': n_niches,
                  'selected_grids': selected_grid_list}
    
    return df, param_dict

def simulate_unstructured(n_patients_condA = 10,
                          n_patients_condB = 10,
                        num_grids_x = 5,
                        num_grids_y = 5,
                        n_niches_A = 1000,
                        n_niches_B = 1000,
                        scale = 1,
                        n_regionsA = 1,
                        n_regionsB = [1],
                        ratio = 0.2,
                        sample_size_A = None,
                        sample_size_B = None, 
                        da_vec_A = ['A', 'C', 'E'],
                        da_vec_B = [['B', 'D']],
                        hex_A = '#e41a1c',
                        hex_B = '#377eb8',
                        random_state_list_A = [58, 322, 1426, 65, 651, 417, 2788, 576, 213, 1828],
                        random_state_list_B = [51, 1939, 2700, 1831, 804, 2633, 2777, 2053, 948, 420], 
                        fig_id = 'ACE_BD_region1_ratio0.2',
                        save_directory = 'data/simulated'):
    
    df_condA = pd.DataFrame()
    param_condA = {'seed':[], 'num_grids_x': [], 'num_grids_y': [], 'n_regions': [], 'da_vec': [], 'n_niches': [], 'selected_grids': []}
    run = 0
    for i in range(0, n_patients_condA):
        if run < int(n_patients_condA*ratio):
            df, param_dict = spatial_regions_unstructured(scale = scale, num_grids_x = num_grids_x, hex = hex_A, num_grids_y = num_grids_y, n_regions = n_regionsA, n_niches = n_niches_A, da_vec = da_vec_A, seed = random_state_list_A[i], save_directory = save_directory, filename_save=None, sample_size=sample_size_A)
        else:
            df, param_dict = spatial_regions_unstructured(scale = scale, num_grids_x = num_grids_x, hex = hex_A, num_grids_y = num_grids_y, n_regions = 0, n_niches = n_niches_A, da_vec = da_vec_A, seed = random_state_list_A[i], save_directory = save_directory, filename_save=None, sample_size=sample_size_A)
        df.to_csv(os.path.join(os.path.join(save_directory, f'spatial_condA{i}_{fig_id}.csv')))
        df['Patient_ID'] = i
        df_condA = pd.concat([df_condA, df], axis = 0)
        param_condA['seed'].append(param_dict['seed'])
        param_condA['num_grids_x'].append(param_dict['num_grids_x'])
        param_condA['num_grids_y'].append(param_dict['num_grids_y'])
        param_condA['n_regions'].append(param_dict['n_regions'])
        param_condA['da_vec'].append(param_dict['da_vec'])
        param_condA['n_niches'].append(param_dict['n_niches'])
        param_condA['selected_grids'].append(param_dict['selected_grids'])
        run +=1
        
    df_condB = pd.DataFrame()
    param_condB = {'seed':[], 'num_grids_x': [], 'num_grids_y': [], 'n_regions': [], 'da_vec': [], 'n_niches': [], 'selected_grids': []}
    run = 0
    for i in range(0, n_patients_condB):
        if run < int(n_patients_condB*ratio):
            df, param_dict = spatial_regions_unstructured(scale = scale, num_grids_x = num_grids_x, hex = hex_B, num_grids_y = num_grids_y, n_regions = n_regionsB, n_niches = n_niches_B, da_vec = da_vec_B, seed = random_state_list_B[i], save_directory = save_directory, filename_save=None, sample_size=sample_size_B)
        else:
            df, param_dict = spatial_regions_unstructured(scale = scale, num_grids_x = num_grids_x, hex = hex_B, num_grids_y = num_grids_y, n_regions = 0, n_niches = n_niches_B, da_vec = da_vec_B, seed = random_state_list_B[i], save_directory = save_directory, filename_save=None, sample_size=sample_size_B)
        df.to_csv(os.path.join(os.path.join(save_directory, f'spatial_condB{i}_{fig_id}.csv')))
        df['Patient_ID'] = i
        df_condB = pd.concat([df_condB, df], axis = 0)
        param_condB['seed'].append(param_dict['seed'])
        param_condB['num_grids_x'].append(param_dict['num_grids_x'])
        param_condB['num_grids_y'].append(param_dict['num_grids_y'])
        param_condB['n_regions'].append(param_dict['n_regions'])
        param_condB['da_vec'].append(param_dict['da_vec'])
        param_condB['n_niches'].append(param_dict['n_niches'])
        param_condB['selected_grids'].append(param_dict['selected_grids'])
        run+=1

    ##aggregates counts with metadata 
    adata_simulated = []
    adata_run = anndata.read_h5ad(os.path.join('data', 'simulated', 'adata_simulated_expression_groups_large.h5ad'))
    
    for cond in ['A', 'B']:
        if cond == 'A':
            sample_size = sample_size_A.copy()
            n_patients_cond = n_patients_condA
        else:
            sample_size = sample_size_B.copy()
            n_patients_cond = n_patients_condB
        for i in range(0, n_patients_cond): #can change this if we want class imbalance        
            expression_df = pd.DataFrame(adata_run.X, index = adata_run.obs_names, columns = adata_run.var_names)
            location = pd.read_csv(os.path.join(os.path.join(save_directory, f'spatial_cond{cond}{i}_{fig_id}.csv')), index_col=0)
            sample = list(adata_run.obs['group'][adata_run.obs['group'] == 'Group1'].sample(sample_size['A']).index)
            A = expression_df.loc[sample, :]
            #adata_run = adata_run[~np.isin(adata_run.obs_names, sample)]
            A.index = location[location['group'] == 'A'].index
            sample = list(adata_run.obs['group'][adata_run.obs['group'] == 'Group2'].sample(sample_size['B']).index)
            B = expression_df.loc[sample, :]
            #adata_run = adata_run[~np.isin(adata_run.obs_names, sample)]
            B.index = location[location['group'] == 'B'].index
            sample = list(adata_run.obs['group'][adata_run.obs['group'] == 'Group3'].sample(sample_size['C']).index)
            C = expression_df.loc[sample, :]
            #adata_run = adata_run[~np.isin(adata_run.obs_names, sample)]
            C.index = location[location['group'] == 'C'].index
            sample = list(adata_run.obs['group'][adata_run.obs['group'] == 'Group4'].sample(sample_size['D']).index)
            D = expression_df.loc[sample, :]
            #adata_run = adata_run[~np.isin(adata_run.obs_names, sample)]
            D.index = location[location['group'] == 'D'].index
            sample = list(adata_run.obs['group'][adata_run.obs['group'] == 'Group5'].sample(sample_size['E']).index)
            E = expression_df.loc[sample, :]
            #adata_run = adata_run[~np.isin(adata_run.obs_names, sample)]
            E.index = location[location['group'] == 'E'].index

            expression_df = pd.concat([A, B, C, D, E], axis = 0)
            expression_df = expression_df.loc[location.index]
            adata = anndata.AnnData(expression_df)
            adata.obsm['spatial'] = location.loc[:, ['x', 'y']].values
            adata.obs['cell_cluster'] = location.loc[:, 'group']
            adata.obs['cell_cluster'] = adata.obs['cell_cluster'].astype('category')
            adata.obs['label'] = [i for i in range(0, len(expression_df))]
            adata.obs['DA_group'] = location.loc[:, 'DA_group']
            adata.obs['DA_group_center'] = location.loc[:, 'DA_group_center']
            adata.obs['condition'] = cond
            adata.obs['Patient_ID'] = f'{cond}{i}'
            adata.obs['Patient_ID'] = adata.obs['Patient_ID'].astype('category')
            adata.obs_names = [f'{cond}{i}_{j}' for j in adata.obs_names] #make unique obs names
            adata.obs['ground_labels'] = 0
            adata.obs['ground_labels'][np.where((np.isin(adata.obs['DA_group_center'], ['_'.join(da_vec_A)])) | (np.isin(adata.obs['DA_group_center'], ['_'.join(da_vec_B)]) ))[0]] = 1
            adata_simulated.append(adata)
    adata_simulated = anndata.concat(adata_simulated)
    return adata_simulated