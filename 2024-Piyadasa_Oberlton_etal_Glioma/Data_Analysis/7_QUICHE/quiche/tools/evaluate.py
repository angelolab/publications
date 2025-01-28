import pandas as pd
import numpy as np
import os
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
from ark.utils.plot_utils import cohort_cluster_plot
import ark.settings as settings
import quiche as qu
from sketchKH import *
from scipy.spatial import cKDTree
import anndata
import random
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
from muon import MuData
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
import pandas as pd
import numpy as np
import squidpy as sq
from tqdm import tqdm
import os
## pairwise enrichment
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns
import pertpy as pt
import scanpy as sc
from statsmodels.stats.multitest import multipletests
import statsmodels.api as sm
from rpy2.robjects import pandas2ri
pandas2ri.activate()
import xarray as xr
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
from ark.analysis import (neighborhood_analysis, spatial_analysis_utils,
                          visualize)

import sys
sys.path.insert(0, os.path.abspath('../../'))
from SOTIP.sotip import *

def pairwise_enrichment(adata = None,
                        patient_key = 'Patient_ID',
                        cluster_key = 'cell_type',
                        n_neighs = None,
                        radius = None,
                        delaunay = True, 
                        condition_key = 'condition',
                        condition = 'cancer_core',
                        coord_type = 'generic'):
    adata_subset = adata[adata.obs[condition_key] == condition].copy()
    z_score_dict = {}
    count_dict = {}
    if n_neighs is None:
        n_neighs = 6
    for id in list(adata_subset.obs[patient_key].cat.categories):
        adata_run = adata_subset[np.isin(adata_subset.obs[patient_key], id)].copy()
        sq.gr.spatial_neighbors(adata_run, library_key = patient_key, delaunay = delaunay, coord_type = coord_type, n_neighs = n_neighs, radius = radius)
        sq.gr.nhood_enrichment(adata_run, cluster_key=cluster_key, show_progress_bar=False)
        #sq.pl.nhood_enrichment(adata_run, cluster_key=cluster_key, cmap = 'RdBu_r')
        z_score_dict[id] = adata_run.uns[f'{cluster_key}_nhood_enrichment']['zscore']
        count_dict[id] = adata_run.uns[f'{cluster_key}_nhood_enrichment']['count']
        
    return z_score_dict, count_dict

def pairwise_significance(adata,
                        z_score_dict_cond1 = None,
                        z_score_dict_cond2 = None,
                        count_dict_cond1 = None,
                        count_dict_cond2 = None,
                        cluster_key = 'mask_name',
                        save_directory = os.path.join('figures', 'simulated'),
                        ylim = None,
                        condition_list = None,
                        filename_save = 'simulated'):
    cell_types = list(adata.obs[cluster_key].cat.categories)
    
    # list of combinations for cell types, including self-interactions
    cell_combinations = [(cell_types[i], cell_types[j]) for i in range(len(cell_types)) for j in range(len(cell_types))]

    # min and max z-score for shared axis lims
    Z_min = np.floor(np.min([np.min(list(z_score_dict_cond1.values())), np.min(list(z_score_dict_cond2.values()))]))-50
    Z_max = np.ceil(np.max([np.max(list(z_score_dict_cond1.values())), np.max(list(z_score_dict_cond2.values()))]))+50

    fig, axes = plt.subplots(ncols=len(cell_combinations)//len(cell_types), nrows=len(cell_types), figsize=(2.5*len(cell_types), 2.25*len(cell_types)))

    lfc_arr =[]
    pval_arr = []
    cell_type_list  = []
    for idx, (cell_type1, cell_type2) in enumerate(cell_combinations):
        row = idx // len(cell_types)
        col = idx % len(cell_types)
        if col >= row:  # lower triangle
            cell_type_list.append(cell_type1+'__'+cell_type2)

            # z-scores for the given cell types from each matrix
            Z_A = [matrix[cell_types.index(cell_type1), cell_types.index(cell_type2)] for matrix in list(z_score_dict_cond1.values())]
            Z_B = [matrix[cell_types.index(cell_type1), cell_types.index(cell_type2)] for matrix in list(z_score_dict_cond2.values())]
            # counts
            C_A = [matrix[cell_types.index(cell_type1), cell_types.index(cell_type2)] for matrix in list(count_dict_cond1.values())]
            C_B = [matrix[cell_types.index(cell_type1), cell_types.index(cell_type2)] for matrix in list(count_dict_cond2.values())]
            
            # Wilcoxon rank-sum test for counts
            _, p_value = stats.ranksums(C_A, C_B)
            pval_arr.append(p_value)
            
            # log fold change
            mean_A = np.mean(C_A)
            mean_B = np.mean(C_B)
            log_fold_change = np.log2((mean_A / mean_B))
            lfc_arr.append(log_fold_change)
                        
            # plot
            df = pd.DataFrame({f'{condition_list[0]}': Z_A, f'{condition_list[1]}': Z_B})
            g = sns.violinplot(x='variable', y='value', data=df.melt(), ax=axes[row, col], inner='point', hue='variable', palette='Set1')
            g.set_xlabel('', fontsize=12)
            g.set_ylabel('Z-score', fontsize=12)
            g.tick_params(labelsize=10)
            if ylim is not None:
                g.set_ylim(Z_min, Z_max)
            try:
                axes[row, col].get_legend().remove()
            except: 
                pass
            axes[row, col].set_title(f'{cell_type1}__{cell_type2}')
        else:
            # Hide axes for upper triangle
            axes[row, col].axis('off')

    # Adjust layout
    plt.tight_layout()
    if filename_save is not None:
        plt.savefig(os.path.join(save_directory, filename_save+'.pdf'), bbox_inches = 'tight')
    else:
        plt.show()
    #multiple hypothesis testing
    adj_pval_arr = multipletests(pval_arr, method='fdr_bh')[1] 
    scores_df = pd.DataFrame(lfc_arr, index = cell_type_list, columns = ['logFC'])
    scores_df['pval'] = pval_arr
    scores_df['adj_pval'] = adj_pval_arr
    return scores_df

def evaluate_pairwise(adata = None,
                        patient_key = 'Patient_ID',
                        cluster_key = 'mask_name',
                        n_neighs = None,
                        condition_key = 'condition',
                        radius = None,
                        delaunay = True, 
                        coord_type = 'generic',
                        save_directory = os.path.join('figures', 'simulated'),
                        ylim = None,
                        condition_list = ['cancer_core', 'cancer_border'],
                        sig_threshold = 0.05, 
                        filename_save = 'simulated'):
    
    z_score_dict_cond1, count_dict_cond1 = pairwise_enrichment(adata = adata, patient_key = patient_key, cluster_key = cluster_key, condition_key = condition_key,
                                                               condition = condition_list[0], n_neighs = n_neighs, radius = radius, coord_type = coord_type, delaunay = delaunay)
    z_score_dict_cond2, count_dict_cond2 = pairwise_enrichment(adata = adata, patient_key = patient_key, cluster_key = cluster_key, condition_key = condition_key,
                                                               condition = condition_list[1], n_neighs = n_neighs, radius = radius, coord_type = coord_type, delaunay = delaunay)
    
    scores_df = pairwise_significance(adata, z_score_dict_cond1 = z_score_dict_cond1, z_score_dict_cond2 = z_score_dict_cond2,
                                    count_dict_cond1 = count_dict_cond1, count_dict_cond2 = count_dict_cond2, cluster_key = cluster_key,
                                    condition_list = condition_list, ylim = ylim, save_directory = save_directory, filename_save = filename_save)
    
    scores_df = scores_df[scores_df['pval'] <= sig_threshold]
    
    return scores_df

def graphcompass_interaction(adata,
                             patient_key = 'Patient_ID',
                             cluster_key = 'cell_cluster',
                             spatial_key = 'X_spatial',
                             coord_type = 'generic',
                            n_neighs = None,
                            radius = None,
                            delaunay = True):
    
    if n_neighs is None:
        n_neighs = 6
    
    adata.obs[patient_key] = pd.Categorical(adata.obs[patient_key])
    adata.obs[cluster_key] = pd.Categorical(adata.obs[cluster_key])
        
    sq.gr.spatial_neighbors(adata, library_key = patient_key, spatial_key = spatial_key, coord_type=coord_type,
                            radius=radius, n_neighs = n_neighs, delaunay = delaunay)
    
    adata.obs['library_id'] = adata.obs[patient_key].copy()
    cell_type_levels = adata.obs[cluster_key].cat.categories

    count_list = []
    for i, name in tqdm(enumerate(adata.obs_names)):
        row, col = adata.obsp['spatial_connectivities'][i, :].nonzero()
        count = adata.obs[cluster_key][col].value_counts()
        count_list.append(count)

    neighborhood_composition = pd.DataFrame(count_list, index=adata.obs_names)
    adata.uns['neighborhood_composition'] = neighborhood_composition

    zscore_list = []
    count_list = []
    celltype_names = []

    for i in adata.obs.library_id.cat.categories:
        adata_sub = adata[adata.obs.library_id == i].copy() #for each patient 
        adata_sub.obs[cluster_key] = adata_sub.obs[cluster_key].astype('category')
        ##compute pairwise enrichment between all cell types
        zscore, count = sq.gr.nhood_enrichment(adata_sub,
                                                cluster_key=cluster_key,
                                                copy=True,
                                                show_progress_bar=False)

        ct_labels = adata_sub.obs[cluster_key].cat.categories
        del adata_sub
        
        celltype_names.append(ct_labels)
        zscore_list.append(zscore)
        count_list.append(count)

    cell_type_combinations = pd.DataFrame()

    cell_types = pd.Series(cell_type_levels)
    cell_type_combinations['cell_type'] = cell_types.repeat(len(cell_types))
    cell_type_combinations['interactor'] = np.repeat([cell_types], len(cell_type_levels), axis=0).flatten()
    cell_type_combinations['combination'] = cell_type_combinations['cell_type'].add('__' + cell_type_combinations['interactor'])
    cell_type_combinations['dummy'] = cell_type_combinations['combination'].factorize()[0]
    n_celltypes = len(adata.obs[cluster_key].cat.categories)

    arr = np.array(cell_type_combinations.dummy).reshape(n_celltypes, n_celltypes)
    celltypexcelltype = pd.DataFrame(arr, index=adata.obs[cluster_key].cat.categories,
                                columns=adata.obs[cluster_key].cat.categories)

    df_list = []
    sample_id = adata.obs.library_id.cat.categories
    quant_data = zscore_list
    for i in range(len(quant_data)):
        df = pd.DataFrame()
        a = quant_data[i]
        upper = a[np.triu_indices(a.shape[0])] #take upper triangular portion of pairwise enrichment 
        values = np.array(upper.reshape(-1,1)) #flatten it into n x 1 array
        df['values'] = values.ravel()
        dummy_vars = celltypexcelltype.loc[celltypexcelltype.index.isin(celltype_names[i]),celltype_names[i]]
        dummy_vars = np.array(dummy_vars)
        dummy_vars = np.array(dummy_vars[np.triu_indices(dummy_vars.shape[0])].reshape(-1,1)).ravel()
        df['interaction_id'] = dummy_vars
        df.index = np.repeat(sample_id[i], df.shape[0]) 
        df_list.append(df)
        
    df_long = pd.concat(df_list) #concatenate by 0

    # replace dummy factors with cell type labels
    label_dict = {'interaction_id': dict(zip(cell_type_combinations['dummy'], cell_type_combinations['combination']))}
    df_long.replace(label_dict, inplace=True)

    # melt the matrix  
    df_wide = df_long.pivot(columns='interaction_id', values='values')
    df_wide[np.isnan(df_wide)] = 0

    adata.uns['celltype_interaction'] = df_wide
    return adata
    
def evaluate_graphcompass(adata,
                        patient_key = 'Patient_ID',
                        condition_key = 'Status',
                        cluster_key = 'mask_name',
                        spatial_key = 'spatial',
                        coord_type = 'generic',
                        n_neighs = None,
                        radius = None,
                        delaunay = True,
                        feature_key = 'spatial_nhood',
                        sig_threshold = 0.05):
    
    adata = graphcompass_interaction(adata, patient_key = patient_key, cluster_key = cluster_key, spatial_key = spatial_key,
                                     coord_type = coord_type, n_neighs = n_neighs, radius = radius, delaunay = delaunay)
    
    group = adata.obs[['library_id', patient_key, condition_key]].copy()
    group.drop_duplicates(inplace=True)
    group.set_index('library_id', inplace=True)

    interaction_mat = adata.uns['celltype_interaction']
    interaction_mat['condition'] = group[condition_key].astype('category').cat.reorder_categories(['cancer_border', 'cancer_core'], ordered=True)
    interaction_mat['subject_id'] = group[patient_key].astype('category')
    interaction_mat = pd.melt(interaction_mat, id_vars = ['condition','subject_id'])

    model = sm.formula.ols(formula='value ~  subject_id + interaction_id * condition', data=interaction_mat).fit()
    signif_coefs = model.params[model.pvalues <= sig_threshold]
    sig_niches = pd.DataFrame({'coef_label': signif_coefs.index, 'pval': model.pvalues[model.pvalues <= sig_threshold]})
    sig_niches = sig_niches[sig_niches['coef_label'].str.contains('condition')]
    sig_niches.index = sig_niches.index.str.extract(r'interaction_id\[T\.([^\]]+)\]:condition')[0].values
    sig_niches = sig_niches[~sig_niches.index.isna()]
    sig_niches.drop(columns = 'coef_label', inplace = True)

    mdata = MuData({'expression': adata, feature_key: anndata.AnnData(pd.DataFrame(adata.uns['neighborhood_composition']))})
    return mdata, sig_niches

def compute_kmeans_inertia(neighbor_mat_data, min_k=2, max_k=10, seed=42):
    """For a given neighborhood matrix, cluster and compute inertia using k-means clustering
       from the range of k=min_k to max_k

    Args:
        neighbor_mat_data (pandas.DataFrame):
            neighborhood matrix data with only the desired fovs
        min_k (int):
            the minimum k we want to generate cluster statistics for, must be at least 2
        max_k (int):
            the maximum k we want to generate cluster statistics for, must be at least 2
        seed (int):
            the random seed to set for k-means clustering

    Returns:
        xarray.DataArray:
            contains a single dimension, `cluster_num`, which indicates the inertia
            when `cluster_num` was set as k for k-means clustering
    """

    # create array we can store the results of each k for clustering
    coords = [np.arange(min_k, max_k + 1)]
    dims = ["cluster_num"]
    stats_raw_data = np.zeros(max_k - min_k + 1)
    cluster_stats = xr.DataArray(stats_raw_data, coords=coords, dims=dims)

    # iterate over each k value
    pb_format = '{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}]'
    for n in tqdm(range(min_k, max_k + 1), bar_format=pb_format):
        cluster_fit = KMeans(n_clusters=n, random_state=seed, n_init='auto').fit(neighbor_mat_data)
        cluster_stats.loc[n] = cluster_fit.inertia_

    return cluster_stats

def compute_kmeans(adata_niche,
                    unidentified_idx,
                    n_clusters = 2,
                    random_state = 0,
                    key_added = 'kmeans_cluster'):

    bool_idx = ~np.isin(adata_niche.obs_names, unidentified_idx)
    
    df = adata_niche[~np.isin(adata_niche.obs_names, unidentified_idx)].to_df()
    cluster_fit = KMeans(n_clusters=n_clusters, random_state=random_state, n_init='auto').fit_predict(df)
    cluster_fit = cluster_fit + 1 #add 1 to ensure labels start at 1 and not 0. this ensures our plotting functions will still work as currently written
    cluster_df = pd.DataFrame(cluster_fit, index = df.index, columns = [key_added])
    if key_added in adata_niche.obs.columns:
        adata_niche.obs.drop(columns = key_added, inplace = True)

    adata_niche.obs = pd.merge(adata_niche.obs, cluster_df, left_index = True, right_index = True, how = 'left')
    adata_niche.obs[key_added][bool_idx] = np.round(adata_niche.obs[key_added][bool_idx])
    adata_niche.obs[key_added][~bool_idx] = np.nan
    adata_niche.obs[key_added] = pd.Categorical(adata_niche.obs[key_added])

    return adata_niche

def compute_microenv_statistics(adata_niche,
                                fov_key = 'Patient_ID',
                                condition_key = 'condition',
                                labels_key = 'kmeans_cluster'):

    cell_sum = adata_niche.obs[[fov_key, labels_key]].groupby(by = fov_key).count()
    cell_sum = cell_sum.rename(columns={labels_key: 'cells_in_image'})     

    # create df will all fovs and all cluster rows
    all_df = []
    nclusters = len(adata_niche.obs[labels_key].unique())
    for fov in list(adata_niche.obs[fov_key].unique()):

        df = pd.DataFrame({
            fov_key: [fov] * nclusters,
            labels_key: list(range(1, nclusters+1))
        })
        
        all_df.append(df)
    all_df = pd.concat(all_df)

    # get cluster counts per image
    cluster_prop = adata_niche.obs[[fov_key, labels_key, 'label']].groupby(by=[fov_key, labels_key]).count().reset_index()

    # merge, replace nan with zero 
    cluster_prop = all_df.merge(cluster_prop, on=[fov_key, labels_key], how='right')
    if cluster_prop.isna().any().any():
        cluster_prop.fillna(0, inplace=True)

    cluster_prop = cluster_prop.rename(columns={'label': 'cells_in_cluster'})

    # # calculate proportions
    cluster_prop = cluster_prop.merge(cell_sum, on=[fov_key])
    cluster_prop['proportion'] = cluster_prop.cells_in_cluster / cluster_prop.cells_in_image
    cluster_prop = pd.merge(cluster_prop, adata_niche.obs[[fov_key, condition_key]].drop_duplicates(), on = [fov_key], how = 'left')

    return cluster_prop

def compute_microenv_significance(cluster_stats,
                                fov_key = 'Patient_ID',
                                condition_key = 'condition',
                                labels_key = 'kmeans_cluster',
                                condition_list = ['cancer_core', 'cancer_border']):
                                
    n_clusters = len(cluster_stats[labels_key].unique())
    z_score = StandardScaler(with_mean = True, with_std = True).fit_transform(cluster_stats['proportion'].values.reshape(-1,1))
    cluster_stats['normalized_value'] = z_score
    melted_df = cluster_stats.groupby([fov_key, labels_key])['normalized_value'].mean().unstack().reset_index()
    melted_df = pd.merge(melted_df, cluster_stats[[fov_key, condition_key]], on = fov_key, how = 'left')
    melted_df = melted_df.melt(ignore_index = False, id_vars = [fov_key, condition_key])

    _, axes = plt.subplots(1, n_clusters, figsize = (3*n_clusters, 2), gridspec_kw={'hspace': 0.45, 'wspace': 0.5, 'bottom':0.15})
    unique_clusters = list(melted_df['variable'].unique())
    pval_arr = []
    for i, ax in zip(range(0, n_clusters), axes.flat):
        df_run = melted_df[melted_df['variable'] == unique_clusters[i]].copy()
        group1 = df_run[df_run[condition_key] == condition_list[0]]['value'].values
        group2 = df_run[df_run[condition_key] == condition_list[1]]['value'].values
        _, p_value = stats.ranksums(group1, group2)
        pval_arr.append(p_value)
        g = sns.violinplot(x=condition_key, y='value', data=df_run, ax=ax, inner='point', palette='Set1', order = [condition_list[0], condition_list[1]])
        g.set_xlabel('', fontsize=10)
        g.set_ylabel('Z-score (avg. freq)', fontsize=10)
        g.tick_params(labelsize=10)
        g.set_title('cluster ' + str(unique_clusters[i]), fontsize = 10)

    #multiple hypothesis testing
    adj_pval_arr = multipletests(pval_arr, method='fdr_bh')[1] 
    scores_df = pd.DataFrame(pval_arr, index = unique_clusters, columns = ['pval'])
    scores_df['adj_pval'] = adj_pval_arr
    return scores_df

def evaluate_kmeans(adata,
                    n_clusters = 2,
                    random_state = 0,
                    fov_key = 'Patient_ID',
                    condition_key = 'condition',
                    labels_key = 'mask_name',
                    save_directory = os.path.join('figures', 'simulated'),
                    condition_list = ['cancer_core', 'cancer_border'],
                    filename_save = 'simulated',
                    radius = 200,
                    n_neighbors = None,
                    spatial_key = 'spatial', 
                    coord_type = 'generic',
                    delaunay = True,
                    min_cells = 3,
                    feature_key = 'spatial_nhood',
                    sig_threshold = 0.05,
                    nlargest = 3,
                    n_jobs = -1):

    adata = qu.tl.compute_spatial_neighbors(adata, radius = radius, n_neighbors = n_neighbors, spatial_key = spatial_key, delaunay = delaunay, fov_key = fov_key, coord_type = coord_type)
    adata_niche, cells_nonn = qu.tl.compute_niche_composition(adata, labels_key = labels_key, min_cells = min_cells)
    neighbor_inertia = compute_kmeans_inertia(adata_niche[~np.isin(adata_niche.obs_names, cells_nonn)].to_df(), min_k=2, max_k=10, seed=42)
    visualize.visualize_neighbor_cluster_metrics(neighbor_inertia, metric_name="inertia")
    adata_niche = compute_kmeans(adata_niche, cells_nonn, n_clusters = n_clusters, random_state=random_state, key_added = f'kmeans_cluster')    

    df = adata_niche.obs.groupby([f'kmeans_cluster', labels_key])['label'].count().unstack()
    df = df.div(df.sum(1), axis = 0)

    annotations = qu.tl.compute_niche_abundance(df, nlargest = nlargest)
    annotations = dict(zip(annotations.index, annotations.index.astype('str') + '_' + annotations.values))
    adata_niche.obs[f'kmeans_cluster_labeled'] = adata_niche.obs[f'kmeans_cluster'].map(annotations)

    sc.pl.matrixplot(adata_niche, var_names = adata_niche.var_names, groupby = f'kmeans_cluster_labeled', colorbar_title = 'avg. frequency', cmap = 'Purples', vmin = 0, vmax = 1, save = filename_save+'.pdf')
    cluster_stats = compute_microenv_statistics(adata_niche, fov_key = fov_key, condition_key = condition_key, labels_key = f'kmeans_cluster_labeled')
    scores_df = compute_microenv_significance(cluster_stats, fov_key = fov_key, condition_key = condition_key, labels_key = f'kmeans_cluster_labeled', condition_list = condition_list)
    scores_df = scores_df[scores_df['pval'] <= sig_threshold]
    mdata = MuData({'expression': adata, feature_key: adata_niche})
    return mdata, scores_df

def compute_MED_sotip(adata, spatial_key = 'spatial', n_neighbors = 20, labels_key = 'mask_name', patient_key = 'Patient_ID', threshold = 0):
    sc.pp.neighbors(adata, n_neighbors = n_neighbors, use_rep = 'X')
    sc.tl.umap(adata)
    sc.tl.paga(adata, groups = labels_key)
    sc.pl.paga(adata, threshold = threshold)
    sc.pl.paga_compare(adata, basis='X_umap', threshold = threshold)
    G = MED_multi(adata, use_cls = labels_key, nn = n_neighbors, copy = False, spatial_var = spatial_key, batch_obs = patient_key)
    composition_df = pd.DataFrame(G, index = adata.obs_names, columns = adata.obs[labels_key].cat.categories)
    adata_niche = anndata.AnnData(composition_df.div(composition_df.sum(axis = 1), axis = 0))
    return adata, adata_niche, G

def plot_matrix(mat, filename_save = None):
    fig, ax = plt.subplots(figsize=(2,2))  
    sns.heatmap(mat,annot=True,ax=ax)
    plt.gca().set_aspect('equal', adjustable='box') 
    if filename_save is not None:
        plt.savefig(filename_save+'pdf', bbox_inches = 'tight')

def compute_ME_sotip(adata, 
                     method = 'paga_guided_umap',
                     labels_key = 'mask_name',
                     n_jobs = -1,
                     n_neighbors = 800):
    
    #compute ground distance
    gd = get_ground_distance(adata, method = method, cls_key = labels_key) 
    gd_df = pd.DataFrame(gd, index = list(adata.obs[labels_key].cat.categories), columns = list(adata.obs[labels_key].cat.categories)) 
    plot_matrix(gd_df, filename_save='ground_distance_example')
    adata.uns['GD'] = gd
    #pairwise distances using pyemd
    adata_phEMD = MED_phEMD_mp(adata.copy(), GD_method = method, MED_knn = 10, CT_obs = labels_key, ifspatialplot = False, OT_method = 'pyemd',
                               ME_precompyted = True, GD_precomputed = True, mp = n_jobs)
    adata.obsp['ME_EMD_mat'] = adata_phEMD.obsm['X_ME_EMD_mat']
    sc.pp.neighbors(adata, n_neighbors = n_neighbors, use_rep = 'X')
    knn_indices, knn_dists, forest = qu.tl.compute_neighbors_umap(adata.obsp['ME_EMD_mat'], n_neighbors = n_neighbors, metric = 'precomputed')
    adata.obsp['distances'], adata.obsp['connectivities'] = qu.tl.compute_connectivities_umap(knn_indices, knn_dists, adata.shape[0], n_neighbors)
    adata.uns['neighbors_EMD'] = adata.uns['neighbors'].copy()
    return adata

def merge_cls_paga(adata,
                   cls_key = 'SOTIP_cluster',
                   neighbor_key = 'neighbors_EMD',
                   thresh = 0.5,
                   min_cls = 2,
                   paga_plot = True):
    
    if f'{cls_key}_merge_colors' in adata.uns:
        del adata.uns[f'{cls_key}_merge_colors']
    adata.obs[f'{cls_key}_merge'] = adata.obs[cls_key].copy()
    merge_cls_list = []
    while True:
        sc.tl.paga(adata,groups=f'{cls_key}_merge',neighbors_key=neighbor_key)
        
        cur_conn = adata.uns['paga']['connectivities'].toarray()
        cur_ME_cls = adata.obs[f'{cls_key}_merge'].copy()
        if len(cur_ME_cls.cat.categories)<=min_cls:
            break
        merge_cls_idx = np.unravel_index(np.argmax(cur_conn),cur_conn.shape)
        # print(merge_cls_idx)
        if cur_conn[merge_cls_idx]<thresh:
            break
        merge_cls_i = cur_ME_cls.cat.categories[merge_cls_idx[0]]
        merge_cls_j = cur_ME_cls.cat.categories[merge_cls_idx[1]]

        new_ME_cls = merge_cls(cur_ME_cls,merge_cls_i,merge_cls_j)
        new_ME_cls = new_ME_cls.cat.remove_unused_categories()
        merge_cls_list.append(new_ME_cls.copy())
        adata.obs[f'{cls_key}_merge'] = new_ME_cls
        
    adata.uns['merge_cls_list'] = merge_cls_list
    if paga_plot == True:
        sc.pl.paga(adata,color=[f'{cls_key}_merge'],threshold=0)

def merge_cls(old_cls,cls_i,cls_j):
    old_cls[old_cls==cls_i] = cls_j
    print(f'merged {cls_i} to {cls_j}')
    return old_cls

def cluster_sotip(adata, thresh = 0.2, n_clusters = 5, key_added = 'SOTIP_cluster'):
    sc.tl.umap(adata, neighbors_key = 'neighbors_EMD')
    adata.obsm['X_umap_EMD'] = adata.obsm['X_umap']
    sc.tl.leiden(adata, neighbors_key='neighbors_EMD', key_added = key_added)
    sc.tl.paga(adata, groups= key_added, neighbors_key='neighbors_EMD')
    merge_cls_paga(adata, thresh = thresh, min_cls = n_clusters, paga_plot = True, cls_key = key_added)

    return adata

def _add_metadata_sotip(adata, adata_niche):
    adata_niche.obs = adata.obs.copy()
    return adata_niche

def evaluate_sotip(adata,
                    n_clusters = 2,
                    n_neighbors = 30,
                    k_sim = 50,
                    fov_key = 'Patient_ID',
                    condition_key = 'condition',
                    labels_key = 'mask_name',
                    save_directory = os.path.join('figures', 'simulated'),
                    condition_list = ['cancer_core', 'cancer_border'],
                    filename_save = 'simulated',
                    spatial_key = 'spatial', 
                    sig_threshold = 0.05,
                    thresh = 0.2,
                    feature_key = 'spatial_nhood',
                    gd_method = 'paga_guided_umap',
                    nlargest = 3,
                    n_jobs = 8):

    adata.obs[fov_key] = pd.Categorical(adata.obs[fov_key])
    adata.obs[labels_key] = pd.Categorical(adata.obs[labels_key])

    adata, adata_niche, G = compute_MED_sotip(adata, spatial_key = spatial_key, n_neighbors = n_neighbors, labels_key = labels_key, patient_key = fov_key)
    adata = compute_ME_sotip(adata, method = gd_method, labels_key = labels_key, n_neighbors = k_sim, n_jobs = n_jobs)
    adata = cluster_sotip(adata, thresh = thresh, n_clusters = n_clusters)
    adata_niche = _add_metadata_sotip(adata, adata_niche)

    df = adata_niche.obs.groupby([f'SOTIP_cluster_merge', labels_key])['label'].count().unstack()
    df = df.div(df.sum(1), axis = 0)

    annotations = qu.tl.compute_niche_abundance(df, nlargest = nlargest)
    annotations = dict(zip(annotations.index, annotations.index.astype('str') + '_' + annotations.values))
    adata_niche.obs['SOTIP_cluster_merge_labeled'] = adata_niche.obs['SOTIP_cluster_merge'].map(annotations)

    sc.pl.matrixplot(adata_niche, var_names = adata_niche.var_names, groupby = f'SOTIP_cluster_merge_labeled', colorbar_title = 'avg. frequency', cmap = 'Purples', vmin = 0, vmax = 1, save = filename_save+'.pdf')
    cluster_stats = compute_microenv_statistics(adata_niche, fov_key = fov_key, condition_key = condition_key, labels_key = f'SOTIP_cluster_merge_labeled')
    scores_df = compute_microenv_significance(cluster_stats, fov_key = fov_key, condition_key = condition_key, labels_key = f'SOTIP_cluster_merge_labeled', condition_list = condition_list)
    scores_df = scores_df[scores_df['pval'] <= sig_threshold]
    mdata = MuData({'expression': adata, feature_key: adata_niche})
    return mdata, scores_df

def run_milo(adata,
             n_neighbors = 100,
             design = '~condition',
             model_contrasts = 'conditionA-conditionB',
             patient_key = 'Patient_ID',
             cluster_key = 'mask_name',
             feature_key = 'rna',
             sig_threshold = 0.05,
             prop = 0.1,
             **kwargs):
    import pertpy as pt 
    import scanpy as sc
    milo = pt.tl.Milo()
    mdata = milo.load(adata)
    sc.tl.pca(mdata[feature_key])
    sc.pp.neighbors(mdata[feature_key], n_neighbors = n_neighbors, use_rep = 'X_pca')
    # milo.make_nhoods(mdata[feature_key], prop = prop)
    mdata[feature_key].uns["nhood_neighbors_key"] = None
    mdata = qu.tl.build_milo_graph(mdata, feature_key=feature_key)
    mdata = milo.count_nhoods(mdata, sample_col = patient_key)
    milo.da_nhoods(mdata,
                design=design,
                model_contrasts = model_contrasts)
    milo.annotate_nhoods(mdata, anno_col = cluster_key, feature_key = feature_key)
    mdata = MuData({'expression': mdata[feature_key], 'milo': mdata['milo']})
    scores_df = pd.DataFrame(mdata['milo'].var.groupby('nhood_annotation')['SpatialFDR'].median()[mdata['milo'].var.groupby('nhood_annotation')['SpatialFDR'].median() < sig_threshold])
    scores_df.columns = ['pval']

    return mdata, scores_df