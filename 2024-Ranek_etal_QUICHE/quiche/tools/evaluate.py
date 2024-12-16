import pandas as pd
import numpy as np
import os
import multiprocessing as mp
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import matplotlib.pyplot as plt
import seaborn as sns
import quiche as qu
from sketchKH import *
import anndata
import cellcharter as cc
from muon import MuData
import squidpy as sq
from tqdm import tqdm
from scipy import stats
import pertpy as pt
import scanpy as sc
from statsmodels.stats.multitest import multipletests
import statsmodels.api as sm
import xarray as xr
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
from ark.analysis import visualize
import pertpy as pt 
import scanpy as sc

def pairwise_enrichment(adata = None,
                        patient_key = 'Patient_ID',
                        cluster_key = 'cell_cluster',
                        n_neighs = None,
                        radius = None,
                        delaunay = True, 
                        condition_key = 'condition',
                        condition = 'cancer_core',
                        coord_type = 'generic'):
    """Performs pairwise cell type enrichment (https://pubmed.ncbi.nlm.nih.gov/28783155) using squidpy framework: https://squidpy.readthedocs.io/en/stable/index.html

    Parameters
    adata: anndata.AnnData (default = None)
        annotated data object containing preprocessed single-cell data 
    patient_key: str (default = 'Patient_ID')
        string referring to column in adata.obs with fov-level information
    cluster_key: str (default = 'cell_cluster')
        string referring to column in adata.obs with cell phenotype information
    n_neighs: int (default = None)
        number of nearest neighbors
    radius: int (default = None)
        radius for bounding local niche detection
    delaunay: bool (default = False)
        boolean referring to whether delaunay triangulation should be used
    condition_key: str (default = 'condition') 
        string in adata.obs containing condition information
    condition: str (default = 'cancer_core')
        string referring to specific condition for testing across all samples
    coord_type: str (default = 'generic')
        string referring to spatial coordinated
    ----------

    Returns
    z_score_dict: dictionary
        dictionary containing cell type x cell type Z-scores
    count_dict: dictionary
        dictionary containing cell type x cell type counts
    ----------
    """
    adata_subset = adata[adata.obs[condition_key] == condition].copy()
    z_score_dict = {}
    count_dict = {}
    if n_neighs is None:
        n_neighs = 6
    for id in list(adata_subset.obs[patient_key].cat.categories):
        adata_run = adata_subset[np.isin(adata_subset.obs[patient_key], id)].copy()
        sq.gr.spatial_neighbors(adata_run, library_key = patient_key, delaunay = delaunay, coord_type = coord_type, n_neighs = n_neighs, radius = radius)
        sq.gr.nhood_enrichment(adata_run, cluster_key=cluster_key, show_progress_bar=False)
        z_score_dict[id] = adata_run.uns[f'{cluster_key}_nhood_enrichment']['zscore']
        count_dict[id] = adata_run.uns[f'{cluster_key}_nhood_enrichment']['count']
        
    return z_score_dict, count_dict

def pairwise_significance(adata,
                        z_score_dict_cond1 = None,
                        z_score_dict_cond2 = None,
                        count_dict_cond1 = None,
                        count_dict_cond2 = None,
                        cluster_key = 'cell_cluster',
                        save_directory = os.path.join('figures', 'simulated'),
                        ylim = None,
                        condition_list = None,
                        filename_save = 'simulated'):
    """Computes niches according to khop spatial neighborhood

    Parameters
    adata: anndata.AnnData (default = None)
        annotated data object containing preprocessed single-cell data 
    z_score_dict_cond1: dictionary (default = None)
        dictionary containing cell type x cell type Z-scores for condition 1
    z_score_dict_cond2: dictionary (default =  None)
        dictionary containing cell type x cell type Z-scores for condition 2
    count_dict_cond1: dictionary (default = None)
        dictionary containing cell type x cell type counts for condition 1
    count_dict_cond2: dictionary (default = None)
        dictionary containing cell type x cell type counts for condition 2
    cluster_key: str (default = 'cell_cluster')
        string referring to the column in adata.obs that contains cell phenotype labels
    condition_list: list (default = None) 
        list containing condition information
    save_directory: str (default =  os.path.join('figures', 'simulated'))
        string specifying path for saving plots
    ylim: list (default = None)
        y-axis limits for plotting
    filename_save: str (default = 'simulated)
        str specifying name of plot for saving
    ----------

    Returns
    scores_df: pd.DataFrame
        dataframe containing pairwise significance scores
    ----------
    """
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
            axes[row, col].axis('off')

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
                        cluster_key = 'cell_cluster',
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
    """Performs pairwise cell type enrichment (https://pubmed.ncbi.nlm.nih.gov/28783155) using squidpy framework: https://squidpy.readthedocs.io/en/stable/index.html

    Parameters
    adata: anndata.AnnData (default = None)
        annotated data object containing preprocessed single-cell data 
    patient_key: str (default = 'Patient_ID')
        string referring to column in adata.obs with fov-level information
    cluster_key: str (default = 'cell_cluster')
        string referring to column in adata.obs with cell phenotype information
    n_neighs: int (default = None)
        number of nearest neighbors
    radius: int (default = None)
        radius for bounding local niche detection
    delaunay: bool (default = False)
        boolean referring to whether delaunay triangulation should be used
    condition_key: str (default = 'condition') 
        string in adata.obs containing condition information
    condition: str (default = 'cancer_core')
        string referring to specific condition for testing across all samples
    condition_list: list (default = None) 
        list containing condition information
    coord_type: str (default = 'generic')
        string referring to spatial coordinated
    save_directory: str (default =  os.path.join('figures', 'simulated'))
        string specifying path for saving plots
    ylim: list (default = None)
        y-axis limits for plotting
    filename_save: str (default = 'simulated)
        str specifying name of plot for saving
    sig_threshold: float (default = 0.05)
        threshold for significance
    ----------

    Returns
    scores_df: pd.DataFrame
        dataframe containing pairwise significance scores
    ----------
    """
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
    """Performs spatial enrichment analysis using GraphCompass: https://github.com/theislab/graphcompass

    Parameters
    adata: anndata.AnnData (default = None)
        annotated data object containing preprocessed single-cell data 
    patient_key: str (default = 'Patient_ID')
        string referring to column in adata.obs with fov-level information
    cluster_key: str (default = 'cell_cluster')
        string referring to column in adata.obs with cell phenotype information
    radius: int (default = 200)
        integer referring to the radius in pixels for bounding local niche detection
    spatial_key: str (default = 'spatial') 
        string in adata.obsm containing cell centroid coordinates
    coord_type: str (default = 'generic')
        string referring to spatial coordinated
    n_neighs: int (default = None)
        number of nearest neighbors
    delaunay: bool (default = False)
        boolean referring to whether delaunay triangulation should be used
    ----------

    Returns
    adata: anndata.AnnData
        annotated data object containing spatial enrichment analysis
    ----------
    """
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
        
    df_long = pd.concat(df_list)

    label_dict = {'interaction_id': dict(zip(cell_type_combinations['dummy'], cell_type_combinations['combination']))}
    df_long.replace(label_dict, inplace=True)

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
                        condition_list = ['cancer_border', 'cancer_core'],
                        feature_key = 'spatial_nhood',
                        sig_threshold = 0.05):
    """Performs spatial enrichment analysis using GraphCompass: https://github.com/theislab/graphcompass

    Parameters
    adata: anndata.AnnData (default = None)
        annotated data object containing preprocessed single-cell data 
    patient_key: str (default = 'Patient_ID')
        string referring to column in adata.obs with fov-level information
    cluster_key: str (default = 'cell_cluster')
        string referring to column in adata.obs with cell phenotype information
    radius: int (default = 200)
        integer referring to the radius in pixels for bounding local niche detection
    spatial_key: str (default = 'spatial') 
        string in adata.obsm containing cell centroid coordinates
    coord_type: str (default = 'generic')
        string referring to spatial coordinated
    n_neighs: int (default = None)
        number of nearest neighbors
    delaunay: bool (default = False)
        boolean referring to whether delaunay triangulation should be used
    condition_key: str (default = 'condition') 
        string in adata.obs containing condition information
    feature_key: str (default = 'spatial_nhood')
        string specifying where spatial enrichment analysis should be stored in mdata object
    condition_list: list (default = None) 
        list containing condition information
    sig_threshold: float (default = 0.05)
        threshold for significance
    ----------

    Returns
    mdata: mudata object
        annotated data object containing spatial enrichment analysis
    sig_niches: pd.DataFrame
        dataframe containing pairwise significance scores
    ----------
    """ 
    adata = graphcompass_interaction(adata, patient_key = patient_key, cluster_key = cluster_key, spatial_key = spatial_key,
                                     coord_type = coord_type, n_neighs = n_neighs, radius = radius, delaunay = delaunay)
    
    group = adata.obs[['library_id', patient_key, condition_key]].copy()
    group.drop_duplicates(inplace=True)
    group.set_index('library_id', inplace=True)

    interaction_mat = adata.uns['celltype_interaction']
    interaction_mat['condition'] = group[condition_key].astype('category').cat.reorder_categories(condition_list, ordered=True)
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

def compute_kmeans_inertia(neighbor_mat_data,
                           min_k = 2,
                           max_k = 10,
                           seed = 42):
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
    """Performs unsupervised clustering using KMeans++: https://scikit-learn.org/stable/modules/generated/sklearn.cluster.KMeans.html

    Parameters
    adata_niche: anndata.AnnData (default = None)
        annotated data object containing niche information (niche x cell type)
    unidentified_idx: list (default = None)
        list containing cells that don't pass the min cell cutoff
    n_clusters: int (default = 2)
        integer referring to the number of spatial clusters
    random_state: int (default = 0)
        integer referring to the random state parameter
    key_added: str (default = 'kmeans_cluster')
        string referring to the column in adata_niche.obs to save spatial cluster information
    ----------

    Returns
    adata_niche: anndata.AnnData
        annotated data object containing niche information and annotations
    ----------
    """
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
    """Computes the proportion of spatial clusters in each sample

    Parameters
    adata_niche: anndata.AnnData (default = None)
        annotated data object containing niche information (niche x cell type)
    fov_key: str (default = 'Patient_ID')
        string referring to the column in adata_niche.obs with patient-level annotations
    condition_key: str (default = 'condition')
        string referring to the column in adata_niche.obs with condition-level annotations
    labels_key: str (default = 'kmeans_cluster')
        string referring to the column in adata_niche.obs with spatial cluster annotations
    ----------

    Returns
    cluster_prop: pd.DataFrame
        dataframe containing the proportion of spatial clusters in each sample
    ----------
    """
    cell_sum = adata_niche.obs[[fov_key, labels_key]].groupby(by = fov_key).count()
    cell_sum = cell_sum.rename(columns={labels_key: 'cells_in_image'})     

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

    # calculate proportions
    cluster_prop = cluster_prop.merge(cell_sum, on=[fov_key])
    cluster_prop['proportion'] = cluster_prop.cells_in_cluster / cluster_prop.cells_in_image
    cluster_prop = pd.merge(cluster_prop, adata_niche.obs[[fov_key, condition_key]].drop_duplicates(), on = [fov_key], how = 'left')

    return cluster_prop

def compute_microenv_significance(cluster_stats,
                                fov_key = 'Patient_ID',
                                condition_key = 'condition',
                                labels_key = 'kmeans_cluster',
                                condition_list = ['cancer_core', 'cancer_border']):
    """Wilcoxon rank sum test for spatial clustering approaches

    Parameters
    cluster_stats: pd.DataFrame (default = None)
        dataframe containing spatial cluster proportion information for each sample
    fov_key: str (default = 'Patient_ID')
        string referring to the column with patient-level annotations
    condition_key: str (default = 'condition')
        string referring to the column with condition-level annotations
    labels_key: str (default = 'kmeans_cluster')
        string referring to the column with spatial cluster annotations
    condition_list: list (default = ['cancer_core', 'cancer_border'])
        list containing conditions
    ----------

    Returns
    scores_df: pd.DataFrame
        dataframe containing spatial cluster significance scores
    ----------
    """                              
    n_clusters = len(cluster_stats[labels_key].unique())
    melted_df = cluster_stats.groupby([fov_key, labels_key])['proportion'].mean().unstack().reset_index()
    melted_df = pd.merge(melted_df, cluster_stats[[fov_key, condition_key]].drop_duplicates(), on = fov_key, how = 'left')
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
        g.set_ylabel('(avg. freq)', fontsize=10)
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
    """Performs spatial enrichment analysis using KMeans++ clustering: https://scikit-learn.org/stable/modules/generated/sklearn.cluster.KMeans.html

    Parameters
    adata: anndata.AnnData (default = None)
        annotated data object containing preprocessed single-cell data 
    fov_key: str (default = 'Patient_ID')
        string referring to column in adata.obs with fov-level information
    labels_key: str (default = 'cell_cluster')
        string referring to column in adata.obs with cell phenotype information
    n_clusters: int (default = 2)
        integer referring to the number of spatial clusters
    random_state: int (default = 0)
        integer referring to the random state parameter
    radius: int (default = 200)
        integer referring to the radius in pixels for bounding local niche detection
    spatial_key: str (default = 'spatial') 
        string in adata.obsm containing cell centroid coordinates
    coord_type: str (default = 'generic')
        string referring to spatial coordinated
    n_neighbors: int (default = None)
        number of nearest neighbors
    delaunay: bool (default = False)
        boolean referring to whether delaunay triangulation should be used
    condition_key: str (default = 'condition') 
        string in adata.obs containing condition information
    feature_key: str (default = 'spatial_nhood')
        string specifying where spatial enrichment analysis should be stored in mdata object
    sig_threshold: float (default = 0.05)
        threshold for significance
    condition_list: list (default = ['cancer_core', 'cancer_border'])
        list containing conditions
    nlargest: int (default = 3)
        integer referring to the number of cell types for labeling spatial clusters
    n_jobs: int (default = -1)
        number of tasks for parallelization
    ----------

    Returns
    mdata: mudata object
        annotated data object containing spatial enrichment analysis
    sig_niches: pd.DataFrame
        dataframe containing spatial cluster significance scores
    ----------
    """
    adata = qu.tl.compute_spatial_neighbors(adata, radius = radius, n_neighbors = n_neighbors, spatial_key = spatial_key, delaunay = delaunay, fov_key = fov_key, coord_type = coord_type)
    adata_niche, cells_nonn = qu.tl.compute_niche_composition(adata, labels_key = labels_key, min_cells = min_cells)
    neighbor_inertia = compute_kmeans_inertia(adata_niche[~np.isin(adata_niche.obs_names, cells_nonn)].to_df(), min_k=2, max_k=10, seed=42)
    visualize.visualize_neighbor_cluster_metrics(neighbor_inertia, metric_name="inertia")
    adata_niche = compute_kmeans(adata_niche, cells_nonn, n_clusters = n_clusters, random_state=random_state, key_added = f'kmeans_cluster')    

    df = adata_niche.obs.groupby([f'kmeans_cluster', labels_key])['label'].count().unstack()
    df = df.div(df.sum(1), axis = 0)

    annotations = qu.tl.compute_niche_abundance_fov(df, nlargest = nlargest)
    annotations = dict(zip(annotations.index, annotations.index.astype('str') + '_' + annotations.values))
    adata_niche.obs[f'kmeans_cluster_labeled'] = adata_niche.obs[f'kmeans_cluster'].map(annotations)

    sc.pl.matrixplot(adata_niche, var_names = adata_niche.var_names, groupby = f'kmeans_cluster_labeled', colorbar_title = 'avg. frequency', cmap = 'Purples', vmin = 0, vmax = 1, save = filename_save+'.pdf')
    cluster_stats = compute_microenv_statistics(adata_niche, fov_key = fov_key, condition_key = condition_key, labels_key = f'kmeans_cluster_labeled')
    scores_df = compute_microenv_significance(cluster_stats, fov_key = fov_key, condition_key = condition_key, labels_key = f'kmeans_cluster_labeled', condition_list = condition_list)
    scores_df = scores_df[scores_df['pval'] <= sig_threshold]
    mdata = MuData({'expression': adata, feature_key: adata_niche})
    mdata[feature_key].obs['pval'] = pd.Series(mdata[feature_key].obs['kmeans_cluster_labeled']).map(dict(zip(list(scores_df.index), scores_df.pval)))
    return mdata, scores_df

def evaluate_cell_charter(adata,
                    random_state = 0,
                    n_clusters = None,
                    fov_key = 'Patient_ID',
                    condition_key = 'condition',
                    batch_correct = False,
                    sig_threshold = 0.05,
                    n_jobs = -1, 
                    max_runs = 5,
                    n_layers=3,
                    condition_list = ['cancer_core', 'cancer_border']):
    """Performs spatial enrichment analysis using cellcharter: https://pubmed.ncbi.nlm.nih.gov/38066188/, https://cellcharter.readthedocs.io/en/latest/index.html

    Parameters
    adata: anndata.AnnData (default = None)
        annotated data object containing preprocessed single-cell data 
    fov_key: str (default = 'Patient_ID')
        string referring to column in adata.obs with fov-level information
    n_clusters: int (default = 2)
        integer referring to the number of spatial clusters
    random_state: int (default = 0)
        integer referring to the random state parameter
    condition_key: str (default = 'condition') 
        string in adata.obs containing condition information
    sig_threshold: float (default = 0.05)
        threshold for significance
    condition_list: list (default = ['cancer_core', 'cancer_border'])
        list containing conditions
    batch_correct: bool (default = False)
        boolean referring to whether batch effect correction should be performed
    max_runs: int (default = 5)
        integer referring to number of runs for automatically determining number of clusters
    n_layers: int (default = 3)
        integer for number of aggregation layers
    n_jobs: int (default = -1)
        number of tasks for parallelization
    ----------

    Returns
    mdata: mudata object
        annotated data object containing spatial enrichment analysis
    sig_niches: pd.DataFrame
        dataframe containing spatial cluster significance scores
    ----------
    """
    if n_jobs == -1:
        n_jobs = mp.cpu_count()
    
    scaler = StandardScaler()
    X = scaler.fit_transform(adata.X.copy())
    adata.X = X.astype('float32') 

    conditions = list(adata.obs[fov_key].unique())

    if batch_correct:

        early_stopping_kwargs = {
            "early_stopping_metric": "val_unweighted_loss",
            "threshold": 0,
            "patience": 20,
            "reduce_lr": True,
            "lr_patience": 13,
            "lr_factor": 0.1,
        }

        trvae = sca.models.TRVAE(adata=adata,
                                condition_key=fov_key,
                                conditions=conditions,
                                hidden_layer_sizes=[128, 128],
                                recon_loss='mse')
        
        adata = remove_sparsity(adata)

        trvae.train(
            n_epochs=500,
            alpha_epoch_anneal=200,
            early_stopping_kwargs=early_stopping_kwargs,
            recon_loss='mse')
        
        adata.obsm['X_embed'] = sc.AnnData(trvae.get_latent()).X
    else:
        sc.tl.pca(adata, n_comps = 50)
        adata.obsm['X_embed'] = adata.obsm['X_pca']

    sq.gr.spatial_neighbors(adata, library_key=fov_key, coord_type='generic', delaunay=True)
    cc.gr.remove_long_links(adata)
    cc.gr.aggregate_neighbors(adata, n_layers=n_layers, use_rep='X_embed')

    if n_clusters is None:
        model_params = {'random_state': random_state, 'trainer_params': {'accelerator':'cpu', 'enable_progress_bar': False, 'devices':n_jobs}}
        models = cc.tl.ClusterAutoK(n_clusters=(2,10), model_class=cc.tl.GaussianMixture, model_params=model_params, max_runs=max_runs)
        models.fit(adata, use_rep='X_cellcharter')
        pred= models.predict(adata, use_rep='X_cellcharter')
        if 'spatial_cluster' in adata.obs:
            del adata.obs['spatial_cluster']
        new_column = pd.DataFrame({'spatial_cluster': pred}, index=adata.obs.index).astype('category')
        adata.obs = pd.concat([adata.obs, new_column], axis=1)
    else:
        gmm = cc.tl.Cluster(
            n_clusters=n_clusters,
            random_state=random_state,
            trainer_params=dict(accelerator='cpu', devices=n_jobs))

        gmm.fit(adata, use_rep='X_cellcharter')
        pred = gmm.predict(adata, use_rep='X_cellcharter')
        if 'spatial_cluster' in adata.obs:
            del adata.obs['spatial_cluster']
        new_column = pd.DataFrame({'spatial_cluster': pred}, index=adata.obs.index).astype('category')
        adata.obs = pd.concat([adata.obs, new_column], axis=1)

    cluster_stats = compute_microenv_statistics(adata, fov_key = fov_key, condition_key = condition_key, labels_key = f'spatial_cluster')
    scores_df = compute_microenv_significance(cluster_stats, fov_key = fov_key, condition_key = condition_key, labels_key = f'spatial_cluster', condition_list = condition_list)
    mdata = MuData({'expression': adata})
    mdata['expression'].obs['pval'] = pd.Series(mdata['expression'].obs['spatial_cluster']).map(dict(zip(list(scores_df.index), scores_df.pval)))
    mdata['expression'].obs['pval'] = mdata['expression'].obs['pval'].values.astype('float')
    scores_df = scores_df[scores_df['pval'] <= sig_threshold]
    return mdata, scores_df

def run_milo(adata,
             n_neighbors = 100,
             design = '~condition',
             model_contrasts = 'conditionA-conditionB',
             patient_key = 'Patient_ID',
             cluster_key = 'cell_cluster',
             feature_key = 'rna',
             sig_threshold = 0.05,
             prop = 0.1,
             **kwargs):
    """Performs differential cell type abundance analysis using Milo: https://www.nature.com/articles/s41587-021-01033-z, https://pertpy.readthedocs.io/en/latest/usage/tools/pertpy.tools.Milo.html

    Parameters
    adata: anndata.AnnData (default = None)
        annotated data object containing preprocessed single-cell data 
    n_neighbors: int (default = 10)
        number of nearest neighbors
    design: str (default = '~condition')
        string referring to design
    model_contrasts: str (default = 'conditionA-conditionB')
        string for condition-specific testing
    patient_key: str (default = 'Patient_ID') 
        string in adata.obs containing patient-level information
    cluster_key: str (default = 'cell_cluster')
        string in adata.obs containing cell cluster annotations
    feature_key: str (default = 'rna')
        str referring to expression information
    sig_threshold: float (default = 0.05)
        threshold for significance
    ----------

    Returns
    mdata: mudata object
        annotated data object containing cell type abundance analysis
    scores_df: pd.DataFrame
        dataframe containing differential cell type abundance scores
    ----------
    """
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
    scores_df = pd.DataFrame(mdata['milo'].var.groupby('nhood_annotation')['SpatialFDR'].median()[mdata['milo'].var.groupby('nhood_annotation')['SpatialFDR'].median() <= sig_threshold])
    scores_df.columns = ['pval']
    mdata['milo'].var[mdata['expression'].obs.columns] = mdata['expression'].obs.values
    return mdata, scores_df