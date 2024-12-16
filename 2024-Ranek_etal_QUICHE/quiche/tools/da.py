import scanpy as sc
import pertpy as pt
import pandas as pd
import numpy as np
from scipy.sparse import csr_matrix
import scipy
from rpy2.robjects import pandas2ri
import anndata
import quiche as qu
from sketchKH import sketch
import logging
from muon import MuData
import matplotlib.pyplot as plt
pandas2ri.activate()

def run_quiche(adata,
            radius = 200,
            labels_key = 'cell_cluster',
            spatial_key = 'spatial',
            fov_key = 'Patient_ID',
            patient_key = 'Patient_ID',
            n_neighbors = None,
            delaunay = True,
            min_cells = 5,
            khop = None,
            coord_type = 'generic',
            sketch_key = 'Patient_ID',
            test_key = 'Patient_ID',
            sketch_size = None,
            frequency_seed = 0,
            gamma = 1,
            k_sim = 100,
            design = '~condition',
            model_contrasts = 'conditionA-conditionB',
            nlargest = 4,
            min_perc = 0.1,
            label_scheme = 'normal',
            merge = False,
            feature_key = 'spatial_nhood',
            annotation_key = 'quiche_niche',
            sig_key = 'SpatialFDR',
            sig_threshold = 0.05,
            n_jobs = -1,
            **kwargs):
    """Performs spatial enrichment analysis using QUICHE

    Parameters
    adata: anndata.AnnData (default = None)
        annotated data object containing preprocessed single-cell data 
    radius: int (default = 200)
        integer referring to the radius in pixels for bounding local niche detection
    labels_key: str (default = 'cell_cluster')
        string referring to column in adata.obs with cell phenotype information
    spatial_key: str (default = 'spatial') 
        string in adata.obsm containing cell centroid coordinates
    fov_key: str (default = 'Patient_ID')
        string referring to column in adata.obs with fov-level information
    n_neighbors: int (default = None)
        number of nearest neighbors for spatial proximity detection
    delaunay: bool (default = False)
        boolean referring to whether delaunay triangulation should be used
    min_cells: int (default = 0)
        integer referring to the minimum number of nearest neighbors for a niche cell type proportion vector to be considered
    khop: int (default = None)
        integer referring to number of hops*n_neighbors for local niche detection
    coord_type: str (default = 'generic')
        string referring to spatial coordinated
    sketch_key: str (default = 'Patient_ID')
        string in adata.obs containing patient-level information for downsampling
    test_key: str (default = 'Patient_ID')
        string in adata.obs for condition-specific testing
    sketch_size: int (default = None)
        integer referring to the number of niches to select per patient sample See: https://dl.acm.org/doi/10.1145/3535508.3545539, https://github.com/CompCy-lab/SketchKH.
        if None: defaults to the minimum within each sample
    frequency_seed: int (default = 0)
        integer referring to the random state parameter in downsampling
    gamma: int (default = 1)
        scale parameter for the normal distribution standard deviation in random Fourier frequency feature computation in downsampling
    k_sim: int (default = 100)
        number of nearest neighbors for niche similarity graph constrcution
    design: str (default = '~condition')
        string referring to design
    model_contrasts: str (default = 'conditionA-conditionB')
        string for condition-specific testing
    annotation_key: str (default = 'quiche_niche')
        string specifying column for niche annotations
    nlargest: int (default = 3)
        integer referring to the number of cell types for labeling spatial clusters
    label_scheme: str (default = 'normal')
        string specifying how niche neighborhoods should be annotated
    merge: bool (default = False)
        boolean specifying whether niches neighborhoods should be merged according to logFCs
    min_perc: float (default = 0.1)
        minimum proportion for labeling niche neighborhoods
    feature_key: str (default = 'spatial_nhood')
        string specifying where spatial enrichment analysis should be stored in mdata object
    sig_key: str (default = 'SpatialFDR')
        string in mdata specifying significance scores (either PValue or SpatialFDR)
    sig_threshold: float (default = 0.05)
        threshold for significance
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
    try:
        adata.obs[fov_key] = adata.obs[fov_key].astype('str').astype('category')
    except:
        pass
        
    logging.info('computing spatial neighbors')
    if khop is not None:
        niche_df, _ = qu.tl.spatial_niches_khop(adata, radius = radius, p = 2, k = n_neighbors, khop = khop, min_cell_threshold = 0, labels_key = labels_key, spatial_key = spatial_key, fov_key = fov_key, n_jobs = n_jobs)
        adata_niche = anndata.AnnData(niche_df)
        adata_niche.obs = adata.obs.loc[niche_df.index, :]
    else:
        adata = qu.tl.compute_spatial_neighbors(adata, radius = radius, n_neighbors = n_neighbors, spatial_key = spatial_key, delaunay = delaunay, fov_key = fov_key, coord_type = coord_type)
        adata_niche, cells_nonn = qu.tl.compute_niche_composition(adata, labels_key = labels_key, min_cells = min_cells)

    adata_niche = adata_niche[np.where(pd.DataFrame(adata_niche.X).sum(1) != 0)[0], :].copy()
    adata = adata[np.where(pd.DataFrame(adata_niche.X).sum(1) != 0)[0], :].copy()
    if sketch_size is None:
        logging.info('skipping distribution-focused downsampling')
        adata_niche_subsample = adata_niche.copy()
    else:
        logging.info('performing distribution-focused downsampling')
        _, adata_niche_subsample = sketch(adata_niche, sample_set_key = sketch_key, gamma = gamma, num_subsamples = sketch_size, frequency_seed = frequency_seed, n_jobs = n_jobs)
    logging.info('computing between-patient niche similarity')
    adata_niche_subsample = qu.tl.construct_niche_similarity_graph(adata_niche_subsample, k = k_sim, n_jobs = n_jobs)
    logging.info('testing for differential spatial enrichment across conditions')
    mdata = quicheDA(adata_niche_subsample, design = design, model_contrasts = model_contrasts, patient_key = test_key)
    if label_scheme == 'fov_norm':
        condition = design.strip('~').split('+')[-1]
        if adata.obs[condition].nunique() <= 2:
            cell_type_abundance = adata_niche_subsample.obs.groupby([condition, labels_key]).size().unstack()
            norm_freq = cell_type_abundance.div(cell_type_abundance.sum(1), axis = 0)
            df = mdata['spatial_nhood'].to_df()
            df.index = mdata['spatial_nhood'].obs[condition]
        else:
            cell_type_abundance = pd.DataFrame(adata_niche_subsample.obs.groupby([labels_key]).size()).transpose()
            norm_freq = cell_type_abundance.div(cell_type_abundance.sum(1), axis = 0)     
            df = mdata['spatial_nhood'].to_df()
            df.index = [0]*len(df.index)
        annotations = compute_niche_abundance_fov_norm(df, nlargest = nlargest, min_perc = min_perc, norm_freq = norm_freq)
    elif label_scheme == 'normal':
        annotations = label_niches(mdata, nlargest = nlargest, min_perc = min_perc)
    elif label_scheme == 'fov':
        annotations = compute_niche_abundance_fov(mdata[feature_key].to_df(), nlargest = nlargest, min_perc = min_perc)
    elif label_scheme == 'neighborhood_norm':
        condition = design.strip('~').split('+')[-1]
        if adata.obs[condition].nunique() <= 2:
            cell_type_abundance = adata_niche_subsample.obs.groupby([condition, labels_key]).size().unstack()
            norm_freq = cell_type_abundance.div(cell_type_abundance.sum(1), axis = 0)      
        else:
            cell_type_abundance =  pd.DataFrame(adata_niche_subsample.obs.groupby([labels_key]).size()).transpose()
            norm_freq = cell_type_abundance.div(cell_type_abundance.sum(1), axis = 0)              
        annotations = compute_niche_abundance_neighborhood_norm(mdata, feature_key = feature_key,  anno_col = labels_key, nlargest = nlargest, min_perc = min_perc, norm_freq = norm_freq, condition = condition)
    else:
        annotations = compute_niche_abundance_neighborhood(mdata, feature_key = feature_key,  anno_col = labels_key, nlargest = nlargest, min_perc = min_perc)
    try:
        mdata['milo'].var[annotation_key] = annotations.values
    except:
        mdata['milo'].var[annotation_key] = annotations
    try:
        mdata['milo'].var[annotation_key].loc[np.isin(mdata['milo'].var['index_cell'], cells_nonn)] = 'unidentified'
    except:
        pass

    if merge == True:
        mdata, _ = label_niches_aggregate(mdata, annotation_key = annotation_key, aggregate_key =annotation_key, lfc_pct = 1)

    mdata = MuData({'expression': adata, feature_key: mdata[feature_key], 'quiche': mdata['milo']})
    mdata['quiche'].var[mdata['spatial_nhood'].obs.columns] = mdata['spatial_nhood'].obs.values
    mdata[feature_key].obs[annotation_key] = mdata['quiche'].var[annotation_key].values
    scores_df = pd.DataFrame(mdata['quiche'].var.groupby(annotation_key)[sig_key].median())
    scores_df.columns = ['pval']
    scores_df['logFC'] = mdata['quiche'].var.groupby(annotation_key)['logFC'].mean()
    scores_df = scores_df[scores_df['pval'] <= sig_threshold]
    scores_df = scores_df.iloc[np.where((scores_df['logFC'] <= -1) | (scores_df['logFC'] >= 1))[0]]
    ids = list(set(scores_df.index).intersection(set(list(mdata['quiche'].var[annotation_key].value_counts()[mdata['quiche'].var[annotation_key].value_counts() >= min_cells].index))))
    scores_df = scores_df.loc[ids]
    return mdata, scores_df

def quicheDA(adata,
            design = '~condition',
            model_contrasts = 'conditionA-conditionB',
            patient_key = 'Patient_ID',
            solver = 'edger', 
            feature_key = 'spatial_nhood'):
    """Performs condition-specific testing using QUICHE
    
    Parameters
    adata: anndata.AnnData (default = None)
        annotated data object containing preprocessed single-cell data 
    design: str (default = '~condition')
        string referring to design
    model_contrasts: str (default = 'conditionA-conditionB')
        string for condition-specific testing
    patient_key: str (default = 'Patient_ID')
        string referring to column in adata.obs with patient-level information
    solver: str (default = 'edger')
        string specifying solver
    feature_key: str (default = 'spatial_nhood')
        string specifying where spatial enrichment analysis should be stored in mdata object
    """
    milo = pt.tl.Milo()
    mdata = milo.load(adata, feature_key = feature_key)
    mdata[feature_key].uns["nhood_neighbors_key"] = None
    mdata = qu.tl.build_milo_graph(mdata, feature_key = feature_key)
    mdata = milo.count_nhoods(mdata, sample_col = patient_key, feature_key = feature_key)
    milo.da_nhoods(mdata, design = design, model_contrasts = model_contrasts, feature_key = feature_key, solver = solver)
    return mdata

def percent_change(before, after):
    perc_change = ((after - before)/ before)*100
    return perc_change
    
def label_niches_aggregate(mdata, key = 'milo', annotation_key = 'original_annotations', aggregate_key = 'mapped_annotations', lfc_pct = 0.5, logFC_cutoff = 0.85):
    second_round_annotations = {}
    for group in np.unique(mdata[key].var[annotation_key]):
        try:
            comparison_groups = np.unique([i for i in np.unique(mdata[key].var[annotation_key]) if np.isin(i.split('__'), group.split('__')).all()])
            idx_group = np.where(mdata[key].var[annotation_key] == group)[0]
            logfc_group = mdata[key].var.iloc[idx_group, :].loc[:, 'logFC'].mean()
            logfc_comparisons_before = []
            logfc_comparisons_after = []
            for comparison_group in comparison_groups:
                idx_comparison_before = np.where(mdata[key].var[annotation_key] == comparison_group)[0]
                logfc_comparisons_before.append(mdata[key].var.iloc[idx_comparison_before, :].loc[:, 'logFC'].mean())

                tmp_mapping = pd.DataFrame(mdata[key].var[annotation_key]).copy()
                tmp_mapping[np.isin(tmp_mapping[annotation_key], comparison_group)] = group
                idx_comparison_after = np.where(tmp_mapping == group)[0]
                logfc_comparisons_after.append(mdata[key].var.iloc[idx_comparison_after, :].loc[:, 'logFC'].mean())

            log_fc_group_list = np.repeat(logfc_group, len(logfc_comparisons_after))
            
            perc_change = [percent_change(log_fc_group_list[i], logfc_comparisons_after[i]) for i in range(0, len(log_fc_group_list))]
            idx_group_comparison = np.where(comparison_groups != group)[0]
            idx_group_original = np.where(comparison_groups == group)[0]
            if (np.array(logfc_comparisons_after).max() > -logFC_cutoff) and (np.array(logfc_comparisons_after).max() < logFC_cutoff): #if niche becomes not significant then keep original
                second_round_annotations[group] = group
                second_round_annotations[group] = comparison_groups[np.where(np.array(perc_change) == np.array(perc_change).max())[0]][0]
            else:
                print(group, np.array(logfc_comparisons_after).max(), comparison_groups, perc_change, np.abs(np.array(perc_change).max()))
                second_round_annotations[group] = group
        except:
            second_round_annotations[group] = group #worse then keep   
    mdata[key].var[aggregate_key] = pd.Series(mdata[key].var[annotation_key]).map(second_round_annotations)
    return mdata, second_round_annotations

def top_labels_with_condition(row, nlargest, min_perc = 0.1):
    top_labels = row.nlargest(nlargest)
    top_labels = top_labels.index[top_labels > min_perc]
    sorted_labels = '__'.join(sorted(top_labels))
    return sorted_labels

def label_niches(mdata, nlargest = 3, min_perc = 0.1):
    knn_mat = mdata['spatial_nhood'].obsp['connectivities']
    df_prop = mdata['spatial_nhood'].to_df()
    annotations = []
    for cell_idx in range(knn_mat.shape[0]):
        neighbor_indices = knn_mat[cell_idx].nonzero()[1]
        avg_abundances = df_prop.iloc[neighbor_indices].mean(axis=0)
        label = top_labels_with_condition(avg_abundances, nlargest=nlargest, min_perc=min_perc)
        annotations.append(label)
    return annotations

def compute_niche_abundance_neighborhood(mdata, feature_key = 'spatial_nhood', nlargest = 3, anno_col = 'cell_cluster', min_perc = 0.1):
    anno_dummies = pd.get_dummies(mdata[feature_key].obs[anno_col])
    anno_count = mdata[feature_key].obsm["nhoods"].T.dot(csr_matrix(anno_dummies.values))
    anno_count_dense = anno_count.toarray()
    anno_sum = anno_count_dense.sum(1)
    anno_frac = np.divide(anno_count_dense, anno_sum[:, np.newaxis])
    anno_frac_dataframe = pd.DataFrame(anno_frac, columns=anno_dummies.columns, index=mdata["milo"].var_names)
    annotations = anno_frac_dataframe.apply(top_labels_with_condition, nlargest = nlargest, min_perc = min_perc, axis=1)
    return annotations

def compute_niche_abundance_neighborhood_norm(mdata, feature_key = 'spatial_nhood', nlargest = 3, anno_col = 'cell_cluster', min_perc = 0.1, norm_freq = None, condition = 'condition'):
    anno_dummies = pd.get_dummies(mdata[feature_key].obs[anno_col])
    anno_count = mdata[feature_key].obsm["nhoods"].T.dot(csr_matrix(anno_dummies.values))
    anno_count_dense = anno_count.toarray()
    anno_sum = anno_count_dense.sum(1)
    anno_frac = np.divide(anno_count_dense, anno_sum[:, np.newaxis])
    anno_frac_dataframe = pd.DataFrame(anno_frac, columns=anno_dummies.columns, index=mdata["milo"].var_names)
    if norm_freq.shape[0] == 1:
        anno_frac_dataframe.index = [0]*len(anno_frac_dataframe.index)
    else:
        anno_frac_dataframe.index = mdata[feature_key].obs[condition].values
    annotations = anno_frac_dataframe.apply(top_labels_with_condition_norm, nlargest = nlargest, min_perc = min_perc, axis=1, norm_freq = norm_freq)
    return annotations

def compute_niche_abundance_fov(df, nlargest = 3, min_perc = 0.1):
    annotations = df.apply(top_labels_with_condition, nlargest = nlargest, min_perc = min_perc, axis=1)
    return annotations

def compute_niche_abundance_fov_norm(df, nlargest = 3, min_perc = 0.1, norm_freq = None):
    annotations = df.apply(top_labels_with_condition_norm, nlargest = nlargest, min_perc = min_perc, axis=1, norm_freq = norm_freq)
    return annotations

def top_labels_with_condition_norm(row, nlargest, min_perc = 0.1, norm_freq = None):
    top_labels = row.div(norm_freq.loc[row.name]).nlargest(nlargest)
    top_labels = top_labels.index[top_labels > min_perc]
    sorted_labels = '__'.join(sorted(top_labels))
    return sorted_labels

def compute_functional_expression(mdata = None,
                                sig_niches = None,
                                labels_key = 'cell_cluster',
                                annot_key = 'quiche_niche',
                                fov_key = 'fov',
                                segmentation_label_key = 'label',
                                patient_key = 'Patient_ID',
                                min_cell_count = 3,
                                foldchange_key = 'logFC',
                                markers = None):
    """Computes the expression of cell types within outcome-associated niches
    
    Parameters
    mdata: mudata (default = None)
        annotated data object containing QUICHE enrichment analysis
    sig_niches: list (default = None)
        list containing niches of interest
    labels_key: str (default = 'cell_cluster')
        string in adata.obs containing cell phenotype information
    fov_key: str (default = 'fov')
        string referring to column in adata.obs with sample (fov) information
    annot_key: str (default = 'quiche_niche')
        string specifying column in mdata with predicted niches
    segmentation_label_key: str (default = 'label')
        string specifying column in mdata with cell segmentation cell labels
    patient_key: str (default = 'Patient_ID')
        string specifying column in mdata with patient-level information
    min_cell_count: int (default = 3)
        integer specifying the minimum number of nearest neighbors within a niche to compute expression analysis
    foldchange_key: str (default = 'logFC')
        string specifying column in mdata with predicted logFCs
    markers: list (default = None)
        list of functional markers
    """
    #subset connectivities matrix according to downsampled cells
    idx = mdata['expression'].obs_names.get_indexer(mdata['spatial_nhood'].obs_names)
    conn_mat = mdata['expression'].obsp['spatial_connectivities'][idx, :] #downsampled cells x all cells
    #subset data according to significant niches of interest
    sig_bool = np.isin(mdata['quiche'].var[annot_key], sig_niches)
    conn_mat = conn_mat[sig_bool, :]
    niche_list = mdata['quiche'].var[annot_key][sig_bool]
    nn_array = [row.nonzero()[1] for row in conn_mat] #nn around cell, where elements == 1
    #create a dictionary mapping cell types to indices
    cell_clusters_indices = {cell_type: np.where(mdata['expression'].obs[labels_key] == cell_type)[0] for cell_type in mdata['expression'].obs[labels_key].unique()}
    #return the functional marker expression of cell types of interest within each significant niche in the original undownsampled FOV
    func_df = []
    for i in range(0, len(nn_array)):
        niche = niche_list[i]
        nn = nn_array[i]
        for cell_type in niche.split('__'):
            idx_cell_type = cell_clusters_indices[cell_type]
            idx_cell_type_nn = list(set(nn).intersection(set(idx_cell_type)))
            if len(idx_cell_type_nn) >= min_cell_count:
                exp = mdata['expression'][idx_cell_type_nn, :][:, markers].to_df()
                exp[annot_key] = niche
                exp[labels_key] = cell_type
                exp[segmentation_label_key] = mdata['expression'].obs[segmentation_label_key]
                exp['quiche_niche_cell_type'] = niche + ':' + cell_type
                exp[fov_key] = mdata['expression'][idx_cell_type_nn].obs[fov_key].unique()[0]
                exp[patient_key] = mdata['expression'][idx_cell_type_nn].obs[patient_key].unique()[0]
                func_df.append(exp)
    func_df = pd.concat(func_df, axis = 0)
    adata_func = anndata.AnnData(func_df.drop(columns = [annot_key, labels_key, segmentation_label_key, 'quiche_niche_cell_type', fov_key, patient_key]))
    adata_func.obs = func_df.loc[:, [annot_key, labels_key, 'quiche_niche_cell_type', segmentation_label_key, fov_key, patient_key]]
    adata_func.obs = pd.merge(adata_func.obs, pd.DataFrame(mdata['quiche'].var.groupby([annot_key])[foldchange_key].mean()), on = [annot_key])
    return adata_func

def relabel_niches_celltype(adata, annotation_key, cluster_key, sig_niches):
    adata.obs['niche_cell_type'] = pd.DataFrame(list(adata.obs[annotation_key].values)) + ' niche: ' + pd.DataFrame(list(adata.obs[cluster_key].values))
    niche_dict = dict(zip(list(adata.obs[annotation_key].values), list(adata.obs['niche_cell_type'].values)))
    groups = [i.split(':')[0] for i in pd.Series(sig_niches).map(niche_dict)]
    mapped_niches = list(dict.fromkeys([i for g in groups for i in np.unique(adata.obs['niche_cell_type']) if g in i]).keys())
    sorted_mapped_niches = [item for group in groups for item in mapped_niches if item.startswith(group)]
    return sorted_mapped_niches

def create_single_positive_table(marker_vals, threshold_list):
    """ Determine whether a cell is positive for a marker based on the provided threshold.
    Args:
        marker_vals (pd.DataFrame): dataframe containing the marker intensity values
        threshold_list (list): list of functional markers and their pre-determined thresholds

    Returns:
        pd.DataFrame:
            contains the marker intensities as well as the single positive marker data
    """
    # create binary functional marker table, append to anndata table
    for marker, threshold in threshold_list:
        marker_vals[marker] = (marker_vals[marker].values >= threshold).astype('int')

    return marker_vals

def compute_patient_proportion(mdata,
                                niches = None,
                                feature_key = 'quiche',
                                annot_key = 'quiche_niche',
                                patient_key = 'Patient_ID',
                                design_key = 'Relapse',
                                patient_niche_threshold = 3):
    """Computes niche level metadata
    
    Parameters
    mdata: mudata (default = None)
        annotated data object containing QUICHE enrichment analysis
    niches: list (default = None)
        list containing niches of interest
    design_key: str (default = 'Relapse')
        string specifying column in mdata object that condition-specific testing was performed on
    patient_key: str (default = 'Patient_ID')
        string referring to column in adata.obs with patient-level information
    annot_key: str (default = 'quiche_niche')
        string referring to the column with niche annotations
    feature_key: str (default = 'quiche')
        string specifying where spatial enrichment analysis is stored
    patient_niche_threshold: int (default = 3)
        integer referring to the minimum number of niches per patient samples to be considered +
    """
    count_df = mdata[feature_key].var.groupby([annot_key, patient_key, design_key]).size()[mdata[feature_key].var.groupby([annot_key, patient_key, design_key]).size() > patient_niche_threshold].reset_index()
    avg_counts = count_df.groupby([annot_key, design_key])[0].mean().reset_index()
    avg_counts.columns = [annot_key, design_key, 'avg_niche_abundance']
    patient_count_df = pd.DataFrame(count_df.groupby(annot_key)[patient_key].nunique())
    niche_list = list(set(niches).intersection(set(patient_count_df.index)))
    patient_count_df = patient_count_df.loc[niche_list].sort_values(by = patient_key)

    cov_count_df = pd.DataFrame(count_df.groupby([annot_key, design_key])[[patient_key]].size())
    cov_count_df = cov_count_df.loc[niche_list].sort_values(by = design_key)
    cov_count_df = cov_count_df.reset_index()
    cov_count_df.columns = [annot_key, design_key, 'patient_count']
    
    patient_ids_df = count_df.groupby([annot_key, design_key])[patient_key].apply(lambda x: list(x.unique())).reset_index()
    patient_ids_df.columns = [annot_key, design_key, 'patient_ids']

    cov_count_df = pd.merge(cov_count_df, patient_ids_df, on=[annot_key, design_key])

    cov_count_df['med_logFC'] = mdata[feature_key].var.groupby(annot_key)['logFC'].median().loc[cov_count_df[annot_key]].values
    cov_count_df['mean_logFC'] = mdata[feature_key].var.groupby(annot_key)['logFC'].mean().loc[cov_count_df[annot_key]].values
    cov_count_df['pval'] = mdata[feature_key].var.groupby(annot_key)['SpatialFDR'].median().loc[cov_count_df[annot_key]].values
    
    cov_count_df = pd.merge(cov_count_df, avg_counts, on=[annot_key, design_key])

    cov_total = mdata[feature_key].var.groupby([design_key])[patient_key].nunique()
    cov_total = pd.DataFrame(cov_total)
    cov_total.columns = ['patient_cov']
    cov_count_df = pd.merge(cov_count_df, cov_total, on=[design_key])
    cov_count_df['prop_cov'] = cov_count_df['patient_count'].div(cov_count_df['patient_cov'])
    return cov_count_df