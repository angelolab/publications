import scanpy as sc
import pertpy as pt
import pandas as pd
import numpy as np
import os
from scipy.sparse import csr_matrix
from abc import abstractmethod
from sklearn.preprocessing import OneHotEncoder, LabelEncoder, label_binarize, StandardScaler
from sklearn.model_selection import StratifiedShuffleSplit, StratifiedKFold, GridSearchCV, ShuffleSplit
from sklearn.metrics import f1_score, roc_auc_score, accuracy_score, precision_score, balanced_accuracy_score, confusion_matrix
import scipy
from rpy2.robjects import pandas2ri
from scipy.sparse import csr_matrix
import anndata
import quiche as qu
from sketchKH import sketch
import logging
from muon import MuData
import pertpy as pt 
import scanpy as sc
logging.basicConfig(level=logging.INFO)
pandas2ri.activate()

def miloDA(adata,
            n_neighbors = 100,
            design = '~condition',
            model_contrasts = 'conditionA-conditionB',
            patient_key = 'Patient_ID',
            cluster_key = 'mask_name',
            feature_key = 'expression',
            sig_threshold = 0.05,
            gamma = 1,
            sketch_size = 1000,
            frequency_seed = 0,
            n_jobs = -1,
            **kwargs):
    _, adata_subsample = sketch(adata, sample_set_key = patient_key, gamma = gamma, num_subsamples = sketch_size, frequency_seed = frequency_seed, n_jobs = n_jobs)
    milo = pt.tl.Milo()
    mdata = milo.load(adata_subsample, feature_key = feature_key)
    sc.tl.pca(mdata[feature_key])
    sc.pp.neighbors(mdata[feature_key], n_neighbors = n_neighbors, use_rep = 'X_pca')
    mdata[feature_key].uns["nhood_neighbors_key"] = None
    mdata = qu.tl.build_milo_graph(mdata, feature_key = feature_key)
    mdata = milo.count_nhoods(mdata, sample_col = patient_key, feature_key = feature_key)
    milo.da_nhoods(mdata, design = design, model_contrasts = model_contrasts, feature_key = feature_key)
    milo.annotate_nhoods(mdata, anno_col = cluster_key, feature_key = feature_key)
    mdata = MuData({'expression': mdata[feature_key], 'milo': mdata['milo']})
    mdata['milo'].var[patient_key] = mdata['expression'].obs[patient_key].values
    mdata['milo'].var[design.split('~')[1]] = mdata['expression'].obs[design.split('~')[1]].values
    scores_df = pd.DataFrame(mdata['milo'].var.groupby('nhood_annotation')['SpatialFDR'].median()[mdata['milo'].var.groupby('nhood_annotation')['SpatialFDR'].median() < sig_threshold])
    scores_df.columns = ['pval']
    return mdata, scores_df

def quicheDA(adata,
                design = '~condition',
                model_contrasts = 'conditionA-conditionB',
                patient_key = 'Patient_ID',
                feature_key = 'spatial_nhood'):
    milo = pt.tl.Milo()
    mdata = milo.load(adata, feature_key = feature_key)
    mdata[feature_key].uns["nhood_neighbors_key"] = None
    mdata = qu.tl.build_milo_graph(mdata, feature_key = feature_key)
    mdata = milo.count_nhoods(mdata, sample_col = patient_key, feature_key = feature_key)
    milo.da_nhoods(mdata, design = design, model_contrasts = model_contrasts, feature_key = feature_key)
    return mdata

def run_quiche(adata,
               radius = 200,
               labels_key = 'mask_name',
               spatial_key = 'spatial',
               fov_key = 'Patient_ID',
               patient_key = 'Patient_ID',
               n_neighbors = None,
               delaunay = True,
               min_cells = 3,
               k_sim = 100,
               design = '~condition',
               model_contrasts = 'conditioncancer_core-conditioncancer_border',
               sketch_size = None,
               frequency_seed = 0,
               sig_threshold = 0.05,
               nlargest = 4,
               min_perc = 0.1,
               gamma = 1,
               feature_key = 'spatial_nhood',
               annotation_key = 'quiche_niche',
               coord_type = 'generic',
               n_jobs = -1,
               **kwargs):
        
    logging.info('computing spatial neighbors')
    adata = qu.tl.compute_spatial_neighbors(adata, radius = radius, n_neighbors = n_neighbors, spatial_key = spatial_key, delaunay = delaunay, fov_key = fov_key, coord_type = coord_type)
    adata_niche, cells_nonn = qu.tl.compute_niche_composition(adata, labels_key = labels_key, min_cells = min_cells)
    adata_niche = adata_niche[np.where(pd.DataFrame(adata_niche.X).sum(1) != 0)[0], :].copy()
    if sketch_size is None:
        logging.info('skipping distribution-focused downsampling')
        adata_niche_subsample = adata_niche.copy()
    else:
        logging.info('performing distribution-focused downsampling')
        _, adata_niche_subsample = sketch(adata_niche, sample_set_key = patient_key, gamma = gamma, num_subsamples = sketch_size, frequency_seed = frequency_seed, n_jobs = n_jobs)
    logging.info('computing between-patient niche similarity')
    adata_niche_subsample = qu.tl.construct_niche_similarity_graph(adata_niche_subsample, k = k_sim, n_jobs = n_jobs)
    logging.info('testing for differential spatial enrichment across conditions')
    mdata = quicheDA(adata_niche_subsample, design = design, model_contrasts = model_contrasts, patient_key = patient_key)
    annotations = compute_niche_abundance(mdata[feature_key].to_df(), nlargest = nlargest, min_perc = min_perc)
    
    mdata['milo'].var[annotation_key] = annotations.values
    mdata['milo'].var[annotation_key].loc[np.isin(mdata['milo'].var['index_cell'], cells_nonn)] = 'unidentified'

    mdata = MuData({'expression': adata, feature_key: mdata[feature_key], 'quiche': mdata['milo']})
    mdata[feature_key].obs[annotation_key] = mdata['quiche'].var[annotation_key].values
    mdata['quiche'].var[design.split('~')[1]] = mdata[feature_key].obs[design.split('~')[1]].values
    mdata['quiche'].var[patient_key] = mdata[feature_key].obs[patient_key].values

    scores_df = pd.DataFrame(mdata['quiche'].var.groupby(annotation_key)['SpatialFDR'].median(), columns = ['pval'])
    scores_df['logFC'] = mdata['quiche'].var.groupby(annotation_key)['logFC'].median()
    scores_df = scores_df[scores_df['pval'] < sig_threshold]

    scores_df = pd.DataFrame(mdata['quiche'].var.groupby(annotation_key)['SpatialFDR'].median())
    scores_df.columns = ['pval']
    scores_df['logFC'] = mdata['quiche'].var.groupby(annotation_key)['logFC'].median()
    scores_df = scores_df[scores_df['pval'] < sig_threshold]
    ids = list(set(scores_df.index).intersection(set(list(mdata['quiche'].var[annotation_key].value_counts()[mdata['quiche'].var[annotation_key].value_counts() >= min_cells].index))))
    scores_df = scores_df.loc[ids]

    return mdata, scores_df

# def annotate_nhoods(
#     mdata,
#     anno_col: str,
#     feature_key = "rna",
#     nlargest = 2
# ):
#     """Assigns a categorical label to neighbourhoods, based on the most frequent label among cells in each neighbourhood. This can be useful to stratify DA testing results by cell types or samples.

#     Args:
#         mdata: MuData object
#         anno_col: Column in adata.obs containing the cell annotations to use for nhood labelling
#         feature_key: If input data is MuData, specify key to cell-level AnnData object. Defaults to 'rna'.

#     Returns:
#         None. Adds in place:
#         - `milo_mdata['milo'].var["nhood_annotation"]`: assigning a label to each nhood
#         - `milo_mdata['milo'].var["nhood_annotation_frac"]` stores the fraciton of cells in the neighbourhood with the assigned label
#         - `milo_mdata['milo'].varm['frac_annotation']`: stores the fraction of cells from each label in each nhood
#         - `milo_mdata['milo'].uns["annotation_labels"]`: stores the column names for `milo_mdata['milo'].varm['frac_annotation']`

#     """
#     try:
#         sample_adata = mdata["milo"]
#     except KeyError:
#         print(
#             "milo_mdata should be a MuData object with two slots: feature_key and 'milo' - please run milopy.count_nhoods(adata) first"
#         )
#         raise
#     adata_run = mdata[feature_key].copy()

#     # Check value is not numeric
#     if pd.api.types.is_numeric_dtype(adata_run.obs[anno_col]):
#         raise ValueError(
#             "adata.obs[anno_col] is not of categorical type - please use milopy.utils.annotate_nhoods_continuous for continuous variables"
#         )

#     anno_dummies = pd.get_dummies(adata_run.obs[anno_col])
#     anno_count = adata_run.obsm["nhoods"].T.dot(csr_matrix(anno_dummies.values))
#     anno_count_dense = anno_count.toarray()
#     anno_sum = anno_count_dense.sum(1)
#     anno_frac = np.divide(anno_count_dense, anno_sum[:, np.newaxis])
#     anno_frac_dataframe = pd.DataFrame(anno_frac, columns=anno_dummies.columns, index=sample_adata.var_names)

#     sample_adata.varm["frac_annotation"] = anno_frac_dataframe.values
#     sample_adata.uns["annotation_labels"] = anno_frac_dataframe.columns
#     sample_adata.uns["annotation_obs"] = anno_col
#     sample_adata.var["nhood_annotation"] = anno_frac_dataframe.apply(top_labels_with_condition, nlargest = nlargest, axis=1) #anno_frac_dataframe.idxmax(1)
#     sample_adata.var["nhood_annotation_frac"] = anno_frac_dataframe.max(1)

# def top_labels_with_condition(row, nlargest):
#     top_labels = row.nlargest(nlargest)
#     top_labels = top_labels.index[top_labels > 0.02]
#     sorted_labels = '__'.join(sorted(top_labels))
#     return sorted_labels

# class BaseLabelPropagation:
#     """Class for performing label propagation

#     Parameters
#     W: ndarray
#         adjacency matrix to compute label propagation on
#     ----------

#     Returns
#     ----------
#     """
#     def __init__(self, W):
#         self.W_norm = self._normalize(W)
#         self.n_nodes = np.shape(W)[0]
#         self.indicator_labels = None
#         self.n_classes = None
#         self.labeled_mask = None
#         self.predictions = None

#     @staticmethod
#     @abstractmethod
#     def _normalize(W):
#         raise NotImplementedError("_normalize must be implemented")

#     @abstractmethod
#     def _propagate(self):
#         raise NotImplementedError("_propagate must be implemented")


#     def _encode(self, labels):
#         # Get the number of classes
#         classes = np.unique(labels)
#         classes = classes[classes != -1] #-1 are unlabeled nodes so we'll exclude them
#         self.n_classes = np.shape(classes)[0]
#         # One-hot encode labeled data instances and zero rows corresponding to unlabeled instances
#         unlabeled_mask = (labels == -1)
#         labels = labels.copy()
#         labels[unlabeled_mask] = 0
#         onehot_encoder = OneHotEncoder(sparse_output=False)
#         self.indicator_labels = labels.reshape(len(labels), 1)
#         self.indicator_labels = onehot_encoder.fit_transform(self.indicator_labels)
#         self.indicator_labels[unlabeled_mask, 0] = 0

#         self.labeled_mask = ~unlabeled_mask

#     def fit(self, labels, max_iter, tol):
#         """Fits semisupervised label propagation model

#         Parameters
#         labels: ndarray
#             labels for every node, where -1 indicates unlabeled nodes
#         max_iter: int (default = 10000)
#             maximum number of iterations before stopping prediction
#         tol: float (default = 1e-3)
#             float referring to the error tolerance between runs. If unchanging, stop prediction
#         """
#         self._encode(labels)

#         self.predictions = self.indicator_labels.copy()
#         prev_predictions = np.zeros((self.n_nodes, self.n_classes), dtype = np.float)

#         for i in range(max_iter):
#             # Stop iterations if the system is considered at a steady state
#             variation = np.abs(self.predictions - prev_predictions).sum().item()

#             if variation < tol:
#                 print(f"The method stopped after {i} iterations, variation={variation:.4f}.")
#                 break

#             prev_predictions = self.predictions
#             self._propagate()

#     def predict(self):
#         return self.predictions

#     def predict_labels(self):
#         """
#         Returns
#         predicted_labels: ndarray
#             array of predicted labels according to the maximum probability
#         predicted_scores: ndarray
#             array of probability scores with dimensions n x nclasses
#         uncertainty: ndarray
#             array 1 - max of predictions

#         ----------
#         """
#         predicted_labels = np.argmax(self.predictions, axis = 1)
#         predicted_scores = self.predictions
#         uncertainty = 1 - np.max(predicted_scores, 1)

#         return predicted_labels, predicted_scores, uncertainty

# class LabelPropagation(BaseLabelPropagation):
#     def __init__(self, W):
#         super().__init__(W)

#     @staticmethod
#     def _normalize(W):
#         """ Computes row normalized adjacency matrix: D^-1 * W"""
#         d = W.sum(axis=0).getA1()
#         d = 1/d
#         D = scipy.sparse.diags(d)

#         return D @ W

#     def _propagate(self):
#         self.predictions = self.W_norm @ self.predictions

#         # Put back already known labels
#         self.predictions[self.labeled_mask] = self.indicator_labels[self.labeled_mask]

#     def fit(self, labels, max_iter = 500, tol = 1e-3):
#         super().fit(labels, max_iter, tol)

# def annotate_nhoods_label_prop(
#     mdata,
#     n_splits = 20,
#     anno_col = None,
#     feature_key = "rna",
#     train_size = 0.6,
#     nlargest = 2,
#     label_key = 'nhood_annotation', 
# ):
#     """Assigns a categorical label to neighbourhoods, based on the most frequent label among cells in each neighbourhood. This can be useful to stratify DA testing results by cell types or samples.

#     Args:
#         mdata: MuData object
#         anno_col: Column in adata.obs containing the cell annotations to use for nhood labelling
#         feature_key: If input data is MuData, specify key to cell-level AnnData object. Defaults to 'rna'.

#     Returns:
#         None. Adds in place:
#         - `milo_mdata['milo'].var["nhood_annotation"]`: assigning a label to each nhood
#         - `milo_mdata['milo'].var["nhood_annotation_frac"]` stores the fraciton of cells in the neighbourhood with the assigned label
#         - `milo_mdata['milo'].varm['frac_annotation']`: stores the fraction of cells from each label in each nhood
#         - `milo_mdata['milo'].uns["annotation_labels"]`: stores the column names for `milo_mdata['milo'].varm['frac_annotation']`

#     """
#     try:
#         sample_adata = mdata["milo"]
#     except KeyError:
#         print(
#             "milo_mdata should be a MuData object with two slots: feature_key and 'milo' - please run milopy.count_nhoods(adata) first"
#         )
#         raise
#     adata_run = mdata[feature_key].copy()

#     # Check value is not numeric
#     if pd.api.types.is_numeric_dtype(adata_run.obs[anno_col]):
#         raise ValueError(
#             "adata.obs[anno_col] is not of categorical type - please use milopy.utils.annotate_nhoods_continuous for continuous variables"
#         )

#     anno_dummies = pd.get_dummies(adata_run.obs[anno_col])
#     # anno_count = adata_run.obsm["nhoods"].T.dot(csr_matrix(anno_dummies.values))
#     # anno_count_dense = anno_count.toarray()
#     # anno_sum = anno_count_dense.sum(1)
#     # anno_frac = np.divide(anno_count_dense, anno_sum[:, np.newaxis])
#     anno_frac = adata_run.X.copy()
#     anno_frac_dataframe = pd.DataFrame(anno_frac, columns=anno_dummies.columns, index=sample_adata.var_names)

#     sample_adata.varm["frac_annotation"] = anno_frac_dataframe.values
#     sample_adata.uns["annotation_labels"] = anno_frac_dataframe.columns
#     sample_adata.uns["annotation_obs"] = anno_col
#     og_df = pd.DataFrame(mdata['rna'].X, columns = mdata['rna'].var_names)
#     y_max = og_df.apply(top_labels_with_condition, nlargest = nlargest, axis=1) #anno_frac_dataframe.idxmax(1)
#     y_max[np.isin(list(y_max.values), list(y_max.value_counts()[y_max.value_counts() < n_splits].index))] = 'unidentified'
#     le = LabelEncoder()
#     y = le.fit_transform(y_max).astype(int)
#     # y.loc[y.index[y.value_counts() < 10]] = -1. #if barely any cells in this category, we say unconfident annotation and learn their labels based on neighbors
#     n_nodes = np.shape(y)[0]
#     X = mdata['rna'].X
#     # X = anno_frac_dataframe.values
#     sss = StratifiedShuffleSplit(n_splits = n_splits, train_size = train_size, random_state = 0)
#     sss.get_n_splits(X, y = y)
#     i=0
#     predicted_scores_total = 0
#     # sc.pp.neighbors(adata_run, n_neighbors = 50)
#     for train_index, test_index in sss.split(X, y):
#         print(i)
#         i=i+1
#         y_t = np.full(n_nodes, -1.)
#         y_t[train_index] = y[train_index].copy()
#         label_propagation = LabelPropagation(mdata['rna'].obsp['connectivities'])
#         label_propagation.fit(y_t)
#         _, predicted_scores, uncertainty = label_propagation.predict_labels()
#         predicted_scores_total += predicted_scores
#     predicted_scores_avg = predicted_scores_total / n_splits
#     predicted_labels_avg = np.argmax(predicted_scores_avg, axis = 1)   
#     predicted_labels_avg = le.inverse_transform(predicted_labels_avg)

#     sample_adata.var[label_key] = predicted_labels_avg
#     sample_adata.var[f'{label_key}_frac'] = anno_frac_dataframe.max(1)
#     return mdata, predicted_scores_avg, uncertainty

# def top_labels_with_condition(row, nlargest):
#     top_labels = row.nlargest(nlargest)
#     top_labels = top_labels.index[top_labels > 0.10]
#     sorted_labels = '__'.join(sorted(top_labels))
#     return sorted_labels

# def label_niches_abundance(mdata, feature_key = 'rna', anno_col = 'cell_cluster', nlargest = 3, annotation_key = 'original_annotations'):
#     sample_adata = mdata["milo"]

#     adata_run = mdata[feature_key].copy()
#     anno_dummies = pd.get_dummies(adata_run.obs[anno_col])
#     anno_count = adata_run.obsm["nhoods"].T.dot(csr_matrix(anno_dummies.values))
#     anno_count_dense = anno_count.toarray()
#     anno_sum = anno_count_dense.sum(1)
#     anno_frac = np.divide(anno_count_dense, anno_sum[:, np.newaxis])
#     anno_frac_dataframe = pd.DataFrame(anno_frac, columns=anno_dummies.columns, index=sample_adata.var_names)
#     idx_empty = np.where(anno_frac_dataframe.sum(1) == 1)[0]
#     largest_group_annotation = anno_frac_dataframe.apply(top_labels_with_condition, nlargest = nlargest, axis=1) #anno_frac_dataframe.idxmax(1)
#     largest_group_annotation[idx_empty] = 'empty'
#     mdata['milo'].var[annotation_key] = largest_group_annotation.values
#     return mdata

# def percent_change_directional(before, after):
#     if (before < 0 and after < before) or (before > 0 and after > before):
#         perc_change = abs(((after - before) / abs(before)) * 100)
#     elif before == after:
#         perc_change = 0
#     else: 
#         perc_change = -abs(((after - before) / abs(before)) * 100)
#     return perc_change
    
# def label_niches_aggregate(mdata, annotation_key = 'original_annotations', aggregate_key = 'mapped_annotations', lfc_pct = 15):
#     second_round_annotations = {}
#     for group in np.unique(mdata['milo'].var[annotation_key]):
#         try:
#             comparison_groups = np.unique([i for i in np.unique(mdata['milo'].var[annotation_key]) if np.isin(i.split('__'), group.split('__')).all()])
#             idx_group = np.where(mdata['milo'].var[annotation_key] == group)[0]
#             logfc_group = mdata['milo'].var.iloc[idx_group, :].loc[:, 'logFC'].mean()
#             logfc_comparisons_before = []
#             logfc_comparisons_after = []
#             for comparison_group in comparison_groups:
#                 idx_comparison_before = np.where(mdata['milo'].var[annotation_key] == comparison_group)[0]
#                 logfc_comparisons_before.append(mdata['milo'].var.iloc[idx_comparison_before, :].loc[:, 'logFC'].mean())

#                 tmp_mapping = pd.DataFrame(mdata['milo'].var[annotation_key]).copy()
#                 tmp_mapping[np.isin(tmp_mapping[annotation_key], comparison_group)] = group
#                 idx_comparison_after = np.where(tmp_mapping == group)[0]
#                 logfc_comparisons_after.append(mdata['milo'].var.iloc[idx_comparison_after, :].loc[:, 'logFC'].mean())

#             log_fc_group_list = np.repeat(logfc_group, len(logfc_comparisons_after))
            
#             perc_change = [percent_change_directional(log_fc_group_list[i], logfc_comparisons_after[i]) for i in range(0, len(log_fc_group_list))]
#             idx_group_comparison = np.where(comparison_groups != group)[0]
#             if np.array(perc_change).max() >= lfc_pct:
#                 second_round_annotations[group] = comparison_groups[np.where(np.array(perc_change) == np.array(perc_change).max())[0]][0]
#             elif (abs(np.array(perc_change))[idx_group_comparison] < lfc_pct).any():
#                 second_round_annotations[group] = comparison_groups[np.where(np.abs(np.array(perc_change)) == np.abs(np.array(perc_change)[idx_group_comparison]).min())[0]][0]
#             else:
#                 second_round_annotations[group] = group
            
#             #compare logFC of group to all comparison groups 
#         except:
#             second_round_annotations[group] = group #worse then keep   
#     mdata['milo'].var[aggregate_key] = pd.Series(mdata['milo'].var[annotation_key]).map(second_round_annotations)
#     return mdata, second_round_annotations

def top_labels_with_condition(row, nlargest, min_perc = 0.1):
    top_labels = row.nlargest(nlargest)
    top_labels = top_labels.index[top_labels > min_perc]
    sorted_labels = '__'.join(sorted(top_labels))
    return sorted_labels

# def label_niches_abundance(mdata, feature_key = 'rna', anno_col = 'cell_cluster', nlargest = 3, annotation_key = 'original_annotations'):
#     adata_run = mdata[feature_key].copy()
#     largest_group_annotation = adata_run.to_df().apply(top_labels_with_condition, nlargest = nlargest, axis=1)
#     mdata['milo'].var[annotation_key] = largest_group_annotation.values
#     return mdata

def compute_niche_abundance(df, nlargest = 3, min_perc = 0.1):
    annotations = df.apply(top_labels_with_condition, nlargest = nlargest, min_perc = min_perc, axis=1)
    return annotations

# def compute_functional_expression(mdata = None,
#                                 sig_niches = None,
#                                 labels_key = 'cell_cluster',
#                                 annot_key = 'quiche_niche',
#                                 fov_key = 'fov',
#                                 patient_key = 'Patient_ID',
#                                 min_cell_count = 3,
#                                 markers = None):

#     #subset connectivities matrix according to downsampled cells
#     idx = mdata['expression'].obs_names.get_indexer(mdata['spatial_nhood'].obs_names)
#     conn_mat = mdata['expression'].obsp['spatial_connectivities'][idx, :] #downsampled cells x all cells
#     #subset data according to significant niches of interest
#     sig_bool = np.isin(mdata['quiche'].var[annot_key], sig_niches)
#     conn_mat = conn_mat[sig_bool, :]
#     niche_list = mdata['quiche'].var[annot_key][sig_bool]
#     nn_array = [row.nonzero()[1] for row in conn_mat] #nn around cell, where elements == 1
#     print(len(np.unique(niche_list)))
#     #create a dictionary mapping cell types to indices
#     cell_clusters_indices = {cell_type: np.where(mdata['expression'].obs[labels_key] == cell_type)[0] for cell_type in mdata['expression'].obs[labels_key].unique()}
#     #return the functional marker expression of cell types of interest within each significant niche in the original undownsampled FOV
#     func_df = []
#     for i in range(0, len(nn_array)):
#         niche = niche_list[i]
#         nn = nn_array[i]
#         for cell_type in niche.split('__'):
#             idx_cell_type = cell_clusters_indices[cell_type]
#             idx_cell_type_nn = list(set(nn).intersection(set(idx_cell_type)))
#             if len(idx_cell_type_nn) >= min_cell_count:
#                 exp = mdata['expression'][idx_cell_type_nn, :][:, markers].to_df()
#                 exp['niche'] = niche
#                 exp[labels_key] = cell_type
#                 exp['label'] = mdata['expression'].obs['combined_label']
#                 exp['niche_cell_type'] = niche + ':' + cell_type
#                 exp[fov_key] = mdata['expression'][idx_cell_type_nn].obs[fov_key].unique()[0]
#                 exp[patient_key] = mdata['expression'][idx_cell_type_nn].obs[patient_key].unique()[0]
#                 func_df.append(exp)
#     func_df = pd.concat(func_df, axis = 0)

#     # adata_func = anndata.AnnData(func_df.drop(columns = ['niche', labels_key, 'label', 'niche_cell_type', fov_key, patient_key]))
#     # adata_func.obs = func_df.loc[:, ['niche', labels_key, 'niche_cell_type', 'label', fov_key, patient_key]]
#     adata_func = anndata.AnnData(func_df.drop(columns = ['niche', labels_key, 'label', 'niche_cell_type', fov_key]))
#     adata_func.obs = func_df.loc[:, ['niche', labels_key, 'niche_cell_type', 'label', fov_key]]
#     return adata_func

def compute_functional_expression(mdata = None,
                                sig_niches = None,
                                labels_key = 'cell_cluster',
                                annot_key = 'quiche_niche',
                                fov_key = 'fov',
                                segmentation_label_key = 'label',
                                patient_key = 'Patient_ID',
                                min_cell_count = 3,
                                markers = None):
    print(len(sig_niches))
    #subset connectivities matrix according to downsampled cells
    idx = mdata['expression'].obs_names.get_indexer(mdata['spatial_nhood'].obs_names)
    conn_mat = mdata['expression'].obsp['spatial_connectivities'][idx, :] #downsampled cells x all cells
    #subset data according to significant niches of interest
    sig_bool = np.isin(mdata['quiche'].var[annot_key], sig_niches)
    conn_mat = conn_mat[sig_bool, :]
    niche_list = mdata['quiche'].var[annot_key][sig_bool]
    nn_array = [row.nonzero()[1] for row in conn_mat] #nn around cell, where elements == 1

    print(len(np.unique(niche_list)))

    assert len(sig_niches) == len(np.unique(niche_list))

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
                exp['quiche_niche'] = niche
                exp[labels_key] = cell_type
                exp[segmentation_label_key] = mdata['expression'].obs[segmentation_label_key]
                exp['quiche_niche_cell_type'] = niche + ':' + cell_type
                exp[fov_key] = mdata['expression'][idx_cell_type_nn].obs[fov_key].unique()[0]
                exp[patient_key] = mdata['expression'][idx_cell_type_nn].obs[patient_key].unique()[0]
                func_df.append(exp)
    func_df = pd.concat(func_df, axis = 0)



    adata_func = anndata.AnnData(func_df.drop(columns = ['quiche_niche', labels_key, segmentation_label_key, 'quiche_niche_cell_type', fov_key, patient_key]))
    adata_func.obs = func_df.loc[:, ['quiche_niche', labels_key, 'quiche_niche_cell_type', segmentation_label_key, fov_key, patient_key]]
    adata_func.obs = pd.merge(adata_func.obs, pd.DataFrame(mdata['quiche'].var.groupby(['quiche_niche'])['logFC'].mean()), on = ['quiche_niche'])
    return adata_func

def get_niche_expression_ind(mdata, adata, sig_niches, nn_dict, annotation_key, cluster_key, fov_key, patient_key,feature_key):
    #gives only the cell types that make up the max cell types of the niche for all of the niches
    subset_df = pd.DataFrame()
    for niche in sig_niches:
        niche_df = mdata[feature_key].obs[np.isin(mdata['milo'].var[annotation_key], niche)].copy() #subset data based on niche of interest
        fov_list = np.unique(niche_df[fov_key])
        for fov in fov_list: #for every FOV that contains the niche of interest
            cells = niche_df[niche_df[fov_key] == fov].index
            adata_subset = adata[adata.obs[fov_key] == fov].copy()
            nn_subset = nn_dict[fov][np.where(np.isin(adata_subset.obs_names, cells))[0]] #find the nn of index cells that are enriched 
            for nn in nn_subset:
                idx = [i for i in range(0, len(list(adata_subset[nn].obs[cluster_key].values))) if list(adata_subset[nn].obs[cluster_key].values)[i] in niche.split('__')] #return the cell types within each niche that make up the niche
                func_markers = adata_subset[nn][idx].to_df()
                func_markers[patient_key] = adata_subset[nn][idx].obs[patient_key].values
                func_markers[fov_key] = fov
                func_markers[annotation_key] = niche
                subset_df = pd.concat([subset_df, func_markers], axis = 0)

    other_metadata = pd.merge(subset_df.reset_index().loc[:, ['index', annotation_key]], adata.obs.reset_index(), on = 'index')
    adata_runner = anndata.AnnData(subset_df.drop(columns = [patient_key, fov_key, annotation_key]))
    adata_runner.obs = other_metadata.copy()
    return adata_runner

def binarize_functional_expression(adata_func, threshold_list):
    adata_func_binary = adata_func.to_df()
    for marker, threshold in threshold_list:
        adata_func_binary[marker] = (adata_func_binary[marker].values >= threshold).astype('int')

    adata_func_binary = anndata.AnnData(adata_func)
    adata_func_binary.obs = adata_func.obs
    return adata_func_binary

def relabel_niches_celltype(adata, annotation_key, cluster_key, sig_niches):
    adata.obs['niche_cell_type'] = pd.DataFrame(list(adata.obs[annotation_key].values)) + ' niche: ' + pd.DataFrame(list(adata.obs[cluster_key].values))
    niche_dict = dict(zip(list(adata.obs[annotation_key].values), list(adata.obs['niche_cell_type'].values)))
    groups = [i.split(':')[0] for i in pd.Series(sig_niches).map(niche_dict)]
    mapped_niches = list(dict.fromkeys([i for g in groups for i in np.unique(adata.obs['niche_cell_type']) if g in i]).keys())
    sorted_mapped_niches = [item for group in groups for item in mapped_niches if item.startswith(group)]
    return sorted_mapped_niches
