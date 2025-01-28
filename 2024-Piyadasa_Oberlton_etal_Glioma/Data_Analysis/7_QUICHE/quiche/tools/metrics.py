import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
import matplotlib.pyplot as plt
from sklearn.metrics import accuracy_score, roc_curve, auc, f1_score, recall_score, precision_score
import quiche as qu
import re
import itertools
import os
import networkx as nx
from itertools import combinations
#evaluate precision, recall, AUPRC
def evaluate_precision_recall(mdata = None,
                                scores_df = None,
                                thresholds = np.array([1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]),
                                labels_key = 'mask_name',
                                ground_truth_key = 'ground_truth_DA',
                                group = 'immune1',
                                feature_key = 'quiche',
                                **args):
    y_true = _access_ytrue(mdata = mdata, ground_truth_key = ground_truth_key, labels_key = labels_key, group = group, feature_key = feature_key)
    precision = []
    recall = []
    for threshold in thresholds:
        y_pred = np.select([((mdata[feature_key].var['SpatialFDR'] <= threshold) & (mdata[feature_key].var[labels_key] == group))],
                            [1], default=0)

        tp = np.sum((y_pred == 1) & (y_true == 1))
        fp = np.sum((y_pred == 1) & (y_true == 0))
        fn = np.sum((y_pred == 0) & (y_true == 1))

        # Calculate precision and recall
        precision.append(tp / (tp + fp) if (tp + fp) != 0 else 0)
        recall.append(tp / (tp + fn) if (tp + fn) != 0 else 0)
    
    eval_df =  pd.DataFrame(thresholds, columns = ['sig_cutoff'])
    eval_df['precision'] = precision
    eval_df['recall'] = recall
    eval_df = eval_df.melt(id_vars = 'sig_cutoff')
    prc_auc = auc(recall, precision)
    prc_auc = pd.DataFrame([np.nan, 'AUPRC', prc_auc], index = ['sig_cutoff', 'variable', 'value']).transpose()
    eval_df = pd.concat([eval_df, prc_auc], axis = 0)

    return eval_df

def _access_ytrue(mdata = None, ground_truth_key = 'ground_truth_DA', labels_key = 'mask_name', group = 'immune1', feature_key = 'milo'):
    try:
        mdata[feature_key].var['SpatialFDR'][mdata[feature_key].var['SpatialFDR'].isna()] = 1.0
    except:
        pass
    mdata[feature_key].var[ground_truth_key] = mdata['expression'].obs[ground_truth_key]
    mdata[feature_key].var[labels_key] = mdata['expression'].obs[labels_key].values

    idx_true = np.where((mdata['expression'].obs[ground_truth_key] == 1) &(mdata['expression'].obs[labels_key] == group))[0]
    y_true = np.zeros(len(mdata['expression'].obs.index))
    y_true[idx_true] = 1
    return y_true

def plot_precision_recall(scores_df, figsize = (4,4), xlim = [-0.02, 1.02], ylim = [-0.02, 1.02], save_directory = None, filename_save = None):
    fig, ax = plt.subplots(figsize=figsize)
    recall = scores_df[scores_df['variable'] == 'recall']['value'].values
    precision = scores_df[scores_df['variable'] == 'precision']['value'].values
    prc_auc = scores_df[scores_df['variable'] == 'AUPRC']['value'][0]
    ax.plot(recall, precision, lw=2, label=f'AUPRC {prc_auc:0.2f})', color='navy', marker = '.')
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.grid(True)
    ax.set_xlabel('Recall')
    ax.set_ylabel('Precision')
    ax.legend(loc="upper right", bbox_to_anchor=(1.0,0.95))
    if filename_save is not None:
        plt.savefig(os.path.join(save_directory, filename_save + '.pdf'), bbox_inches = 'tight')

def _normalize_niche(niche):
    parts = niche.split("__")
    parts.sort()
    if len(parts) == 1:
        parts = re.sub(r'[\d\.]*_', '', parts[0])
        parts = [parts, parts]
    return "__".join(parts)

def evaluate_jaccard(scores_df = None,
                     y_true = None,
                    cell_types = None,
                    **args):
    pairwise_lcd = list(itertools.combinations_with_replacement(cell_types, 2))
    pairwise_lcd = ["__".join(combination) for combination in pairwise_lcd]
    pairwise_lcd = [_normalize_niche(niche) for niche in pairwise_lcd]
    if len(scores_df.index) == 0:
        jaccard_index = 0
    else:
        y_pred = [_normalize_niche(niche) for niche in list(scores_df.index)]
        y_true = [_normalize_niche(niche) for niche in y_true]

        y_pred_simplified = []
        for y_pred_i in y_pred:
            for lcd_i in pairwise_lcd:
                if lcd_i in y_pred_i:
                    y_pred_simplified.append(lcd_i)

        inter = len(set(y_pred_simplified).intersection(set(y_true)))
        union = len(set(y_pred_simplified).union(set(y_true)))
        jaccard_index = inter / union
    eval_df = pd.DataFrame([jaccard_index], columns = ['jaccard_index'])
    eval_df = eval_df.melt()
    return eval_df

def compute_purity_score(adata_niche,
                        annot_key = 'kmeans_cluster_10',
                        labels_key = 'mask_name',
                        fov_key = 'Patient_ID',
                        condition_key = 'condition'):

    df = adata_niche.obs.copy()
    #index for counting cells
    df['cell_count'] = 1

    count_df = df.groupby([fov_key, annot_key, labels_key]).cell_count.sum()

    total_cells = df.groupby([fov_key, annot_key]).cell_count.sum()

    #frequency of cell types per cluster per patient
    proportion_df = count_df/total_cells

    count_df = count_df.reset_index()
    proportion_df = proportion_df.reset_index()

    #frequency of cells in each condition per cluster
    total_counts = pd.DataFrame(df.groupby([condition_key, annot_key]).cell_count.sum().unstack().sum(0)).transpose()
    purity_score = df.groupby([condition_key, annot_key]).cell_count.sum().unstack() / total_counts.values

    return count_df, proportion_df, purity_score

def evaluate_purity(mdata,
                    scores_df,
                    annot_key = 'kmeans_cluster_10',
                    labels_key = 'mask_name',
                    fov_key = 'Patient_ID',
                    condition_key = 'condition',
                    feature_key = 'spatial_nhood'):            
    if len(scores_df.index) == 0:
        avg_purity = 0 
    else:
        _, _, purity_score = compute_purity_score(mdata[feature_key], annot_key = annot_key, labels_key = labels_key, fov_key = fov_key, condition_key = condition_key)
        avg_purity = purity_score.loc[:, scores_df.index].max().mean()
    eval_df = pd.DataFrame([avg_purity], columns = ['avg_purity'])
    eval_df = eval_df.melt()
    return eval_df

def compute_patient_proportion(mdata,
                                niches = None,
                                feature_key = 'quiche',
                                annot_key = 'quiche_niche',
                                patient_key = 'Patient_ID',
                                design_key = 'Relapse',
                                patient_niche_threshold = 3):
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
    cov_count_df['med_logFC'] = mdata[feature_key].var.groupby(annot_key)['logFC'].median().loc[cov_count_df[annot_key]].values
    cov_count_df['mean_logFC'] = mdata[feature_key].var.groupby(annot_key)['logFC'].mean().loc[cov_count_df[annot_key]].values
    cov_count_df['pval'] = mdata[feature_key].var.groupby(annot_key)['SpatialFDR'].median().loc[cov_count_df[annot_key]].values
    cov_count_df = pd.merge(cov_count_df, avg_counts, on = [annot_key, design_key])

    cov_total = mdata[feature_key].var.groupby([design_key])[patient_key].nunique()
    cov_total = pd.DataFrame(cov_total)
    cov_total.columns = ['patient_cov']
    cov_count_df = pd.merge(cov_count_df, cov_total, on = [design_key])
    cov_count_df['prop_cov'] = cov_count_df['patient_count'].div(cov_count_df['patient_cov'])
    return cov_count_df

def compute_niche_network_centrality(G):
    eigen_cent = nx.eigenvector_centrality(G, weight = 'weight')
    eigen_cent = pd.DataFrame(eigen_cent, index = ['value']).transpose()
    eigen_cent['centrality'] = 'eigenvector'
    between_cent = nx.betweenness_centrality(G,weight='inv_weight')
    between_cent = pd.DataFrame(between_cent, index = ['value']).transpose()
    between_cent['centrality'] = 'betweenness'
    # Sort nodes by betweenness centrality
    centrality_df = pd.concat([eigen_cent, between_cent], axis = 0)
    return centrality_df

def compute_niche_network(cov_count_df = None,
                        colors_dict = None,
                        lineage_dict = None,
                        annot_key = 'quiche_niche'):

    node_dict = {}
    for index, row in cov_count_df.iterrows():
        nodes = row[annot_key].split('__')
        for node in nodes:
            if node in node_dict:
                # increase patient count and accumulate mean_logFC
                node_dict[node]['mean_logFC'].append(row['mean_logFC'])
            else:
                node_dict[node] = {
                    'mean_logFC': [row['mean_logFC']]
                }

    # Calculating average mean_logFC for each node
    for node in node_dict:
        node_dict[node]['mean_logFC'] = sum(node_dict[node]['mean_logFC']) / len(node_dict[node]['mean_logFC'])

    # Now calculate edge weights
    # Use a dictionary to keep track of the edges and their weights
    edge_dict = {}
    for index, row in cov_count_df.iterrows():
        nodes = row[annot_key].split('__')
        edges = combinations(nodes, 2)
        for edge in edges:
            edge = tuple(sorted(edge))  # ensure the edge (node1, node2) is the same as (node2, node1)
            if edge in edge_dict:
                edge_dict[edge] += 1  # increase count
            else:
                edge_dict[edge] = 1  # initialize

    # Now, let's create a graph
    G = nx.Graph()

    # Add nodes with size and color
    for node in node_dict:
        color = colors_dict[lineage_dict[node]]
        G.add_node(node, color=color)

    # Add edges with weight
    for edge in edge_dict:
        G.add_edge(edge[0], edge[1], weight=edge_dict[edge])


    # Add a new edge attribute that is the inverse of the weight
    for u, v, data in G.edges(data=True):
        data['inv_weight'] = 1 / data['weight']

    return G
