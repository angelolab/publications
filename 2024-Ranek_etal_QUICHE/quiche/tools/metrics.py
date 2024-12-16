import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import re
import itertools
import os

def plot_precision_recall(scores_df, figsize = (4,4), xlim = [-0.02, 1.02], ylim = [-0.02, 1.02], save_directory = None, filename_save = None):
    """Plots precision recall scores"""
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

def group_recall(mdata,
                 method_type = 'spatial_cluster',
                 ground_truth_niches = ['A_C_E', 'B_D'],
                 scores_df = None,
                 feature_key = 'spatial_nhood',
                 ground_key = 'DA_group',
                 labels_key = 'quiche_niche',
                 label_type = 1):
    """Computes niche recall scores

    Parameters
    mdata: mudata or anndata.AnnData
        annotated data object containing preprocessed single-cell data with spatial enrichment information
    method_type: str (default = 'other')
        string referring to the spatial enrichment method used. Either pairwise or spatial_cluster. 
    ground_truth_niches: list (default = ['A_C_E', 'B_D'])
        list containing ground truth niches
    scores_df: pd.DataFrame (default = None)
        dataframe containing niches and their significance scores
    feature_key: str (default = 'spatial_nhood)
        string referring to key in mdata object with spatial enrichment information
    labels_key: str (default = 'quiche_niche')
        string referring to the column in mdata that contains niche labels
    ground_key: str (default = 'DA_group') 
        string in mdata containing ground truth niche annotations
    label_type: int (default = 1)
        integer referring to unstructured (1) or structured (2) simulation
    ----------

    Returns
    eval_df: pd.DataFrame
        dataframe containing niche recall scores
    ----------
    """
    if method_type == 'pairwise':
        n1 = ground_truth_niches[0]
        n2 = ground_truth_niches[1]
        niche_val = 0
        pairwise_lcd1 = list(itertools.combinations(n1.split('_'), 2))
        n1 = ['__'.join(list(i)) for i in pairwise_lcd1]
        if label_type == 2:
             scores_df.index = scores_df.index.str.replace('cancer_', '', regex=True)
        if np.isin(scores_df.index, n1).sum() == len(n1):
                niche_val+=1
        pairwise_lcd2 = list(itertools.combinations(n2.split('_'), 2))
        n2 = ['__'.join(list(i)) for i in pairwise_lcd2]                
        if np.isin(scores_df.index, n2).sum() == len(n2):
            niche_val+=1
        niche_val = niche_val/2
    else:
        n1 = ground_truth_niches[0]
        n2 = ground_truth_niches[1]
        count_df = mdata[feature_key].obs.groupby([labels_key, ground_key]).size().unstack().loc[scores_df.index]
        count_df['purity_score_n1'] = count_df[n1] / count_df.sum(1)
        count_df['purity_score_n2'] = count_df[n2] / count_df.sum(1)

        niche_val = 0
        if (count_df['purity_score_n1'] > 0.5).any():
            niche_val+=1
        if (count_df['purity_score_n2'] > 0.5).any():
            niche_val+=1
        niche_val = niche_val/2

    eval_df = pd.DataFrame([niche_val], columns = ['group_recall'])
    eval_df = eval_df.melt()
    return eval_df

def compute_purity_score(adata_niche,
                        annot_key = 'kmeans_cluster_10',
                        labels_key = 'cell_cluster',
                        fov_key = 'Patient_ID',
                        condition_key = 'condition'):
    """Computes niche purity scores

    Parameters
    adata_niche: anndata.AnnData
        annotated data object containing niche-level information
    annot_key: str (default = 'kmeans_cluster_10')
        string referring to the column in adata_niche.obs with predicted niche annotations
    labels_key: str (default = 'cell_cluster')
        string referring to the column in adata_niche.obs with cell phenotype information
    fov_key: pd.DataFrame (default = 'Patient_ID')
        string referring to the column in adata_niche.obs with patient information
    condition_key: str (default = 'condition')
        string referring to key in mdata object with condition information
    ----------

    Returns
    purity_score: pd.DataFrame
        dataframe containing niche purity scores
    ----------
    """
    df = adata_niche.obs.copy()
    df['cell_count'] = 1
    count_df = df.groupby([fov_key, annot_key, labels_key]).cell_count.sum()
    total_cells = df.groupby([fov_key, annot_key]).cell_count.sum()
    proportion_df = count_df/total_cells
    count_df = count_df.reset_index()
    proportion_df = proportion_df.reset_index()
    total_counts = pd.DataFrame(df.groupby([condition_key, annot_key]).cell_count.sum().unstack().sum(0)).transpose()
    purity_score = df.groupby([condition_key, annot_key]).cell_count.sum().unstack() / total_counts.values
    return purity_score

def evaluate_purity(mdata,
                    scores_df,
                    annot_key = 'kmeans_cluster_10',
                    labels_key = 'cell_cluster',
                    fov_key = 'Patient_ID',
                    condition_key = 'DA_group',
                    feature_key = 'spatial_nhood'):      
    """Computes niche purity scores

    Parameters
    adata_niche: anndata.AnnData
        annotated data object containing niche-level information
    scores_df: pd.DataFrame
        dataframe containing niches and their significance scores
    annot_key: str (default = 'kmeans_cluster_10')
        string referring to the column in adata_niche.obs with predicted niche annotations
    labels_key: str (default = 'cell_cluster')
        string referring to the column in adata_niche.obs with cell phenotype information
    fov_key: pd.DataFrame (default = 'Patient_ID')
        string referring to the column in adata_niche.obs with patient information
    condition_key: str (default = 'condition')
        string referring to key in mdata object with condition information
    ----------

    Returns
    eval_df: pd.DataFrame
        dataframe containing niche purity scores
    ----------
    """   
    if len(scores_df.index) == 0:
        avg_purity = 0 
    else:
        _, _, purity_score = compute_purity_score(mdata[feature_key], annot_key = annot_key, labels_key = labels_key, fov_key = fov_key, condition_key = condition_key)
        purity_score = purity_score.drop('random')
        avg_purity = purity_score.loc[:, scores_df.index].max().mean()
    eval_df = pd.DataFrame([avg_purity], columns = ['avg_purity'])
    eval_df = eval_df.melt()
    return eval_df