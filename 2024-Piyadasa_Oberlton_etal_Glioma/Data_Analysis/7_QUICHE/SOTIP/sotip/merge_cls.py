
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import scanpy as sc
import palettable
import anndata as ad
import scanpy.external as sce
import sys
from scipy.stats import *
from sklearn.metrics import *

# cls_key = 'leiden_EMD'
# neighbor_key = 'neighbors_EMD'
# thresh = 0.5
def merge_cls_paga(adata,cls_key = 'leiden_EMD',neighbor_key = 'neighbors_EMD',thresh = 0.5,min_cls=2,paga_plot=True):
    if f'{cls_key}_merge_colors' in adata.uns:
        del adata.uns[f'{cls_key}_merge_colors']
    adata.obs[f'{cls_key}_merge'] = adata.obs[cls_key].copy()
    merge_cls_list = []
    while True:
        
        sc.tl.paga(adata,groups=f'{cls_key}_merge',neighbors_key=neighbor_key)
        if paga_plot:
            sc.pl.paga(adata,color=[f'{cls_key}_merge'],threshold=0)
        
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
        if paga_plot==True:
            del adata.uns[f'{cls_key}_merge_colors'][merge_cls_idx[0]]

        new_ME_cls = new_ME_cls.cat.remove_unused_categories()
        merge_cls_list.append(new_ME_cls.copy())
        adata.obs[f'{cls_key}_merge'] = new_ME_cls
        
    adata.uns['merge_cls_list'] = merge_cls_list



# merge_cls: merge cls_i to cls_j
# input: old_cls(categories), cls_i,cls_j.
def merge_cls(old_cls,cls_i,cls_j):
    # cls_i = '3'
    # cls_j = '8'
    
    # old_cls = cur_ME_cls
    old_cls[old_cls==cls_i] = cls_j
    print(f'merged {cls_i} to {cls_j}')
    return old_cls


