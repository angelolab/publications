import pandas as pd
import numpy as np
from collections import defaultdict
def get_unique_dict(array):
    unique_rst = np.unique(array,return_counts=1)
    
 
    return defaultdict(int,zip(unique_rst[0].tolist(),unique_rst[1].tolist()))

def get_bernoulli_entropy(p):
    from scipy.stats import entropy
    return entropy([p,1-p],base=2)


# adata
# heter_obs = heter_key
# ME_cls_obs = 'leiden'
# batch_obs = 'batch'
def get_EH_pd(adata,heter_obs,ME_cls_obs,batch_obs):
    heter_array = np.array(adata.obs[heter_obs])
    ME_cls_array = np.array(adata.obs[ME_cls_obs])
    batch_array = np.array(adata.obs[batch_obs])
    batch_keys = adata.obs[batch_obs].cat.categories
    pd_dict = {
        'ME_cls_no':[],
        'heterogeneity':[],
        'entropy_norm':[],
        'entropy':[],
        'size':[],
        'batch_0_p':[],
        'batch_1_p':[],
        'batch_0_p_norm':[],
        'batch_1_p_norm':[],
        
        
        # 'cls_count':[]

    }
    batch_unique_dict = get_unique_dict(batch_array)
    batch_unique_array = np.array([batch_unique_dict[batch_keys[0]],batch_unique_dict[batch_keys[1]]])
    for i in adata.obs[ME_cls_obs].cat.categories:
        ME_cls_idx = np.where(ME_cls_array==i)[0]
        cur_heter_array = heter_array[ME_cls_idx]
        cur_batch_array = batch_array[ME_cls_idx]
        cur_batch_unique_dict = get_unique_dict(cur_batch_array)
        cur_batch_unique_array = np.array([cur_batch_unique_dict[batch_keys[0]],cur_batch_unique_dict[batch_keys[1]]])

        cur_batch_array_norm = cur_batch_unique_array/batch_unique_array
        cur_batch_p = cur_batch_unique_array/cur_batch_unique_array.sum()
        cur_batch_p_norm = cur_batch_array_norm/cur_batch_array_norm.sum()

        cur_entropy = get_bernoulli_entropy(cur_batch_p[0])
        cur_entropy_norm = get_bernoulli_entropy(cur_batch_p_norm[0])
        cur_heter_mean = np.mean(cur_heter_array)

        pd_dict['ME_cls_no'].append(i)
        pd_dict['heterogeneity'].append(cur_heter_mean)
        pd_dict['entropy_norm'].append(cur_entropy_norm)
        pd_dict['entropy'].append(cur_entropy)
        pd_dict['size'].append(ME_cls_idx.shape[0])
        
        pd_dict['batch_0_p'].append(cur_batch_p[0])
        pd_dict['batch_1_p'].append(cur_batch_p[1])
        pd_dict['batch_0_p_norm'].append(cur_batch_p_norm[0])
        pd_dict['batch_1_p_norm'].append(cur_batch_p_norm[1])
        
        
    return pd.DataFrame(pd_dict)

