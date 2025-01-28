
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
from sklearn.neighbors import NearestNeighbors
import time
import networkx as nx
import pandas as pd
def MED(adata,use_cls,nn,ME_var_names_np_unique,copy=False,spatial_var='spatial'):
    ## function: calculate mricroenvironment counting
#     use_cls = 'leiden'
#     nn = 8
    if copy:
        adata_use = adata.copy()
    else:
        adata_use = adata
        
    # after running，adata added an obsm called ME

    ME_var_names_np = np.arange(ME_var_names_np_unique.shape[0]).astype('str')
    # a histogram，ordered by (record frequency of every cell)
    adata_use.obsm['ME'] = np.zeros(shape=(adata_use.shape[0],ME_var_names_np.shape[0])) #frequency dim: n_cell*n_cluster



    ME_X = np.zeros(shape=(adata_use.shape[0],ME_var_names_np.shape[0]))
    time_start=time.time()
    spatial_mat = np.array(adata_use.obsm[spatial_var]).astype('int64')
    

    
#     ##############exact knn search##############
    fit = NearestNeighbors(n_neighbors=nn).fit(spatial_mat)
    m = fit.kneighbors(spatial_mat)
    m = m[0], m[1]


    #sort_neighbors
    args = m[0].argsort(axis = 1)
    add = np.arange(m[1].shape[0])*m[1].shape[1]
    I = m[1].flatten()[args+add[:,None]]
    time_end=time.time()
#     ##############knn end##############
    print('knn search time cost',time_end-time_start,'s')
    cls_array = adata_use.obs[use_cls]
    time_start=time.time()


    for i in range(I.shape[0]):

        if i%1000==0:
            time_end=time.time()
            
            print('{0} MEs,time cost {1} s, {2} MEs, {3}s left'.format(i,time_end-time_start,I.shape[0]-i,(I.shape[0]-i)*(time_end-time_start)/1000))
            time_start=time.time()
            
        cur_neighbors = I[i,:]
        
        cur_neighbors_cls = cls_array[cur_neighbors]
        cur_cls_unique,cur_cls_count = np.unique(cur_neighbors_cls,return_counts=1) #counting for each cluster
        cur_cls_idx = [np.where(ME_var_names_np_unique==c)[0][0] for c in cur_cls_unique] #c is string
        ME_X[i,cur_cls_idx] = cur_cls_count
    adata_use.obsm['ME'] = ME_X
    if copy:
        return adata_use
    else:
        return ME_X #microenvironment counting
    


from sklearn.neighbors import NearestNeighbors
import time
import networkx as nx

def MED(adata,
        use_cls,
        nn,
        ME_var_names_np_unique,
        copy=False,
        spatial_var='spatial',
        ME_key_added='ME',#obsm key of ME
        MEidx_key_added='MEidx'#uns key of ME_index
       
       ):

    

    if copy:
        adata_use = adata.copy()
    else:
        adata_use = adata
        

    ME_var_names_np = np.arange(ME_var_names_np_unique.shape[0]).astype('str')
    adata_use.obsm[ME_key_added] = np.zeros(shape=(adata_use.shape[0],ME_var_names_np.shape[0])) #frequency 
    time_start=time.time()
    spatial_mat = np.array(adata_use.obsm[spatial_var]).astype('int64')
    


    
#     ##############exact knn search##############
    fit = NearestNeighbors(n_neighbors=nn).fit(spatial_mat)
    m = fit.kneighbors(spatial_mat)
    m = m[0], m[1]


    #sort_neighbors
    args = m[0].argsort(axis = 1)
    add = np.arange(m[1].shape[0])*m[1].shape[1]
    I = m[1].flatten()[args+add[:,None]]
    
    time_end=time.time()
#     ##############knn end##############
    print('knn search time cost',time_end-time_start,'s')
    cls_array = adata_use.obs[use_cls]

    adata_use.uns[MEidx_key_added] = I
    

    ME_X = MEidx2ME(I,cls_array,ME_var_names_np_unique)
    adata_use.obsm[ME_key_added] = ME_X
    
    if copy:
        return adata_use
    else:
        return ME_X #microenvironment counting
    
    
def MEidx2ME(I,cls_array,ME_var_names_np_unique):
    ME_X = np.zeros(shape=(I.shape[0],ME_var_names_np_unique.shape[0]))
    time_start=time.time()
    for i in range(I.shape[0]):
        
        
        

        if i%1000==0:
            time_end=time.time()
            
            
            print('{0} MEs,time cost {1} s, {2} MEs, {3}s left'.format(i,time_end-time_start,I.shape[0]-i,(I.shape[0]-i)*(time_end-time_start)/1000))
            time_start=time.time()
            
        cur_neighbors = I[i,:]
        
        cur_neighbors_cls = cls_array[cur_neighbors]
        cur_cls_unique,cur_cls_count = np.unique(cur_neighbors_cls,return_counts=1) #counting for each cluster
        cur_cls_idx = [np.where(ME_var_names_np_unique==c)[0][0] for c in cur_cls_unique] #c is string
        ME_X[i,cur_cls_idx] = cur_cls_count
    return ME_X



def MED_multi(
    adata,
    use_cls,
    nn,
    copy=False,
    spatial_var='spatial',
    batch_obs='batch'
):
    unique_batch = np.unique(adata.obs[batch_obs])
    ME_var_names_np_unique = np.array(adata.obs[use_cls].cat.categories)
 
    ME_whole = np.zeros(shape=(adata.shape[0],len(ME_var_names_np_unique)))
    for batch in unique_batch:
        cur_batch_index = np.where(adata.obs[batch_obs]==batch)
        cur_adata = adata[cur_batch_index]
        cur_ME = MED(cur_adata,use_cls,nn,ME_var_names_np_unique,False,spatial_var)
        ME_whole[cur_batch_index] = cur_ME
    adata.obsm['ME'] = ME_whole
    if copy:
        return adata
    else:
        return ME_whole
    


import pyemd
from scipy.spatial.distance import *

def get_cls_center(adata,cls_key,embed_key):


    cat = np.array(adata.obs[cls_key].cat.categories)
    embed_mat = np.array(adata.obsm[embed_key])
    cls_array = np.array(adata.obs[cls_key])

    mean_array = np.zeros(shape=(len(cat),embed_mat.shape[1]))
    for i in range(len(cat)):
        c = cat[i]
        mean_array[i] = np.mean(embed_mat[np.where(cls_array==c)])
    dist_mat = squareform(pdist(mean_array))
    return dist_mat
    
    
def get_ground_distance(adata,method='paga_guided_umap',cls_key=None,embed_key=None,connect_threshold=0.5): 
    
    n_cls = np.array(adata.obs[cls_key].cat.categories).shape[0]
    if method == 'paga_graph':
        adj_mat_use = adata.uns['paga']['connectivities'].toarray()

        adj_mat_use = 1/adj_mat_use

        np.fill_diagonal(adj_mat_use,0)
        G = nx.from_numpy_matrix(adj_mat_use)

        len_path = dict(nx.all_pairs_dijkstra_path_length(G))
        ground_distance_mat = np.inf*np.ones(shape=(len(len_path),len(len_path)))

        for i in range(ground_distance_mat.shape[0]):
            for j in range(ground_distance_mat.shape[1]):
                ground_distance_mat[i,j] = len_path[i][j]
        np.fill_diagonal(ground_distance_mat,0)
    elif method =='paga':
        adj_mat_use = adata.uns['paga']['connectivities'].toarray()
#         adj_mat_use[adj_mat_use==0] = -np.inf
        adj_mat_use = 1-adj_mat_use
        np.fill_diagonal(adj_mat_use,0)
        ground_distance_mat = adj_mat_use
    elif method =='euclidean':   #lower bound of MED(approximation)
        ground_distance_mat = squareform(pdist(adata.uns['paga']['pos']))
    elif method == 'uniform':
#         sz = adata.uns['paga']['pos'].shape[0]
        ground_distance_mat = np.ones(shape=(n_cls,n_cls))
        np.fill_diagonal(ground_distance_mat,0)
    elif method == 'embed':
#         cls_key = ''
        if cls_key is None or embed_key is None:
            print('need cls_ley or embed_key')
            return
        ground_distance_mat = get_cls_center(adata,cls_key,embed_key)
    elif method=='paga_guided_umap':
#         connect_threshold = 0.5
        
        embed_center = adata.uns['paga']['pos']
        embed_center_distmat = squareform(pdist(embed_center))

        adjacent_paga = adata.uns['paga']['connectivities'].toarray()
        adjacent_paga_bool = (adjacent_paga>=connect_threshold).astype('int')



        disconnected_cluster_idxs = np.where(np.sum(adjacent_paga_bool,axis=0)==0)[0]
        for dis_idx in disconnected_cluster_idxs:
            cur_argmax = np.argmax(adjacent_paga[:,dis_idx])
            adjacent_paga_bool[dis_idx,cur_argmax] = 1
            adjacent_paga_bool[cur_argmax,dis_idx] = 1

        paga_guided_dist = np.multiply(embed_center_distmat,adjacent_paga_bool)
        G = nx.from_numpy_matrix(paga_guided_dist)
        ground_distance_mat = -np.inf*np.ones_like(paga_guided_dist)

        
        len_path = dict(nx.all_pairs_dijkstra_path_length(G))
        for i in len_path.keys():
            for j in len_path[i].keys():
                ground_distance_mat[i,j] = len_path[i][j]
        np.fill_diagonal(ground_distance_mat,0)
        max_val = np.max(ground_distance_mat)
        np.nan_to_num(ground_distance_mat,copy=False,neginf=max_val*2)

    adata.uns['GD'] = ground_distance_mat
    return ground_distance_mat
    
    
def MED_phEMD(adata,  
              GD_method='euclidean',
              MED_knn=30,
              CT_obs='clusters',
              ifspatialplot=True,
              OT_method='pyemd',
              ME_precompyted=False,
              GD_precomputed=False
             ):
    
    
    # calculate ground distance, paga requried
    if not GD_precomputed:
        ground_distance_mat = get_ground_distance(adata,method=GD_method)
    else:
        ground_distance_mat = adata.uns['GD']
    # get ME histogram
    if not ME_precompyted:
        MED(adata,CT_obs,MED_knn)
    ME_mat = adata.obsm['ME']


    # calculate EMD matrix
    import time
    ME_mat = adata.obsm['ME']
    ME_dist_EMD = np.zeros(shape=(ME_mat.shape[0],ME_mat.shape[0]))
    time_start_all = time.time()
    time_cost = 0
    time_start = 0
    time_end = 0
    for i in range(ME_mat.shape[0]):
        if i%100==0:
            time_end = time.time()
            time_cost = time_end-time_start
            print(i,time_cost)
            time_start = time.time()

        for j in range(i,ME_mat.shape[0]):
            first_histogram = ME_mat[i,:]
            second_histogram = ME_mat[j,:]
            first_histogram = first_histogram/np.sum(first_histogram)
            second_histogram = second_histogram/np.sum(second_histogram)
            if OT_method=='pyemd':
                cur_dist = pyemd.emd(first_histogram,second_histogram,ground_distance_mat)
    #         print(i,j,cur_dist)
            elif OT_method=='POT':
                first_histogram = first_histogram/np.sum(first_histogram)
                second_histogram = second_histogram/np.sum(second_histogram)
                cur_dist = ot.emd2(first_histogram,second_histogram,ground_distance_mat)
            ME_dist_EMD[i,j] = cur_dist
            ME_dist_EMD[j,i] = cur_dist
    adata.obsm['X_ME_EMD_mat'] = ME_dist_EMD.copy()
    time_end_all = time.time()
    time_cost_all = time_end_all-time_start_all
    print(f'EMD distance matrix cost {time_cost_all}s')
    return adata  
from multiprocessing import Process, Pool
def local_func(para):
    i,j,ME_mat,ground_distance_mat = para
    first_histogram = ME_mat[i,:]
    second_histogram = ME_mat[j,:]
    first_histogram = first_histogram/np.sum(first_histogram)
    second_histogram = second_histogram/np.sum(second_histogram)
    cur_dist = pyemd.emd(first_histogram,second_histogram,ground_distance_mat)
    return cur_dist
def MED_phEMD_mp(adata,
              GD_method='euclidean',
              MED_knn=30,
              CT_obs='clusters',
              ifspatialplot=True,
              OT_method='pyemd',
              ME_precompyted=False,
              GD_precomputed=False,
              mp=100
             ):


    # calculate ground distance, paga requried
    if not GD_precomputed:
        ground_distance_mat = get_ground_distance(adata,method=GD_method)
    else:
        ground_distance_mat = adata.uns['GD']
    # get ME histogram
    if not ME_precompyted:
        MED(adata,CT_obs,MED_knn)
    ME_mat = adata.obsm['ME']
    


    # calculate EMD matrix
    import time
    time_cost = 0
    time_start = time.time()
    ME_mat = adata.obsm['ME']
    ME_dist_EMD = np.zeros(shape=(ME_mat.shape[0],ME_mat.shape[0]))

    pool = Pool(mp)
    
    # prepare parameter list for multiprocessing
    para_list = []
    for i in range(ME_mat.shape[0]):
        for j in range(i,ME_mat.shape[0]):
            para_list.append((i,j,ME_mat,ground_distance_mat))
    # run mp
    rst_list = pool.map(local_func,para_list)
    
    # fill matrix
    k=0
    for i in range(ME_mat.shape[0]):
        for j in range(i,ME_mat.shape[0]):
            ME_dist_EMD[i,j] = rst_list[k]
            ME_dist_EMD[j,i] = rst_list[k]
            k+=1
    

    adata.obsm['X_ME_EMD_mat'] = ME_dist_EMD.copy()
    time_end = time.time()
    time_cost = time_end-time_start
    print(f'EMD distance matrix cost {time_cost}s')

    return adata

from multiprocessing import Process, Pool

def MED_phEMD_mp_v2(adata,
              GD_method='euclidean',
              MED_knn=30,
              CT_obs='clusters',
              ifspatialplot=True,
              OT_method='pyemd',
              ME_precompyted=False,
              GD_precomputed=False,
              mp=100
             ):


    # calculate ground distance, paga requried
    if not GD_precomputed:
        ground_distance_mat = get_ground_distance(adata,method=GD_method)
    else:
        ground_distance_mat = adata.uns['GD']
    # get ME histogram
    if not ME_precompyted:
        MED(adata,CT_obs,MED_knn)
    ME_mat = adata.obsm['ME']



    # calculate EMD matrix
    import time
    time_cost = 0
    time_start = time.time()
    ME_mat = adata.obsm['ME']
    ME_dist_EMD = np.zeros(shape=(ME_mat.shape[0],ME_mat.shape[0]))

    pool = Pool(mp)

    def local_func(para):
        # i,j,ME_mat,ground_distance_mat = para
        i,j = para
        
        first_histogram = ME_mat[i,:]
        second_histogram = ME_mat[j,:]
        first_histogram = first_histogram/np.sum(first_histogram)
        second_histogram = second_histogram/np.sum(second_histogram)
        cur_dist = pyemd.emd(first_histogram,second_histogram,ground_distance_mat)
        return cur_dist
    
    para_list = []
    for i in range(ME_mat.shape[0]):
        for j in range(i,ME_mat.shape[0]):
            para_list.append((i,j))
            
    rst_list = pool.map(local_func,para_list)

    k=0
    for i in range(ME_mat.shape[0]):
        for j in range(i,ME_mat.shape[0]):
            ME_dist_EMD[i,j] = rst_list[k]
            ME_dist_EMD[j,i] = rst_list[k]
            k+=1


    adata.obsm['X_ME_EMD_mat'] = ME_dist_EMD.copy()
    time_end = time.time()
    time_cost = time_end-time_start
    print(f'EMD distance matrix cost {time_cost}s')

    return adata



def cal_ME_heterogeneity(
    adata,
    copy=False,
    key_added='ME_heter',
    ME_key='ME' 
):
    if copy:
        adata_use = adata.copy()
    else:
        adata_use = adata
    ME_mat = adata_use.obsm[ME_key]
    GD_mat = adata_use.uns['GD']
    ME_heter_list = []
    for i in range(ME_mat.shape[0]):  #for each cell
        cur_ME = ME_mat[i,:]
        cur_heter = 0
        for j in range(ME_mat.shape[1]): #for each cluster
            cur_dist_add = np.sum(GD_mat[j]*cur_ME*cur_ME[j]) #edge count*edge value(ground distance)；cell frequency information is included in edge count
            cur_heter+=cur_dist_add
        ME_heter_list.append(cur_heter)
    if copy:
        adata_use.obs[key_added] = ME_heter_list
        return adata_use
    else:
        adata_use.obs[key_added] = ME_heter_list
        return 
    
    
    
def cal_permuted_ME_heterogeneity(
    adata,
    use_cls,
    ME_var_names_np_unique,
    copy=False,
    key_added='ME_heter', 
    MEidx_key='MEidx_2x',
    permute_num=7, 
    n_perms=100
):
    
    if copy:
        adata_use = adata.copy()
    else:
        adata_use = adata
    
    cls_array = adata_use.obs[use_cls]
    
    print(adata_use,MEidx_key)
    MEidx_mat = adata_use.uns[MEidx_key]
    
    heter_perm_mat = np.zeros(shape=(adata_use.shape[0],n_perms))
    for i_perm in range(n_perms):
        chose_MEidx_mat = np.zeros(shape=(MEidx_mat.shape[0],permute_num))
        
        for i in range(MEidx_mat.shape[0]):
            chose_MEidx_mat[i] = np.random.choice(MEidx_mat[i],permute_num,replace=False)

        chose_MEidx_mat = chose_MEidx_mat.astype('int')
        ME_mat = MEidx2ME(chose_MEidx_mat,cls_array,ME_var_names_np_unique)
        GD_mat = adata_use.uns['GD']

        adata_use.obsm['ME_perm'] = ME_mat
        cal_ME_heterogeneity(adata_use,key_added='heter_perm',ME_key='ME_perm')
        heter_perm_mat[:,i_perm] = adata_use.obs['heter_perm'].copy()
    adata_use.obsm[key_added] = heter_perm_mat

    if copy:
        return adata_use
    else:
        return 

    
def prepare_boundary_id(adata):
    exp_1_list = np.array(adata.obs['ME_heter'])
    spatial_mat = adata.obsm['spatial']


    exp_1_list =(exp_1_list-exp_1_list.min())/(exp_1_list.max()-exp_1_list.min())
    spatial_dist = 1/squareform(pdist(spatial_mat))


    polarity_1_max_list = []

    for i in range(len(exp_1_list)):
        dist_array = spatial_dist[i,:]
        polarity_1_array = dist_array*exp_1_list

        polarity_1_array = polarity_1_array[np.logical_not(np.isposinf(polarity_1_array))]

        polarity_1_array = polarity_1_array[np.logical_not(np.isnan(polarity_1_array))]



        polarity_1_max = np.max(polarity_1_array)



        polarity_1_max_list.append(polarity_1_max)

    adata.obs['polarity_1_max_list']=polarity_1_max_list


def test_polar(adata,gene):
    gene_profile = adata[:,test_gene].X.toarray()[:,0]
    soft_dist_profile = adata.obs['polarity_1_max_list']
    pos_idx = (gene_profile>np.mean(gene_profile))
    pos_dist_profile = soft_dist_profile[pos_idx]
    logfc = np.log10(np.mean(pos_dist_profile)/np.mean(soft_dist_profile))
    pvalue = ttest_ind(soft_dist_profile,pos_dist_profile)[1]
    pvalue = -np.log10(pvalue)
    return pvalue,logfc


def plot_ME_cluster(adata,cls_idx,ME_cls_obs,cls_obs):

    pd_dict = {
        'type':[],
        'count':[],
    }
    ME_cluster = np.array(adata[adata.obs[ME_cls_obs]==cls_idx].obsm['ME'])
    #ME_var_names_np = np.unique(adata.obs[cls_obs])
    ME_var_names_np = np.array(adata.obs[cls_obs].cat.categories)
    for i in range(ME_cluster.shape[0]):
        for j in range(ME_cluster.shape[1]):
            cur_val = ME_cluster[i,j]
            cur_type = ME_var_names_np[j]
            pd_dict['type'].append(cur_type)
            pd_dict['count'].append(cur_val)
    pd_df = pd.DataFrame(pd_dict)



    sns.set_style('white')
    fig,ax = plt.subplots(1,1,figsize=(10,3))

    sns.barplot(x='type', y='count', data=pd_df,palette=adata.uns['{0}_colors'.format(cls_obs)])
    for item in ax.get_xticklabels():
        item.set_rotation(45)
    plt.title('ME cluster {0}'.format(cls_idx))


def plot_ME_cluster_stacked(adata,intest_ME,ct_cls_obs,ME_cls_obs,ct_color,width=0.5):
#example: plot_ME_cluster_stacked(adata,intest_ME,'ct_name','leiden_EMD',adata.uns['ct_name_colors'],0.5)
    pd_dict = {
        'ME':[],
        'ct':[],
        'count':[]
    }


    ct_order = np.array(adata.obs[ct_cls_obs].cat.categories)
    for i in intest_ME:
        ME_cluster = np.array(adata[adata.obs[ME_cls_obs]==i].obsm['ME'])
        for j in range(ME_cluster.shape[0]):
            for k in range(ME_cluster.shape[1]):
                cur_val = ME_cluster[j,k]
                cur_ct = ct_order[k]
                pd_dict['ME'].append(i)
                pd_dict['ct'].append(cur_ct)
                pd_dict['count'].append(cur_val)
    pd_df = pd.DataFrame(pd_dict)

    fig, ax = plt.subplots()
    # width = 0.5

    prev_me_mean_list = [0]*len(intest_ME)
    for i in range(len(ct_order)):
        ct = ct_order[i]
        cur_color = ct_color[i]
        cur_pd_ct = pd_df[pd_df['ct']==ct]
        cur_me_mean_list = []
        cur_me_std_list = []
        for ME in intest_ME:
            cur_pd_me = cur_pd_ct[cur_pd_ct['ME']==ME]
            cur_mean = np.mean(cur_pd_me['count'])
            cur_std = np.std(cur_pd_me['count'])
            cur_me_mean_list.append(cur_mean)
            cur_me_std_list.append(cur_std)

 
        ax.bar(intest_ME, cur_me_mean_list, width, label=ct,bottom=prev_me_mean_list,color=cur_color)

        prev_me_mean_list = (np.array(prev_me_mean_list)+np.array(cur_me_mean_list)).tolist()
    ax.set_ylabel('count')
    ax.set_title('ME stacked bar plot')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    plt.show()

