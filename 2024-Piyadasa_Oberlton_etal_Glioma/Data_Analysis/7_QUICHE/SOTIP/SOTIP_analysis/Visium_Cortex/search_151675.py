data_name = '151675'

from sotip import *
import numpy as np
import scanpy as sc
import gc
import os
os.makedirs(f'search_{data_name}', exist_ok=True)

def select_preprocess(adata,HVG,scale,PCA,knn,leiden_res):
    if HVG==False:
        pass
    else:
        sc.pp.highly_variable_genes(adata, flavor=HVG, n_top_genes=2000)
    if scale==False:
        pass
    else:
        sc.pp.scale(adata)

    sc.pp.pca(adata,n_comps=PCA)
    sc.pp.neighbors(adata,n_pcs=PCA,n_neighbors=knn)
    sc.tl.umap(adata)
    sc.tl.leiden(adata,resolution=leiden_res)
    adata_valid = adata[np.logical_not(adata.obs['Region'].isna())]
    ari = adjusted_rand_score(adata_valid.obs['Region'],adata_valid.obs['leiden'])
    return ari

def get_emd_distmat(adata,ME_knn):
    knn = ME_knn
    spatial_var='spatial'
    cls_key='leiden'
    ME_var_names_np_unique = np.array(adata.obs[cls_key].cat.categories) 

    MED(adata,use_cls=cls_key,nn=knn,copy=False,ME_var_names_np_unique=ME_var_names_np_unique,spatial_var=spatial_var) 
    sc.tl.paga(adata,groups=cls_key)
    sc.pl.paga_compare(adata,basis='X_umap',show=True)    
    gd_method = 'paga_guided_umap'
    gd = get_ground_distance(adata,method=gd_method,cls_key=cls_key,embed_key=None,connect_threshold=0.5)  
    heter_key = 'ME_heter_{0}_{1}'.format(cls_key,gd_method)
    cal_ME_heterogeneity(adata,copy=False,key_added=heter_key) 
    adata_phEMD = MED_phEMD_mp(
        adata.copy(),
        GD_method=gd_method,
        MED_knn=knn,
        CT_obs=cls_key,
        ifspatialplot=False,
        OT_method='pyemd',
        ME_precompyted=True,
        GD_precomputed=True,
        mp=200
    )
    return adata_phEMD.obsm['X_ME_EMD_mat']

def get_sotip_ari(adata,n_neighbors,LIBD_cls_num=7):
    knn_indices, knn_dists, forest = sc.neighbors.compute_neighbors_umap( adata.obsp['ME_EMD_mat'], n_neighbors=n_neighbors, metric='precomputed' )
    adata.obsp['distances'], adata.obsp['connectivities'] = sc.neighbors._compute_connectivities_umap(
        knn_indices,
        knn_dists,
        adata.shape[0],
        n_neighbors, # change to neighbors you plan to use
    )
    adata.uns['neighbors_EMD'] = adata.uns['neighbors'].copy()

    sc.tl.leiden(adata,neighbors_key='neighbors_EMD',key_added='leiden_EMD',resolution=1)
    sc.tl.paga(adata,groups='leiden_EMD',neighbors_key='neighbors_EMD')
    merge_cls_paga(adata,thresh=0,min_cls=LIBD_cls_num,paga_plot=False)

    adata_valid = adata[np.logical_not(adata.obs['Region'].isna())]
    cur_ari = adjusted_rand_score(adata_valid.obs['Region'],adata_valid.obs['leiden_EMD_merge'])
    return cur_ari


h5ad_path = 'path/to/SpatialLIBD/h5ad_preprocess'


h5ad_file = f'{h5ad_path}/{data_name}.h5ad'

_hvg = False
_scale = False
# PCA_list = [15,30,50,100,200]
# PCA_list = [15,30,50,100,200]

# PCA_list = [30,50,100,200]

PCA_list = [100]
# cell_cls_knn_list = [100]
# PCA_list = [200]
cell_cls_knn_list = [15,30,50,100]
cell_cls_knn_list = [15]
cell_cls_res_list = [1,2]
ME_knn_list = [10,30,50]
# EMD_neighbors_list = np.arange(400,1100,100)
EMD_neighbors_list = np.arange(100,1100,100)
# EMD_neighbors_list = [500]


# PCA_list = [15,200]
# cell_cls_knn_list = [15,100]
# ME_knn_list = [7,100]
# EMD_neighbors_list = [100,1000]
pd_dict = {
    'PCA':[],
    'Cell_cls_knn':[],
    'Cell_cls_res':[],
    'ME_knn':[],
    'EMD_neighbors':[],
    'leiden_ari':[],
    'sotip_ari':[]
}

for _PCA in PCA_list:
    for _knn in cell_cls_knn_list:
        for leiden_res in cell_cls_res_list:
            for ME_knn in ME_knn_list:
                # adata = sc.read_h5ad('../spagcn_processed_wiHE.h5ad')
                adata = sc.read_h5ad(h5ad_file)
                leiden_ari = select_preprocess(adata,_hvg,_scale,_PCA,_knn,leiden_res)
                emd_distmat = get_emd_distmat(adata,ME_knn)
                adata.obsp['ME_EMD_mat'] = emd_distmat
                for n_neighbors in EMD_neighbors_list:
                    sotip_ari = get_sotip_ari(adata,n_neighbors,LIBD_cls_num=7)
                    pd_dict['PCA'].append(_PCA)
                    pd_dict['Cell_cls_knn'].append(_knn)
                    pd_dict['Cell_cls_res'].append(leiden_res)
                    pd_dict['ME_knn'].append(ME_knn)
                    pd_dict['EMD_neighbors'].append(n_neighbors)
                    pd_dict['leiden_ari'].append(leiden_ari)
                    pd_dict['sotip_ari'].append(sotip_ari)
                    
                    print(f'PCA:{_PCA},CellKnn:{_knn},CellRes:{leiden_res},MEKnn:{ME_knn},EMDNN:{n_neighbors},leidenAri:{leiden_ari},sotipAri:{sotip_ari}')
                    cur_save = f'{data_name},PCA:{_PCA},CellKnn:{_knn},CellRes:{leiden_res},MEKnn:{ME_knn},EMDNN:{n_neighbors}'
              
                    sc.pl.spatial(
                        adata,
                        color=['leiden_EMD_merge'],
                        title=f'PCA:{_PCA},CellKnn:{_knn},CellRes:{leiden_res},MEKnn:{ME_knn},EMDNN:{n_neighbors},leidenAri:{leiden_ari},sotipAri:{sotip_ari}',
                        save=cur_save,
                        show=False
                    )
                    adata.write_h5ad(f'search_{data_name}/{cur_save}.h5ad')
                del adata
                gc.collect()
        pd_df = pd.DataFrame(pd_dict)
        pd_df.to_csv(f'search_{data_name}/sotip_search_PCA{_PCA}_cellknn{_knn}_hvgFalse.csv')
        pd_dict = {
            'PCA':[],
            'Cell_cls_knn':[],
            'Cell_cls_res':[],
            'ME_knn':[],
            'EMD_neighbors':[],
            'leiden_ari':[],
            'sotip_ari':[]
        }
    
pd_df = pd.DataFrame(pd_dict)
pd_df.to_csv(f'search_{data_name}/sotip_search.csv')
                    
                    
