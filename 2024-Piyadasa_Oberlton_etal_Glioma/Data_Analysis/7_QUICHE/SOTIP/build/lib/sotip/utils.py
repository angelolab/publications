
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

def generate_spatial_mat(sz):
#     sz = 256
    spatial_mat = []
    for i in range(sz[0]):
        for j in range(sz[1]):
            spatial_mat.append([i,j])
    spatial_mat = np.array(spatial_mat)
    return spatial_mat



def plot_matrix(mat):
    fig, ax = plt.subplots(figsize=(10,10))  
    sns.heatmap(mat,annot=True,ax=ax)
    plt.gca().set_aspect('equal', adjustable='box') 




def plot_scanpy_spatial_equal(adata,color,size=10,figsize=(10,10),marker='.'):
    n_plots = len(color)
    fig,axes = plt.subplots(1,n_plots,figsize=figsize)
    
    for i in range(n_plots):
        sc.pl.embedding(adata,basis="spatial",color=[color[i]],show=False,size=size,ax=axes[i],marker=marker)

    
        axes[i].set_aspect('equal', adjustable='box') 

    fig.tight_layout()






def load_ST_file(file_fold, count_file='filtered_feature_bc_matrix.h5', load_images=True, file_Adj=None):
    adata_h5 = sc.read_visium(file_fold, load_images=load_images, count_file=count_file)
    adata_h5.var_names_make_unique()

    if load_images is False:
        if file_Adj is None:
            file_Adj = os.path.join(file_fold, "spatial/tissue_positions_list.csv")

        positions = pd.read_csv(file_Adj, header=None)
        positions.columns = [
            'barcode',
            'in_tissue',
            'array_row',
            'array_col',
            'pxl_col_in_fullres',
            'pxl_row_in_fullres',
        ]
        positions.index = positions['barcode']
        adata_h5.obs = adata_h5.obs.join(positions, how="left")
        adata_h5.obsm['spatial'] = adata_h5.obs[['pxl_row_in_fullres', 'pxl_col_in_fullres']].to_numpy()
        adata_h5.obs.drop(columns=['barcode', 'pxl_row_in_fullres', 'pxl_col_in_fullres'], inplace=True)

    print('adata: (' + str(adata_h5.shape[0]) + ', ' + str(adata_h5.shape[1]) + ')')
    return adata_h5



def heter_count_cls(adata):
    from toolz import compose
    cls_count = np.apply_along_axis(compose(len, np.unique), 1, adata.obsm['ME'])
    adata.obs['heter_nct'] = cls_count


def mask_adata_by_obs(adata,mask_list,masked_obs):
    # mask_list: cls list to be highlighted
    # masked_obs: obs to be masked
    # adata
    orig_cat = np.array(adata.obs[masked_obs].cat.categories).tolist()
    orig_color = adata.uns[f'{masked_obs}_colors'].copy()
    mask_array = np.array(adata.obs[masked_obs].copy())
    mask_array = [i if i in mask_list else '-1' for i in mask_array]
    adata.obs[f'{masked_obs}_masked'] =  mask_array
    adata.obs[f'{masked_obs}_masked'] = adata.obs[f'{masked_obs}_masked'].astype('category')

    left_cat = list(filter(lambda a:a in mask_list,orig_cat))
    left_cat_idx = [orig_cat.index(i) for i in left_cat]
    left_cat_color = [orig_color[i] for i in left_cat_idx]
    left_cat_color.append('k')
    left_cat.append('-1')
    adata.obs[f'{masked_obs}_masked'] = adata.obs[f'{masked_obs}_masked'].cat.set_categories(left_cat)
    adata.uns[f'{masked_obs}_masked_colors'] = left_cat_color


def mask_adata_by_list(adata,mask_list,masked_obs,mask_color='k'):
    # mask_list: cls list to be highlighted
    # masked_obs: obs to be masked
    # adata
    orig_cat = np.array(adata.obs[masked_obs].cat.categories).tolist()
    orig_color = adata.uns[f'{masked_obs}_colors'].copy()
    mask_array = np.array(adata.obs[masked_obs].copy())
    mask_array[np.setdiff1d(np.arange(adata.shape[0]),np.array(mask_list))]='-1'
    adata.obs[f'{masked_obs}_masked'] =  mask_array
    adata.obs[f'{masked_obs}_masked'] = adata.obs[f'{masked_obs}_masked'].astype('category')
    
    orig_cat.append('-1')
    adata.obs[f'{masked_obs}_masked'] = adata.obs[f'{masked_obs}_masked'].cat.set_categories(orig_cat)
    
    orig_color.append(mask_color)
    adata.uns[f'{masked_obs}_masked_colors'] = orig_color
    

def plot_ims(adata_pixel,protein,save=None,log=True,title=True,dpi=500):
    sz = int(np.sqrt(adata_pixel.shape[0]))
    img = np.array(adata_pixel[:,protein].X.reshape(sz,sz))

    if log: img = np.log1p(img)
    plt.imshow(img)
    plt.axis('off')
    if title:
        plt.title(protein)
    if save:
        plt.savefig(save,bbox_inches='tight',dpi=dpi)


def search_res(adata, target_num, start=0.4, step=0.1, tol=5e-3, max_run=10):
#example: res_recom = search_res(adata, target_num, start=1, step=0.1, tol=5e-3, max_run=10)
    res=start
    print("Start at res = ", res, "step = ", step)

    sc.tl.leiden(adata,resolution=start)
    y_pred = adata.obs['leiden'].copy()
    old_num=len(set(y_pred))
    print("Res = ", res, "Num of clusters = ", old_num)
    run=0
    while old_num!=target_num:
 
        old_sign=1 if (old_num<target_num) else -1
        sc.tl.leiden(adata,resolution=res+step*old_sign)
        
        y_pred = adata.obs['leiden'].copy()

        new_num=len(set(y_pred))
        print("Res = ", res+step*old_sign, "Num of clusters = ", new_num)
        if new_num==target_num:
            res=res+step*old_sign
            print("recommended res = ", str(res))
            return res
        new_sign=1 if (new_num<target_num) else -1
        if new_sign==old_sign:
            res=res+step*old_sign
            print("Res changed to", res)
            old_num=new_num
        else:
            step=step/2
            print("Step changed to", step)
        if run >max_run:
            print("Exact resolution not found")
            print("Recommended res = ", str(res))
            return res
        run+=1
    print("recommended res = ", str(res))
    return res

