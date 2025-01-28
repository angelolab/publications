import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
import igraph as ig
import scipy
from scipy.sparse import *
from sklearn.neighbors import NearestNeighbors
from sklearn.decomposition import PCA
from itertools import combinations
import multiprocessing as mp
from tqdm import tqdm
from functools import partial
from scipy.spatial.distance import pdist, squareform
import phate
import timeit
import anndata
import logging
import os
import xarray as xr
from skimage.measure import label, regionprops
from skimage import measure
from skimage.draw import line
from skimage.segmentation import find_boundaries
from PIL import Image
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from tqdm import tqdm
import squidpy as sq
from scipy.sparse import coo_matrix
import warnings
import os
import pathlib
from typing import Dict, List, Literal, Optional, Tuple, Union
from matplotlib import gridspec
from matplotlib.axes import Axes

import matplotlib.cm as cm
import matplotlib.colors as colors
from matplotlib import colormaps, patches
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import natsort
import numpy as np
import pandas as pd
from pandas.core.groupby.generic import DataFrameGroupBy
import skimage
import xarray as xr
from mpl_toolkits.axes_grid1 import make_axes_locatable
from skimage import io

import numba as nb
import itertools
import os
import pathlib
import re
import natsort as ns
import numpy as np
import pandas as pd
import skimage.io as io
from alpineer import image_utils, io_utils, misc_utils
from tqdm.notebook import tqdm_notebook as tqdm
import xarray as xr
from ark import settings
from skimage.segmentation import find_boundaries
from skimage.util import img_as_ubyte
from tqdm.auto import tqdm
from ark import settings
from skimage.segmentation import find_boundaries
from ark.utils.data_utils import (
    erode_mask,
    generate_cluster_mask,
    save_fov_mask,
    map_segmentation_labels,
)


def get_igraph(W = None,
               directed: bool = None):
    """Converts adjacency matrix into igraph object

    Parameters
    W: (default = None)
        adjacency matrix
    directed: bool (default = None)
        whether graph is directed or not
    ----------

    Returns
    g: ig.Graph
        graph of adjacency matrix
    ----------
    """
    sources, targets = W.nonzero()
    weights = W[sources, targets]
    if type(weights) == np.matrix:
        weights = weights.A1 #flattens 
    g = ig.Graph(directed = directed)
    g.add_vertices(np.shape(W)[0])
    g.add_edges(list(zip(sources, targets)))
    g.es['weight'] = weights  
    return g

def heat_kernel(dist = None,
                radius = 3):
    """Transforms distances into weights using heat kernel
    Parameters
    dist: np.ndarray (default = None)
        distance matrix (dimensions = cells x k)
    radius: np.int (default = 3)
        defines the per-cell bandwidth parameter (distance to the radius nn)
    ----------
    Returns
    s: np.ndarray
        array containing between cell similarity (dimensions = cells x k)
    ----------
    """         
    sigma = dist[:, [radius]]  # per cell bandwidth parameter (distance to the radius nn)
    s = np.exp(-1 * (dist**2)/ (2.*sigma**2)) # -||x_i - x_j||^2 / 2*sigma_i**2
    return s

def construct_affinity(X = None,
                        k: int = 10,
                        radius: int = 3,
                        n_pcs = None,
                        random_state: int = 0, 
                        n_jobs: int = -1):
    """Computes between cell affinity knn graph using heat kernel
    Parameters
    X: np.ndarray (default = None)
        Data (dimensions = cells x features)
    k: int (default = None)
        Number of nearest neighbors
    radius: int (default = 3)
        Neighbor to compute per cell distance for heat kernel bandwidth parameter
    n_pcs: int (default = None)
        number of principal components to compute pairwise Euclidean distances for between-cell affinity graph construction. If None, uses adata.X
    n_jobs: int (default = -1)
        Number of tasks  
    ----------
    Returns
    W: np.ndarray
        sparse symmetric matrix containing between cell similarity (dimensions = cells x cells)
    ----------
    """
    if n_pcs is not None:
        n_comp = min(n_pcs, X.shape[1])
        pca_op = PCA(n_components=n_comp, random_state = random_state)
        X_ = pca_op.fit_transform(X)
    else:
        X_ = X.copy()

    # find kNN
    knn_tree = NearestNeighbors(n_neighbors=k, algorithm='ball_tree', metric='euclidean', n_jobs=n_jobs).fit(X_)
    dist, nn = knn_tree.kneighbors()  # dist = cells x knn (no self interactions)

    # transform distances using heat kernel
    s = heat_kernel(dist, radius = radius) # -||x_i - x_j||^2 / 2*sigma_i**2
    rows = np.repeat(np.arange(X.shape[0]), k)
    cols = nn.reshape(-1)
    W = scipy.sparse.csr_matrix((s.reshape(-1), (rows, cols)), shape=(X.shape[0], X.shape[0]))

    # make symmetric
    bigger = W.transpose() > W
    W = W - W.multiply(bigger) + W.transpose().multiply(bigger)
    return W

def compute_ground_dist(adata, expression_mat, knn, cluster_key, clusters):
    #euclidean distance between cell subtypes in phate space 
    #first define cell centroids as the mean in PHATE
    #then define euclidean distance
    op = phate.PHATE(knn = knn)
    X_phate = op.fit_transform(expression_mat)
    phate_df = pd.DataFrame(X_phate, index = adata.obs[cluster_key])
    centroids = phate_df.groupby(cluster_key).mean().loc[clusters]
    ground_dist = squareform(pdist(centroids))
    return ground_dist

# def emd_samples(pair_list, ground_dist = None, emd_solver = 'pyemd'):
#     idx_pair, data_pair = pair_list
#     first_prop = data_pair[0]
#     second_prop = data_pair[1]
#     if emd_solver == 'pyemd':
#         emd_dist = pyemd.emd(first_prop, second_prop, ground_dist)
#     elif emd_solver == 'ot_emd':
#         emd_dist = ot.emd2(first_prop,second_prop,ground_dist)
#     elif emd_solver == 'ot_sinkhorn':
#         emd_dist = ot.sinkhorn2(first_prop, second_prop, ground_dist, 0.01)
#     return idx_pair[0], idx_pair[1], emd_dist

# def compute_EMD(prop_arr, ground_dist, n_jobs, emd_solver):
#     if n_jobs == -1:
#         n_jobs = mp.cpu_count()
#     elif n_jobs < -1:
#         n_jobs = mp.cpu_count() + 1 + n_jobs
        
#     n_samples = prop_arr.shape[0]
#     total_pairwise = n_samples * (n_samples - 1) // 2
#     pair_list = _gen_pair_list(n_samples, prop_arr)

#     p = mp.Pool(n_jobs)

#     emd_dist = np.zeros(shape=(n_samples, n_samples))
#     for result in tqdm(p.imap(partial(emd_samples, ground_dist = ground_dist, emd_solver = emd_solver), pair_list), total = total_pairwise, desc="Computing EMD"):
#         emd_dist[result[0], result[1]] = result[2]

#     return emd_dist

# def _gen_pair_list(n_samples, prop_arr):
#     idx_list = list(range(n_samples))
#     idx_pairs = combinations(idx_list, r=2)
#     data_pairs = combinations(prop_arr, r=2)
#     pair_list = zip(idx_pairs, data_pairs)
#     return pair_list

# def test_EMD(prop_arr, ground_dist, n_jobs, emd_solver):
#     if n_jobs == -1:
#         n_jobs = mp.cpu_count()
#     elif n_jobs < -1:
#         n_jobs = mp.cpu_count() + 1 + n_jobs
    
#     n_samples = prop_arr.shape[0]
#     total_pairwise = n_samples * (n_samples - 1) // 2
#     pair_list = _gen_pair_list(n_samples, prop_arr)

#     p = mp.Pool(n_jobs)

#     emd_dist = np.zeros(shape=(n_samples, n_samples))
    
#     start_time = timeit.default_timer()
#     for result in tqdm(p.imap(partial(emd_samples, ground_dist=ground_dist, emd_solver = emd_solver), pair_list), total=total_pairwise, desc="Computing EMD"):
#         emd_dist[result[0], result[1]] = result[2]
    
#     elapsed_time = timeit.default_timer() - start_time
    
#     return elapsed_time

# def compute_EMD_brute(prop_arr, ground_dist_mat):
#     n_niches = prop_arr.shape[0]
#     emd_dist = np.zeros(shape=(n_niches,n_niches))
#     for i in range(n_niches):
#         for j in range(i, n_niches):
#             first_prop = np.ascontiguousarray(prop_arr[i,:])
#             second_prop = np.ascontiguousarray(prop_arr[j,:])
#             emd_dist_ = pyemd.emd(first_prop,second_prop, ground_dist_mat)
#             emd_dist[i,j] = emd_dist_
#             emd_dist[j,i] = emd_dist_
#     return emd_dist

def spatial_neighbors_fov(fov, adata, radius, p, labels_key, spatial_key, fov_key):
    adata_fov = adata[np.where(adata.obs[fov_key] == fov)[0], :].copy()
    spatial_kdTree = cKDTree(adata_fov.obsm[spatial_key])
    nn = spatial_kdTree.query_ball_point(adata_fov.obsm[spatial_key], r=radius, p=p)  # no self interactions

    prop_df = pd.DataFrame()
    for i in range(0, nn.shape[0]):
        counts = adata_fov.obs[labels_key].iloc[nn[i]].value_counts().copy()
        prop_df_ = pd.DataFrame(counts / counts.sum()).transpose()
        prop_df = pd.concat([prop_df, prop_df_], axis=0)
    prop_df.index = adata_fov.obs_names

    return prop_df, nn

def compute_spatial_neighbors_parallel(adata, radius=200, p=2, min_cell_threshold=5, labels_key='cell_cluster', spatial_key='X_spatial', fov_key='fov', n_jobs = -1):
    if n_jobs == -1:
        n_jobs = mp.cpu_count()
    elif n_jobs < -1:
        n_jobs = mp.cpu_count() + 1 + n_jobs

    p = mp.Pool(n_jobs)

    fovs = adata.obs[fov_key].unique().tolist() #should preserve original order
    n_fovs = len(fovs)

    prop_df = []
    nn_list = []
    for result in tqdm(p.imap(partial(spatial_neighbors_fov, adata=adata, radius=radius, p = p, labels_key=labels_key, spatial_key=spatial_key, fov_key=fov_key), fovs), total=n_fovs, desc='computing spatial niches'):
        prop_df.append(result[0])
        nn_list.append(result[1])

    prop_df = pd.concat(prop_df, axis=0)
    prop_df = prop_df.fillna(0).copy()
    prop_df = prop_df.loc[adata.obs_names]

    return prop_df, nn_list

def spatial_niches_radius(adata,
                        radius = 200,
                        p = 2,
                        min_cell_threshold= 0,
                        labels_key = 'cell_cluster',
                        spatial_key = 'X_spatial',
                        fov_key = 'fov',
                        n_jobs = -1):
    niche_df = []
    cells2remove = []
    nn_dict = {}
    for fov in np.unique(adata.obs[fov_key]):
        adata_fov = adata[np.where(adata.obs[fov_key] == fov)[0], :].copy()
        spatial_kdTree = cKDTree(adata_fov.obsm[spatial_key])
        nn = spatial_kdTree.query_ball_point(adata_fov.obsm[spatial_key], r=radius, p = p, workers = n_jobs) #no self interactions
        nn_dict[fov] = nn
        niche_df_fov = pd.DataFrame()
        for i in range(0, nn.shape[0]):
            if len(nn[i]) < min_cell_threshold:
                cells2remove.append(adata_fov.obs_names[i])
            else:
                counts = adata_fov.obs[labels_key].iloc[nn[i]].value_counts().copy()
                niche_df_fov_ = pd.DataFrame(counts / counts.sum()).transpose()
                niche_df_fov = pd.concat([niche_df_fov, niche_df_fov_], axis = 0)
        niche_df_fov.index = adata_fov.obs_names[~np.isin(adata_fov.obs_names, cells2remove)]
        niche_df.append(niche_df_fov)
    try:
        cells2remove = np.concatenate(cells2remove)    
    except:
        pass
    niche_df = pd.concat(niche_df, axis = 0)
    niche_df = niche_df.fillna(0).copy()
    niche_df = niche_df.loc[adata.obs_names[~np.isin(adata.obs_names, list(set(cells2remove)))]]
    return niche_df, nn_dict

def spatial_niches_khop(adata,
                        radius = 200,
                        p = 2,
                        k = 10,
                        khop = 3, 
                        labels_key = 'cell_cluster',
                        spatial_key = 'X_spatial',
                        fov_key = 'fov',
                        min_cell_threshold= 0,
                        n_jobs = -1):
    niche_df = []
    cells2remove = []
    nn_dict = {}
    for fov in np.unique(adata.obs[fov_key]):
        adata_fov = adata[np.where(adata.obs[fov_key] == fov)[0], :].copy()
        if adata_fov.shape[0] < k*khop: #not enough neighbors
            logging.info(f'{fov} is less than khop neighborhood size {k*khop}. Removing FOV.')
            cells2remove.append(list(adata_fov.obs_names))
        else:
            niche_df_fov = pd.DataFrame()
            knn_kdtree = NearestNeighbors(n_neighbors= khop*k, p = p, algorithm='kd_tree', metric='euclidean', n_jobs=n_jobs).fit(adata_fov.obsm[spatial_key])
            dist, nn = knn_kdtree.kneighbors()
            nn_dict[fov] = nn
            for i in range(0, nn.shape[0]):
                niche_nn = nn[i, :][dist[i, :] < radius]
                if len(niche_nn) < min_cell_threshold:
                    cells2remove.append(adata_fov.obs_names[i])
                else:
                    counts = adata_fov.obs[labels_key].iloc[niche_nn].value_counts().copy()
                    niche_df_fov_ = pd.DataFrame(counts / counts.sum()).transpose()
                    niche_df_fov = pd.concat([niche_df_fov, niche_df_fov_], axis = 0)
            niche_df_fov.index = adata_fov.obs_names[~np.isin(adata_fov.obs_names, cells2remove)]
            niche_df.append(niche_df_fov)
    try:
        cells2remove = np.concatenate(cells2remove)    
    except:
        pass
    niche_df = pd.concat(niche_df, axis = 0)
    niche_df = niche_df.fillna(0).copy()
    niche_df = niche_df.loc[adata.obs_names[~np.isin(adata.obs_names, list(set(cells2remove)))]]
    return niche_df, nn_dict

def compute_spatial_niches(adata,
                            radius = 200,
                            p = 2,
                            k = 10,
                            khop = 3, 
                            labels_key = 'cell_cluster',
                            spatial_key = 'X_spatial',
                            fov_key = 'fov',
                            min_cell_threshold = 0, 
                            n_jobs = -1):
    
    if khop is not None:
        niche_df, nn_dict = spatial_niches_khop(adata, radius = radius, p = p, k = k, khop = khop, min_cell_threshold = min_cell_threshold, labels_key = labels_key, spatial_key = spatial_key, fov_key = fov_key, n_jobs = n_jobs)
    else:
        niche_df, nn_dict = spatial_niches_radius(adata, radius = radius, p = p, labels_key = labels_key, min_cell_threshold = min_cell_threshold, spatial_key = spatial_key, fov_key = fov_key, n_jobs = n_jobs)
    adata_niche = anndata.AnnData(niche_df)
    adata_niche.obs = adata.obs.loc[niche_df.index, :]
    return adata_niche, nn_dict

def filter_fovs(adata, patient_key, threshold):
    n_niches = adata.obs[patient_key].value_counts(sort=False)
    adata = adata[~np.isin(adata.obs[patient_key], n_niches[n_niches < threshold].index)]    
    return adata

def construct_niche_similarity_graph(adata, k = 100, n_jobs = -1):
    nn = NearestNeighbors(n_neighbors = k, algorithm = 'kd_tree', n_jobs = n_jobs).fit(adata.X)
    connectivities = nn.kneighbors_graph(mode = 'connectivity')
    distances = nn.kneighbors_graph(mode = 'distance')            
    adata.obsp['connectivities'] = connectivities
    adata.obsp['distances'] = distances
    adata.uns["nhood_neighbors_key"] = None
    adata.uns['neighbors'] = {'connectivities_key': 'connectivities',
                            'distances_key': 'distances',
                            'params': {'n_neighbors': k,
                                    'method': 'umap',
                                    'distance': 'Euclidean'}}
    return adata

def build_milo_graph(mdata, feature_key = 'spatial_nhood'):
    #define idxs as all cells
    mdata[feature_key].obs['nhood_ixs_random'] = 1 
    mdata[feature_key].obs['nhood_ixs_refined'] = 1
    #binary knn graph
    knn_graph = mdata[feature_key].obsp['connectivities'].copy()
    knn_graph[knn_graph!=0] =1 #binarize
    mdata[feature_key].obsm['nhoods'] = knn_graph.copy()
    #distance matrix
    knn_dists = mdata[feature_key].obsp["distances"]
    knn_dists = knn_dists.max(1).toarray().ravel()
    mdata[feature_key].obs["nhood_kth_distance"] = knn_dists
    return mdata

def build_nhood_graph(mdata, basis = "X_umap", feature_key = "spatial_nhood"):
    """Build graph of neighbourhoods used for visualization of DA results

    Args:
        mdata: MuData object
        basis: Name of the obsm basis to use for layout of neighbourhoods (key in `adata.obsm`). Defaults to "X_umap".
        feature_key: If input data is MuData, specify key to cell-level AnnData object. Defaults to 'spatial_nhood'.

    Returns:
        - `milo_mdata['milo'].varp['nhood_connectivities']`: graph of overlap between neighbourhoods (i.e. no of shared cells)
        - `milo_mdata['milo'].var["Nhood_size"]`: number of cells in neighbourhoods
    """
    adata = mdata[feature_key]
    # # Add embedding positions
    mdata["milo"].varm["X_milo_graph"] = adata[adata.obs["nhood_ixs_refined"] == 1].obsm[basis]
    # Add nhood size
    mdata["milo"].var["Nhood_size"] = np.array(adata.obsm["nhoods"].sum(0)).flatten()
    # Add adjacency graph
    mdata["milo"].varp["nhood_connectivities"] = adata.obsm["nhoods"].T.dot(adata.obsm["nhoods"])
    mdata["milo"].varp["nhood_connectivities"].setdiag(0)
    mdata["milo"].varp["nhood_connectivities"].eliminate_zeros()
    mdata["milo"].uns["nhood"] = {
        "connectivities_key": "nhood_connectivities",
        "distances_key": "",
    }

def save_colored_mask(
    fov: str,
    save_dir: str,
    suffix: str,
    data: np.ndarray,
    cmap: colors.ListedColormap,
    norm: colors.BoundaryNorm,
) -> None:
    """Saves the colored mask to the provided save directory.

    Args:
        fov (str):
            The name of the FOV.
        save_dir (str):
            The directory where the colored mask will be saved.
        suffix (str):
            The suffix to append to the FOV name.
        data (np.ndarray):
            The mask to save.
        cmap (colors.ListedColormap):
            The colormap to use for the mask.
        norm (colors.BoundaryNorm):
            The normalization to use for the mask.
    """

    # Create the save directory if it does not exist
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)

    # Create the colored mask
    colored_mask = img_as_ubyte(cmap(norm(data)))

    # Save the image
    image_utils.save_image(
        fname=os.path.join(save_dir, f"{fov}{suffix}"),
        data=colored_mask,
    )

class ClusterMaskData:
    """
    A class containing the cell labels, cluster labels, and segmentation labels for the
    whole cohort. Also contains the mapping from the segmentation label to the cluster
    label for each FOV.
    """

    fov_column: str
    label_column: str
    cluster_column: str
    unique_fovs: List[str]
    cluster_id_column: str
    unassigned_id: int
    n_clusters: int
    mapping: pd.DataFrame

    def __init__(
        self, data: pd.DataFrame, fov_col: str, label_col: str, cluster_col: str
    ) -> None:
        """
        A class containing the cell data, cell label column, cluster column and the mapping from a
        cell label to a cluster.

        Args:
            data (pd.DataFrame):
                A cell table with the cell label column and the cluster column.
            fov_col (str):
                The name of the column in the cell table that contains the FOV ID.
            label_col (str):
                The name of the column in the cell table that contains the cell label.
            cluster_col (str):
                The name of the column in the cell table that contains the cluster label.
        """
        self.fov_column: str = fov_col
        self.label_column: str = label_col
        self.cluster_column: str = cluster_col
        self.cluster_id_column: str = "cluster_id"

        # Extract only the necessary columns: fov ID, segmentation label, cluster label
        mapping_data: pd.DataFrame = data[
            [self.fov_column, self.label_column, self.cluster_column]
        ].copy()

        # Add a cluster_id_column to the column in case the cluster_column is
        # non-numeric (i.e. string), index in ascending order of cell_meta_cluster
        cluster_name_id = pd.DataFrame(
            {self.cluster_column: list(mapping_data[self.cluster_column].cat.categories)})
        
        cluster_name_id.sort_values(by=f'{self.cluster_column}', inplace=True)
        cluster_name_id.reset_index(drop=True, inplace=True)

        cluster_name_id[self.cluster_id_column] = (cluster_name_id.index + 1).astype(np.int32)

        self.cluster_name_id = cluster_name_id

        # merge the cluster_id_column to the mapping_data dataframe
        mapping_data = mapping_data.merge(right=self.cluster_name_id, on=self.cluster_column)

        mapping_data = mapping_data.astype(
            {
                self.fov_column: str,
                self.label_column: np.int32,
                self.cluster_id_column: np.int32,
            }
        )
        self.unique_fovs: List[str] = ns.natsorted(
            mapping_data[self.fov_column].unique().tolist()
        )

        self.unassigned_id: np.int32 = np.int32(
            mapping_data[self.cluster_id_column].max() + 1
        )

        self.n_clusters: int = mapping_data[self.cluster_id_column].max()

        # For each FOV map the segmentation label 0 (background) to the cluster label 0
        cluster0_mapping: pd.DataFrame = pd.DataFrame(
            data={
                self.fov_column: self.unique_fovs,
                self.label_column: np.repeat(0, repeats=len(self.unique_fovs)),
                self.cluster_column: np.repeat(0, repeats=len(self.unique_fovs)),
                self.cluster_id_column: np.repeat(0, repeats=len(self.unique_fovs)),
            }
        )

        mapping_data = pd.concat(objs=[mapping_data, cluster0_mapping]).astype(
            {
                self.fov_column: str,
                self.label_column: np.int32,
                self.cluster_id_column: np.int32,
            }
        )

        # Sort by FOV first, then by segmentation label
        self.mapping = mapping_data.sort_values(by=[self.fov_column, self.label_column])

    def fov_mapping(self, fov: str) -> pd.DataFrame:
        """Returns the mapping for a specific FOV.
        Args:
            fov (str):
                The FOV to get the mapping for.
        Returns:
            pd.DataFrame:
                The mapping for the FOV.
        """
        misc_utils.verify_in_list(requested_fov=[fov], all_fovs=self.unique_fovs)
        fov_data: pd.DataFrame = self.mapping[self.mapping[self.fov_column] == fov]
        return fov_data.reset_index(drop=True)

    @property
    def cluster_names(self) -> List[str]:
        """Returns the cluster names.
        Returns:
            List[str]:
                The cluster names.
        """
        return self.cluster_name_id[self.cluster_column].tolist()





def cohort_cluster_plot(
    fovs: List[str],
    seg_dir: Union[pathlib.Path, str],
    save_dir: Union[pathlib.Path, str],
    cell_data: pd.DataFrame,
    fov_col: str = settings.FOV_ID,
    label_col: str = settings.CELL_LABEL,
    cluster_col: str = settings.CELL_TYPE,
    seg_suffix: str = "_whole_cell.tiff",
    cmap: Union[str, pd.DataFrame] = "viridis",
    style: str = "seaborn-v0_8-paper",
    erode: bool = False,
    display_fig: bool = False,
    fig_file_type: str = "pdf",
    figsize: tuple = (10, 10),
    dpi: int = 600,
) -> None:
    """
    Saves the cluster masks for each FOV in the cohort as the following:
    - Cluster mask numbered 1-N, where N is the number of clusters (tiff)
    - Cluster mask colored by cluster with or without a colorbar (png)
    - Cluster mask colored by cluster (tiff).

    Args:
        fovs (List[str]): A list of FOVs to generate cluster masks for.
        seg_dir (Union[pathlib.Path, str]): The directory containing the segmentation masks.
        save_dir (Union[pathlib.Path, str]): The directory to save the cluster masks to.
        cell_data (pd.DataFrame): The cell data table containing the cluster labels.
        fov_col (str, optional): The column containing the FOV name. Defaults to settings.FOV_ID.
        label_col (str, optional): The column containing the segmentaiton label.
            Defaults to settings.CELL_LABEL.
        cluster_col (str, optional): The column containing the cluster a segmentation label
            belongs to. Defaults to settings.CELL_TYPE.
        seg_suffix (str, optional): The kind of segmentation file to read.
            Defaults to "_whole_cell.tiff".
        cmap (str, pd.DataFrame, optional): The colormap to generate clusters from,
            or a DataFrame, where the user can specify their own colors per cluster.
            The color column must be labeled "color". Defaults to "viridis".
        style (str, optional): Set the matplotlib style image style. Defaults to 
            "seaborn-v0_8-paper".
            View the available styles here: 
            https://matplotlib.org/stable/gallery/style_sheets/style_sheets_reference.html
            Or run matplotlib.pyplot.style.available in a notebook to view all the styles.
        erode (bool, optional): Option to "thicken" the cell boundary via the segmentation label
            for visualization purposes. Defaults to False.
        display_fig (bool, optional): Option to display the cluster mask plots as they are
            generated. Defaults to False. Displaying each figure can use a lot of memory,
            so it's best to try to visualize just a few FOVs, before generating the cluster masks
            for the entire cohort.
        fig_file_type (str, optional): The file type to save figures as. Defaults to 'png'.
        figsize (tuple, optional):
            The size of the figure to display. Defaults to (10, 10).
        dpi (int, optional):
            The resolution of the image to use for saving. Defaults to 300.
    """

    plt.style.use(style)

    if isinstance(seg_dir, str):
        seg_dir = pathlib.Path(seg_dir)

    try:
        io_utils.validate_paths(seg_dir)
    except ValueError:
        raise ValueError(f"Could not find the segmentation directory at {seg_dir.as_posix()}")

    if isinstance(save_dir, str):
        save_dir = pathlib.Path(save_dir)
        if not save_dir.exists():
            save_dir.mkdir(parents=True, exist_ok=True)
    if isinstance(fovs, str):
        fovs = [fovs]

    # Create the subdirectories for the 3 cluster mask files
    for sub_dir in ["cluster_masks", "cluster_masks_colored", "cluster_plots"]:
        (save_dir / sub_dir).mkdir(parents=True, exist_ok=True)

    cmd = ClusterMaskData(
        data=cell_data,
        fov_col=fov_col,
        label_col=label_col,
        cluster_col=cluster_col,
    )
    if isinstance(cmap, pd.DataFrame):
        unique_clusters: pd.DataFrame = cmd.mapping[[cmd.cluster_column,
                                                     cmd.cluster_id_column]].drop_duplicates()
        
        # print(cmap.loc[:, cmd.cluster_column].values)
        # print(pd.Categorical(cmd.mapping[cmd.cluster_column], categories = cmap[cmd.cluster_column])
        # unique_clusters = pd.DataFrame(cmap.loc[:, cmd.cluster_column].values, columns = [cmd.cluster_column])
        # unique_clusters['cluster_id'] = np.arange(0, len(unique_clusters))

        # print(cmd.cluster_names)

        # unique_clusters[cmd.cluster_column] = cmap[cmd.cluster_column]
        
        # print(cmap)
        
        cmap_colors: pd.DataFrame = cmap.merge(
            right=unique_clusters,
            on=cmd.cluster_column
        ).sort_values(by="cluster_id")["color"].values

        print(len(cmap_colors))

        # print(cmap_colors)


        # cmap_colors: pd.DataFrame = cmap.merge(
        #     right=unique_clusters,
        #     on=cmd.cluster_column
        # )["color"].values

        colors_like: list[bool] = [colors.is_color_like(c) for c in cmap_colors]

        if not all(colors_like):
            bad_color_values: np.ndarray = cmap_colors[~np.array(colors_like)]
            raise ValueError(
                ("Not all colors in the provided cmap are valid colors."
                 f"The following colors are invalid: {bad_color_values}"))

        np_colors = colors.to_rgba_array(cmap_colors)

        color_map, norm = create_cmap(np_colors, n_clusters=cmd.n_clusters)

    if isinstance(cmap, str):
        color_map, norm = create_cmap(cmap, n_clusters=cmd.n_clusters)

    # create the pixel cluster masks across each fov
    with tqdm(total=len(fovs), desc="Cluster Mask Generation", unit="FOVs") as pbar:
        for fov in fovs:
            pbar.set_postfix(FOV=fov)

            # generate the cell mask for the FOV
            cluster_mask: np.ndarray = generate_cluster_mask(
                fov=fov,
                seg_dir=seg_dir,
                cmd=cmd,
                seg_suffix=seg_suffix,
                erode=erode,
            )

            # save the cluster mask generated
            save_fov_mask(
                fov,
                data_dir=save_dir / "cluster_masks",
                mask_data=cluster_mask,
                sub_dir=None,
            )

            save_colored_mask(
                fov=fov,
                save_dir=save_dir / "cluster_masks_colored",
                suffix=".tiff",
                data=cluster_mask,
                cmap=color_map,
                norm=norm,
            )
            # print(cmd.cluster_names)
            # print(pd.DataFrame(index = cmd.cluster_names).loc[cmap[cmd.cluster_column].values])
            cluster_labels = ["Background"] + cmd.cluster_names + ["Unassigned"]

            fig = plot_cluster(
                image=cluster_mask,
                fov=fov,
                cmap=color_map,
                norm=norm,
                cbar_visible=True,
                cbar_labels=cluster_labels,
                figsize=figsize,
                dpi=dpi,
            )

            fig.savefig(
                fname=os.path.join(save_dir, "cluster_plots", f"{fov}.{fig_file_type}"),
            )

            if display_fig:
                fig.show(warn=False)
            else:
                plt.close(fig)

            pbar.update(1)


def create_cmap(cmap: Union[np.ndarray, list[str], str],
                n_clusters: int) -> tuple[colors.ListedColormap, colors.BoundaryNorm]:
    """
    Creates a discrete colormap and a boundary norm from the provided colors.

    Args:
        cmap (Union[np.ndarray, list[str], str]): The colormap, or set of colors to use.
        n_clusters (int): The numbe rof clusters for the colormap.

    Returns:
        tuple[colors.ListedColormap, colors.BoundaryNorm]:
            The generated colormap and boundary norm.
    """

    """Creates a colormap and a boundary norm from the provided colors.

    Colors can be of any format that matplotlib accepts.
    See here for color formats: https://matplotlib.org/stable/tutorials/colors/colors.html


    Args:
        colors_array (): The colors to use for the colormap.

    Returns:
        tuple[colors.ListedColormap, colors.BoundaryNorm]: The colormap and the boundary norm
    """

    if isinstance(cmap, np.ndarray):
        if cmap.ndim != 2:
            raise ValueError(
                f"colors_array must be a 2D array, got {cmap.ndim}D array")
        if cmap.shape[0] != n_clusters:
            raise ValueError(
                f"colors_array must have {n_clusters} colors, got {cmap.shape[0]} colors")
        color_map = colors.ListedColormap(colors=_cmap_add_background_unassigned(cmap))
    if isinstance(cmap, list):
        if len(cmap) != n_clusters:
            raise ValueError(
                f"colors_array must have {n_clusters} colors, got {len(cmap)} colors")
    if isinstance(cmap, str):
        try:
            # colorcet colormaps are also supported
            # cmocean colormaps are also supported
            color_map = colormaps[cmap]
        except KeyError:
            raise KeyError(f"Colormap {cmap} not found.")
        colors_rgba: np.ndarray = color_map(np.linspace(0, 1, n_clusters))
        color_map: colors.ListedColormap = colors.ListedColormap(
            colors=_cmap_add_background_unassigned(colors_rgba))

    bounds = [i-0.5 for i in np.linspace(0, color_map.N, color_map.N + 1)]

    norm = colors.BoundaryNorm(bounds, color_map.N)
    return color_map, norm



def plot_cluster(
        image: np.ndarray,
        fov: str,
        cmap: colors.ListedColormap,
        norm: colors.BoundaryNorm,
        cbar_visible: bool = True,
        cbar_labels: list[str] = None,
        dpi: int = 300,
        figsize: tuple[int, int] = None) -> Figure:
    """
    Plots the cluster image with the provided colormap and norm.

    Args:
        image (np.ndarray):
            The cluster image to plot.
        fov (str):
            The name of the clustered FOV.
        cmap (colors.ListedColormap):
            A colormap to use for the cluster image.
        norm (colors.BoundaryNorm):
            A normalization to use for the cluster image.
        cbar_visible (bool, optional):
            Whether or not to display the colorbar. Defaults to True.
        cbar_labels (list[str], optional):
            Colorbar labels for the clusters. Devaults to None, where
            the labels will be automatically generated.
        dpi (int, optional):
            The resolution of the image to use for saving. Defaults to 300.
        figsize (tuple, optional):
            The size of the image to display. Defaults to (10, 10).

    Returns:
        Figure: Returns the cluster image as a matplotlib Figure.
    """
    # Default colorbar labels
    if cbar_labels is None:
        cbar_labels = [f"Cluster {x}" for x in range(1, len(cmap.colors))]

    fig: Figure = plt.figure(figsize=figsize, dpi=dpi)
    fig.set_layout_engine(layout="tight")
    gs = gridspec.GridSpec(nrows=1, ncols=1, figure=fig)
    fig.suptitle(f"{fov}")

    # Image axis
    ax: Axes = fig.add_subplot(gs[0, 0])
    ax.axis("off")
    ax.grid(visible=False)

    ax.imshow(
        X=image,
        cmap=cmap,
        norm=norm,
        origin="upper",
        aspect="equal",
        interpolation="none",
    )

    if cbar_visible:
        # # Manually set the colorbar
        divider = make_axes_locatable(fig.gca())
        cax = divider.append_axes(position="right", size="5%", pad="3%")

        cbar = fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap),
                            cax=cax, orientation="vertical", use_gridspec=True, pad=0.1,
                            shrink=0.9, drawedges=True)
        cbar.ax.set_yticks(
            ticks=np.arange(len(cbar_labels)),
            labels=cbar_labels
        )
        cbar.minorticks_off()

    return fig

def _cmap_add_background_unassigned(cluster_colors: np.ndarray):
    # A pixel with no associated metacluster (gray, #5A5A5A)
    unassigned_color: np.ndarray = np.array([0.0, 0.0, 0.0, 1.0])

    # A pixel assigned to the background (black, #000000)
    background_color: np.ndarray = np.array([0.0, 0.0, 0.0, 1.0])

    return np.vstack([background_color, cluster_colors, unassigned_color])


# def bound_radius(adata,
#                 distances_key = 'spatial_distances',
#                 connectivities_key = 'spatial_connectivities',
#                 radius = 200):

#     dist_mat = adata.obsp[distances_key].todense()
#     connect_mat = adata.obsp[connectivities_key].todense()
#     bool_idx = np.array(dist_mat >= radius)

#     connect_mat[bool_idx] = 0
#     connect_mat = scipy.sparse.csr_matrix(connect_mat)

#     adata.obsp[connectivities_key] = connect_mat
#     return adata

def set_elements_to_zero_lil(sparse_array, mask):
    sparse_array = sparse_array.tolil()
    mask = mask.tolil()
 
    nz = mask.nonzero()
    sparse_array[nz] = 0
 
    return sparse_array.tocsr()

def bound_radius(adata,
                distances_key = 'spatial_distances',
                connectivities_key = 'spatial_connectivities',
                radius = 200):

    dist_mat = adata.obsp[distances_key]
    connect_mat = adata.obsp[connectivities_key]
    bool_idx = dist_mat >= radius

    connect_mat = set_elements_to_zero_lil(connect_mat, bool_idx)
    adata.obsp[connectivities_key] = connect_mat
    return adata


def compute_spatial_neighbors(adata,
                            radius = 200,
                            n_neighbors = None,
                            spatial_key = 'X_spatial',
                            delaunay = True,
                            fov_key = 'fov',
                            coord_type = 'generic'):

    adata.obs[fov_key] = pd.Categorical(adata.obs[fov_key])
    if n_neighbors is not None:
        sq.gr.spatial_neighbors(adata, spatial_key = spatial_key, library_key = fov_key, n_neighs = n_neighbors, coord_type = coord_type)
        adata = bound_radius(adata, distances_key = 'spatial_distances', connectivities_key = 'spatial_connectivities', radius = radius)
    elif delaunay == True:
        sq.gr.spatial_neighbors(adata, spatial_key = spatial_key, library_key = fov_key, delaunay = delaunay, coord_type = coord_type)
        adata = bound_radius(adata, distances_key = 'spatial_distances', connectivities_key = 'spatial_connectivities', radius = radius)
    else:
        sq.gr.spatial_neighbors(adata, spatial_key = spatial_key, library_key = fov_key, radius = radius, coord_type = coord_type)
    
    return adata

def compute_niche_composition(adata,
                            connectivities_key = 'spatial_connectivities',
                            labels_key = 'mask_name',
                            min_cells = 3):

    count_list = []
    for i, name in tqdm(enumerate(adata.obs_names), desc="computing niche composition"):
        row, col = adata.obsp[connectivities_key][i, :].nonzero()
        count = adata.obs[labels_key][col].value_counts()
        count_list.append(count)

    neighborhood_counts = pd.DataFrame(count_list, index=adata.obs_names)
    neighborhood_counts.fillna(0, inplace = True)
    neighborhood_freq = neighborhood_counts.div(neighborhood_counts.sum(axis = 1), axis = 0)

    cells_nonn = list(neighborhood_counts.index[np.where(neighborhood_counts.sum(axis = 1) < min_cells)[0]])
    adata_niche = anndata.AnnData(neighborhood_freq)
    adata_niche.obs = adata.obs.loc[neighborhood_freq.index, :]
    
    return adata_niche, cells_nonn


def compute_connectivities_umap(
    knn_indices,
    knn_dists,
    n_obs,
    n_neighbors,
    set_op_mix_ratio=1.0,
    local_connectivity=1.0,
):
    """\
    This is from umap.fuzzy_simplicial_set [McInnes18]_.

    Given a set of data X, a neighborhood size, and a measure of distance
    compute the fuzzy simplicial set (here represented as a fuzzy graph in
    the form of a sparse matrix) associated to the data. This is done by
    locally approximating geodesic distance at each point, creating a fuzzy
    simplicial set for each such point, and then combining all the local
    fuzzy simplicial sets into a global one via a fuzzy union.
    """
    with warnings.catch_warnings():
        # umap 0.5.0
        warnings.filterwarnings("ignore", message=r"Tensorflow not installed")
        from umap.umap_ import fuzzy_simplicial_set

    X = coo_matrix(([], ([], [])), shape=(n_obs, 1))
    connectivities = fuzzy_simplicial_set(
        X,
        n_neighbors,
        None,
        None,
        knn_indices=knn_indices,
        knn_dists=knn_dists,
        set_op_mix_ratio=set_op_mix_ratio,
        local_connectivity=local_connectivity,
    )

    if isinstance(connectivities, tuple):
        # In umap-learn 0.4, this returns (result, sigmas, rhos)
        connectivities = connectivities[0]

    distances = _get_sparse_matrix_from_indices_distances_umap(
        knn_indices, knn_dists, n_obs, n_neighbors
    )

    return distances, connectivities.tocsr()

def _get_sparse_matrix_from_indices_distances_umap(
    knn_indices, knn_dists, n_obs, n_neighbors
):
    rows = np.zeros((n_obs * n_neighbors), dtype=np.int64)
    cols = np.zeros((n_obs * n_neighbors), dtype=np.int64)
    vals = np.zeros((n_obs * n_neighbors), dtype=np.float64)

    for i in range(knn_indices.shape[0]):
        for j in range(n_neighbors):
            if knn_indices[i, j] == -1:
                continue  # We didn't get the full knn for i
            if knn_indices[i, j] == i:
                val = 0.0
            else:
                val = knn_dists[i, j]

            rows[i * n_neighbors + j] = i
            cols[i * n_neighbors + j] = knn_indices[i, j]
            vals[i * n_neighbors + j] = val

    result = coo_matrix((vals, (rows, cols)), shape=(n_obs, n_obs))
    result.eliminate_zeros()
    return result.tocsr()

def compute_neighbors_umap(X = None,
    n_neighbors = 10,
    random_state = None,
    metric = 'euclidean',
    metric_kwds = None, 
    angular = False,
    verbose = False
):
    """This is from umap.fuzzy_simplicial_set [McInnes18]_.

    Given a set of data X, a neighborhood size, and a measure of distance
    compute the fuzzy simplicial set (here represented as a fuzzy graph in
    the form of a sparse matrix) associated to the data. This is done by
    locally approximating geodesic distance at each point, creating a fuzzy
    simplicial set for each such point, and then combining all the local
    fuzzy simplicial sets into a global one via a fuzzy union.

    Parameters
    ----------
    X: array of shape (n_samples, n_features)
        The data to be modelled as a fuzzy simplicial set.
    n_neighbors
        The number of neighbors to use to approximate geodesic distance.
        Larger numbers induce more global estimates of the manifold that can
        miss finer detail, while smaller values will focus on fine manifold
        structure to the detriment of the larger picture.
    random_state
        A state capable being used as a numpy random state.
    metric
        The metric to use to compute distances in high dimensional space.
        If a string is passed it must match a valid predefined metric. If
        a general metric is required a function that takes two 1d arrays and
        returns a float can be provided. For performance purposes it is
        required that this be a numba jit'd function. Valid string metrics
        include:
            * euclidean
            * manhattan
            * chebyshev
            * minkowski
            * canberra
            * braycurtis
            * mahalanobis
            * wminkowski
            * seuclidean
            * cosine
            * correlation
            * haversine
            * hamming
            * jaccard
            * dice
            * russelrao
            * kulsinski
            * rogerstanimoto
            * sokalmichener
            * sokalsneath
            * yule
        Metrics that take arguments (such as minkowski, mahalanobis etc.)
        can have arguments passed via the metric_kwds dictionary. At this
        time care must be taken and dictionary elements must be ordered
        appropriately; this will hopefully be fixed in the future.
    metric_kwds
        Arguments to pass on to the metric, such as the ``p`` value for
        Minkowski distance.
    angular
        Whether to use angular/cosine distance for the random projection
        forest for seeding NN-descent to determine approximate nearest
        neighbors.
    verbose
        Whether to report information on the current progress of the algorithm.

    Returns
    -------
    **knn_indices**, **knn_dists** : np.arrays of shape (n_observations, n_neighbors)
    """
    with warnings.catch_warnings():
        # umap 0.5.0
        warnings.filterwarnings("ignore", message=r"Tensorflow not installed")
        from umap.umap_ import nearest_neighbors

    knn_indices, knn_dists, forest = nearest_neighbors(
        X,
        n_neighbors,
        random_state=random_state,
        metric=metric,
        metric_kwds=metric_kwds,
        angular=angular,
        verbose=verbose,
    )

    return knn_indices, knn_dists, forest
