import os
import glob
import anndata
import re
import numpy as np
import pandas as pd
from sklearn.cluster import *
from sklearn.metrics import calinski_harabasz_score, silhouette_score, davies_bouldin_score, accuracy_score, r2_score, mean_squared_error, confusion_matrix
from warnings import warn
import seaborn as sns
from sklearn.preprocessing import StandardScaler, LabelEncoder
import matplotlib.pyplot as plt
from sklearn.linear_model import ElasticNet
from sklearn.model_selection import GridSearchCV, LeaveOneGroupOut, StratifiedKFold, GridSearchCV
from tqdm import tqdm_notebook, tqdm
from IPython import get_ipython
import logging
import itertools
from skimage import morphology
from scipy.ndimage import gaussian_filter
from skimage.measure import label
import skimage.io as io
from sklearn.ensemble import RandomForestClassifier
import scipy
from scipy.stats import spearmanr, ttest_ind, ttest_rel, wilcoxon, mannwhitneyu, pearsonr
from alpineer import io_utils
from ark.segmentation import marker_quantification
from alpineer import load_utils, data_utils
from alpineer.io_utils import list_folders
import pickle
from sklearn.cluster import KMeans
from sklearn.pipeline import make_pipeline
from sklearn import preprocessing
import skimage.io as io
from skimage.measure import block_reduce
from skimage.segmentation import find_boundaries
from matplotlib import cm
from matplotlib import colors

def create_combined_channel_mask(chans, channel_dir, percentiles, threshold, smooth_val,
                                 erode_val):
    """
    Creates a mask for the total foreground of all channels in chans
    1. normalize channels in the image by the 99.9 percentile
    2. clip values that were outliers to 1 (>99.9 percentile)
    3. create a combined mask by summing pixel intensities across all channels
    4. apply gaussian blur
    5. create a mask if intensities > threshold 
    6. apply binary erosion 
    """

    normalized_chans = []
    for chan in chans:
        current_img = io.imread(os.path.join(channel_dir, chan + '.tiff')) #normalize data by the 99.9 percentile to account for outliers
        current_img /= percentiles[chan]
        current_img[current_img > 1] = 1
        normalized_chans.append(current_img)

    normalized_chans = np.stack(normalized_chans, axis=0)
    normalized_chans = np.sum(normalized_chans, axis=0) #sum across all channels to generate a combined mask

    smoothed = gaussian_filter(normalized_chans, sigma=smooth_val) #apply Gaussian filter
    mask = smoothed > threshold
    for _ in range(erode_val):
        mask = morphology.binary_erosion(mask)

    return mask

def compute_ECM_percentiles(channel_dir, out_dir, fovs, img_sub_folder, max_image_size, channels, percentile):
    image_data = load_utils.load_imgs_from_tree(channel_dir,
                                                fovs=fovs,
                                                img_sub_folder=img_sub_folder,
                                                max_image_size=max_image_size,
                                                channels=channels)

    # calculate image percentiles, average 99.9th percentile across example fovs for each ECM channel of interest
    percentiles = {}
    for chan in image_data.channels.values:
        current_data = image_data.loc[:, :, :, chan]
        chan_percentiles = []
        for i in range(current_data.shape[0]):
            current_img = current_data.values[i, :, :]
            chan_percentiles.append(np.percentile(current_img[current_img > 0], percentile))
        percentiles[chan] = (np.mean(chan_percentiles))

    # save the percentiles
    pd.DataFrame(percentiles, index=['percentile']).to_csv(os.path.join(out_dir, 'percentiles.csv'))

    # load the percentiles
    percentiles = pd.read_csv(os.path.join(out_dir, 'percentiles.csv'), index_col=0).to_dict(orient='records')[0]
    return percentiles

def kmeans_ECM_tiles(channel_dir, mask_dir, out_dir, channels, crop_size, fov_subset, n_clusters, threshold, random_state):

    tiled_crops = generate_crop_sum_dfs(channel_dir=channel_dir,
                                          mask_dir=mask_dir,
                                          channels=channels,
                                          crop_size=crop_size, fovs=fov_subset, cell_table=None)

    tiled_crops = normalize_by_ecm_area(crop_sums=tiled_crops, crop_size=crop_size,
                                            channels=channels)

    # create a pipeline for normalization and clustering the data
    kmeans_pipe = make_pipeline(preprocessing.PowerTransformer(method='yeo-johnson', standardize=True), #apply transformation and perform clustering (k = 2) (looking for cold and hot collagen)
                                KMeans(n_clusters=n_clusters, random_state=random_state))

    # select subset of data to train on
    no_ecm_mask = tiled_crops.ecm_fraction < threshold
    train_data = tiled_crops[~no_ecm_mask] #train on data where there is a high percentage of ECM
    train_data = train_data.loc[:, channels]

    # fit the pipeline on the data
    kmeans_pipe.fit(train_data.values)

    # save the trained pipeline
    pickle.dump(kmeans_pipe, open(os.path.join(out_dir, 'tile_classification_kmeans_pipe.pkl'), 'wb'))

    # load the model
    kmeans_pipe = pickle.load(open(os.path.join(out_dir, 'tile_classification_kmeans_pipe.pkl'), 'rb'))
    kmeans_preds = kmeans_pipe.predict(tiled_crops[channels].values) #predict on all data

    # get the transformed intermediate data
    transformed_data = kmeans_pipe.named_steps['powertransformer'].transform(tiled_crops[channels].values)
    transformed_df = pd.DataFrame(transformed_data, columns=channels)
    transformed_df['tile_cluster'] = kmeans_preds
    tiled_crops['tile_cluster'] = kmeans_preds
    tiled_crops.loc[no_ecm_mask, 'tile_cluster'] = -1 #data with low fraction gets set to -1 (no ECM)
    transformed_df.loc[no_ecm_mask, 'tile_cluster'] = -1

    # generate average image for each cluster
    cluster_means = transformed_df[~no_ecm_mask].groupby('tile_cluster').mean()
    return tiled_crops, kmeans_pipe, cluster_means, no_ecm_mask

def plot_stiched_ECM(tiled_crops, cluster_channels, no_ecm_mask, percentiles, crop_size, n_examples, channel_dir, mask_dir, out_dir):
    # create a stitched image with example images from each cluster
    for cluster in tiled_crops.tile_cluster.unique():
        if cluster == 'No_ECM':
            continue
        cluster_data = tiled_crops[(~no_ecm_mask) & (tiled_crops.tile_cluster == cluster)]
        cluster_data = cluster_data.sample(n=n_examples, random_state=0)

        stitched_img = np.zeros((crop_size * n_examples, crop_size * (len(cluster_channels) + 1)))
        for i in range(n_examples):
            fov_name = cluster_data.iloc[i]['fov']
            row_start = cluster_data.iloc[i]['row_coord']
            col_start = cluster_data.iloc[i]['col_coord']

            for j, chan in enumerate(cluster_channels):
                img = io.imread(os.path.join(channel_dir, fov_name, chan + '.tiff'))
                img_subset = img[row_start:row_start + crop_size, col_start:col_start + crop_size]
                img_subset = img_subset / percentiles[chan]
                img_subset[img_subset > 1] = 1

                stitched_img[i * crop_size:(i + 1) * crop_size, j * crop_size:(j + 1) * crop_size] = img_subset

            # do the same thing for the ecm mask
            img = io.imread(os.path.join(mask_dir, fov_name, 'total_ecm.tiff'))
            img_subset = img[row_start:row_start + crop_size, col_start:col_start + crop_size]
            stitched_img[i * crop_size:(i + 1) * crop_size, -crop_size:] = img_subset

        io.imsave(os.path.join(out_dir, 'cluster_' + str(cluster) + '.tiff'), stitched_img.astype('float32'),
                    check_contrast=False)

def kmeans_cell_ECM_tiles(channel_dir, mask_dir, cell_table_clusters, channels, crop_size, fov_subset, kmeans_pipe, threshold):
    # generate crops around cells to classify using the trained model
    cell_table_clusters = cell_table_clusters[['fov', 'centroid-0', 'centroid-1', 'label']]

    cell_crops = generate_crop_sum_dfs(channel_dir=channel_dir,
                                            mask_dir=mask_dir,
                                            channels=channels,
                                            crop_size=crop_size, fovs=fov_subset,
                                            cell_table=cell_table_clusters)

    # normalize based on ecm area
    cell_crops = normalize_by_ecm_area(crop_sums=cell_crops, crop_size=crop_size,
                                            channels=channels)

    cell_classifications = kmeans_pipe.predict(cell_crops[channels].values.astype('float64'))
    cell_crops['ecm_cluster'] = cell_classifications

    no_ecm_mask_cell = cell_crops.ecm_fraction < threshold

    cell_crops.loc[no_ecm_mask_cell, 'ecm_cluster'] = -1
    return cell_crops

def relabel_segmentation(labeled_image, labels_dict):
    """Takes a labeled image and translates its labels according to a dictionary.

    Returns the relabeled array (according to the dictionary).

    Args:
        labeled_image (numpy.ndarray):
            2D numpy array of labeled cell objects.
        labels_dict (dict):
            a mapping between labeled cells and their clusters.
    Returns:
        numpy.ndarray:
            The relabeled array.
    """

    img = np.copy(labeled_image)
    unique_cell_ids = np.unique(labeled_image)
    unique_cell_ids = unique_cell_ids[np.nonzero(unique_cell_ids)]

    default_label = max(labels_dict.values()) + 1
    for cell_id in unique_cell_ids:
        img[labeled_image == cell_id] = labels_dict.get(cell_id, default_label)
    return img

def create_cell_overlay(cell_table, seg_folder, fovs, cluster_col):
    cell_subset = cell_table.copy()
    cell_subset['unique_ids'] = pd.factorize(cell_subset[cluster_col])[0] + 1

    categories = cell_subset[[cluster_col, 'unique_ids']].drop_duplicates()[cluster_col].values

    # import viridis colormap from mpl
    num_categories = np.max(cell_subset.unique_ids)
    cm_values = cm.get_cmap('Paired', num_categories)

    # get RGB values from cm_values
    rgb_values = cm_values(np.arange(num_categories))

    # combine with all black for background
    rgb_values = np.vstack((np.array([0, 0, 0, 1]), rgb_values))

    new_cmap = colors.ListedColormap(rgb_values)

    bounds = [i-0.5 for i in np.linspace(0, num_categories+1, num_categories+2)]
    norm = colors.BoundaryNorm(bounds, new_cmap.N + 1)

    for idx, image in enumerate(fovs):
        seg_mask = io.imread(os.path.join(seg_folder, image + '_whole_cell.tiff'))

        edges = find_boundaries(seg_mask, mode='inner')
        seg_mask = np.where(edges == 0, seg_mask, 0)

        # convert string entries in pandas df to unique integers
        cell_subset_plot = cell_subset[cell_subset['fov'] == image]
        labels_dict = dict(zip(cell_subset_plot['label'], cell_subset_plot['unique_ids']))

        # relabel the array
        relabeled_img_array = relabel_segmentation(seg_mask, labels_dict)

        im = plt.imshow(relabeled_img_array, cmap=new_cmap, norm=norm)
        tick_names = ['Empty'] + categories.tolist()
        cbar = plt.colorbar(im, ticks=np.arange(len(tick_names)))
        cbar.set_ticks(cbar.ax.get_yticks())
        cbar.ax.set_yticklabels(tick_names)
        plt.show()

def make_directory(directory: str = None):
    """Creates a directory at the specified path if one doesn't exist.

    Parameters
    ----------
    directory : str
        A string specifying the directory path.

    Returns
    -------
    """
    if not os.path.exists(directory):
        os.makedirs(directory)
        
def standardize(x =  None):
    """Standardizes data by removing the mean and scaling to unit variance.

    Parameters
    x: pd.DataFrame (default = None)
        data matrix (dimensions = cells x features)
    ----------

    Returns
    X: pd.DataFrame
        standardized data matrix (dimensions = cells x features)
    ----------
    """
    scaler = StandardScaler(with_mean = True, with_std = True)
    X = scaler.fit_transform(x)
    return X

def create_channel_mask(img, intensity_thresh, sigma, min_mask_size=0, max_hole_size=100000):
    """Generates a binary mask from a single channel image

    Args:
        img (np.ndarray): image to be masked
        intensity_thresh (float): threshold for the image intensity to use for masking
        sigma (float): sigma for gaussian blur
        min_mask_size (int): minimum size of masked objects to include
        max_hole_size (int): maximum size of holes to leave in masked objects
        """
    # create a binary mask
    img_smoothed = gaussian_filter(img.astype(float), sigma=sigma)
    img_mask = img_smoothed > intensity_thresh

    # if no post-processing return as is
    if min_mask_size == 0:
        return img_mask

    # otherwise, clean up the mask before returning
    label_mask = label(img_mask)
    label_mask = morphology.remove_small_objects(label_mask, min_size=min_mask_size)
    label_mask = morphology.remove_small_holes(label_mask, area_threshold=max_hole_size)

    return label_mask


def create_cell_mask(seg_mask, cell_table, fov_name, cluster_col, cell_types, sigma=10, smooth_thresh=0.3,
                     min_mask_size=0, max_hole_size=100000):
    """Generates a binary from the cells listed in `cell_types`

    args:
        seg_mask (numpy.ndarray): segmentation mask
        cell_table (pandas.DataFrame): cell table containing segmentation IDs and cell types
        fov_name (str): name of the fov to process
        cell_types (list): list of cell types to include in the mask
        sigma (float): sigma for gaussian smoothing
        smooth_thresh (float): threshold for including a pixel in the smoothed mask
        min_mask_size (int): minimum size of a mask to include
        max_hole_size (int): maximum size of a hole to leave without filling

    returns:
        numpy.ndarray: binary mask
    """
    # get cell labels for fov and cell type
    cell_subset = cell_table[cell_table['fov'] == fov_name]
    cell_subset = cell_subset[cell_subset[cluster_col].isin(cell_types)]
    cell_labels = cell_subset['label'].values

    # create mask for cell type
    cell_mask = np.isin(seg_mask, cell_labels)

    # postprocess mask
    cell_mask = create_channel_mask(img=cell_mask, intensity_thresh=smooth_thresh,
                                    sigma=sigma, min_mask_size=min_mask_size,
                                    max_hole_size=max_hole_size)

    return cell_mask


def create_cancer_boundary(img, seg_mask, min_mask_size=3500, max_hole_size=1000,
                           border_size=50, channel_thresh=0.0015):
    """Generate masks representing different tumor regions"""
    img_smoothed = gaussian_filter(img, sigma=10)
    img_mask = img_smoothed > channel_thresh

    # clean up mask prior to analysis
    img_mask = np.logical_or(img_mask, seg_mask)
    label_mask = label(img_mask)
    label_mask = morphology.remove_small_objects(label_mask, min_size=min_mask_size)
    label_mask = morphology.remove_small_holes(label_mask, area_threshold=max_hole_size)

    # define external borders
    external_boundary = morphology.binary_dilation(label_mask)

    for _ in range(border_size):
        external_boundary = morphology.binary_dilation(external_boundary)

    external_boundary = external_boundary.astype(int) - label_mask.astype(int)
    # create interior borders
    interior_boundary = morphology.binary_erosion(label_mask)

    for _ in range(border_size):
        interior_boundary = morphology.binary_erosion(interior_boundary)

    interior_boundary = label_mask.astype(int) - interior_boundary.astype(int)

    combined_mask = np.ones_like(img)
    combined_mask[label_mask] = 4
    combined_mask[external_boundary > 0] = 2
    combined_mask[interior_boundary > 0] = 3

    return combined_mask


def calculate_mask_areas(mask_dir, fovs):
    """Calculate the area of each mask per fov

    Args:
        mask_dir (str): path to directory containing masks for each fov
        fovs (list): list of fovs to calculate mask areas for

    Returns
        pd.DataFrame: dataframe containing the area of each mask per fov
    """
    # get list of masks
    mask_files = io_utils.list_files(os.path.join(mask_dir, fovs[0]))
    mask_names = [os.path.splitext(os.path.basename(x))[0] for x in mask_files]

    # loop through fovs and masks to calculate area
    area_dfs = []
    for fov in fovs:
        mask_areas = []
        for mask_file in mask_files:
            mask = io.imread(os.path.join(mask_dir, fov, mask_file))
            mask_areas.append(np.sum(mask))

        area_df = pd.DataFrame({'compartment': mask_names, 'area': mask_areas,
                                'fov': fov})

        # separately calculate size for non-background compartment
        bg_area = area_df[area_df['compartment'] == 'empty_slide']['area'].values[0]
        foreground_area = mask.shape[0] ** 2 - bg_area
        blank_df = pd.DataFrame({'compartment': ['all'], 'area': [foreground_area], 'fov': [fov]})
        area_df = pd.concat([area_df, blank_df], ignore_index=True)

        area_dfs.append(area_df)

    return pd.concat(area_dfs, axis=0)


def assign_cells_to_mask(seg_dir, mask_dir, fovs):
    """Assign cells an image to the mask they overlap most with

    Args:
        seg_dir (str): path to segmentation directory
        mask_dir (str): path to mask directory, with masks for each FOV in a dedicated folder
        fovs (list): list of fovs to process

    Returns:
        pandas.DataFrame: dataframe with cell assignments to masks
    """

    # extract counts of each mask per cell
    normalized_cell_table, _ = marker_quantification.generate_cell_table(segmentation_dir=seg_dir,
                                                                         tiff_dir=mask_dir,
                                                                         fovs=fovs,
                                                                         img_sub_folder='')
    # drop cell_size column
    normalized_cell_table = normalized_cell_table.drop(columns=['cell_size'])

    # move fov column to front
    fov_col = normalized_cell_table.pop('fov')
    normalized_cell_table.insert(0, 'fov', fov_col)

    # remove all columns after label
    normalized_cell_table = normalized_cell_table.loc[:, :'label']

    # move label column to front
    label_col = normalized_cell_table.pop('label')
    normalized_cell_table.insert(1, 'label', label_col)

    # create new column with name of column max for each row
    normalized_cell_table['mask_name'] = normalized_cell_table.iloc[:, 2:].idxmax(axis=1)

    return normalized_cell_table[['fov', 'label', 'mask_name']]


def identify_cell_bounding_box(row_centroid, col_centroid, crop_size, img_shape):
    """Finds the upper-left hand corner of a bounding box surrounding the cell, corrected for edges.

    Args:
        row_coord (int): row coordinate of the cell centroid
        col_coord (int): column coordinate of the cell centroid
        crop_size (int): size of the bounding box
        img_shape (tuple): shape of the image
    """

    # get the image dimensions
    img_height, img_width = img_shape

    # adjust the centroid to be at least crop_size / 2 away from the bottom right corner
    if row_centroid > img_height - crop_size // 2:
        row_centroid = img_height - crop_size // 2
    if col_centroid > img_width - crop_size // 2:
        col_centroid = img_width - crop_size // 2

    # set new coordinates to be crop_size / 2 up and to the left of the centroid
    col_coord = col_centroid - crop_size // 2
    row_coord = row_centroid - crop_size // 2

    # make sure the coordinates are not negative
    col_coord = max(col_coord, 0)
    row_coord = max(row_coord, 0)

    return int(row_coord), int(col_coord)


def generate_cell_crop_coords(cell_table_fov, crop_size, img_shape):
    """Generates the coordinates for cropping each cell in a fov

    Args:
        cell_table_fov (pd.DataFrame): dataframe containing the location of each cell
        crop_size (int): size of the bounding box
        img_shape (tuple): shape of the image

    Returns:
        pd.DataFrame: dataframe containing the coordinates for cropping each cell
    """
    # get the coordinates for each cell
    cell_coords = cell_table_fov[['centroid-0', 'centroid-1']].values

    # calculate the coordinates for the upper left hand corner of the bounding box
    crop_coords = [identify_cell_bounding_box(row_coord, col_coord, crop_size, img_shape)
                   for row_coord, col_coord in cell_coords]

    # create a dataframe with the coordinates
    crop_coords_df = pd.DataFrame(crop_coords, columns=['row_coord', 'col_coord'])

    # add the label column
    crop_coords_df['id'] = cell_table_fov['label'].values

    return crop_coords_df


def generate_tiled_crop_coords(crop_size, img_shape):
    """Generate coordinates for uniformly tiled crops

    Args:
        crop_size (int): size of the bounding box
        img_shape (tuple): shape of the image

    Returns:
        pd.DataFrame: dataframe containing the coordinates for cropping each tile
    """

    # compute all combinations of start coordinates
    img_height, img_width = img_shape

    row_coords = np.arange(0, img_height, crop_size)
    col_coords = np.arange(0, img_width, crop_size)
    coords = itertools.product(row_coords, col_coords)

    # create a dataframe with the coordinates
    crop_coords_df = pd.DataFrame(coords, columns=['row_coord', 'col_coord'])

    # add a column for the tile combination
    crop_coords_df['id'] = [f'row_{row}_col_{col}' for row, col in zip(crop_coords_df['row_coord'],
                                                                       crop_coords_df['col_coord'])]

    return crop_coords_df


def extract_crop_sums(img_data, crop_size, crop_coords_df):
    """Extracts and sums crops from an image

    Args:
        img_data (np.ndarray): image data for a single fov
        crop_size (int): size of the bounding box around each cell
        crop_coords_df (pd.DataFrame): dataframe containing the coordinates for cropping each tile

    Returns:
        np.ndarray: array of crop sums
    """
    # list to hold crop sums
    crop_sums = []

    for row_coord, col_coord in zip(crop_coords_df['row_coord'], crop_coords_df['col_coord']):
        # crop based on provided coords
        crop = img_data[row_coord:row_coord + crop_size,
                        col_coord:col_coord + crop_size, :]

        # sum the channels within the crop
        crop_sum = crop.sum(axis=(0, 1))

        # add the crop sum to the list
        crop_sums.append(crop_sum)

    return np.array(crop_sums)


def generate_crop_sum_dfs(channel_dir, mask_dir, channels, crop_size, fovs, cell_table):
    """Generates dataframes of summed crops around cells or tiles for each fov

    Args:
        channel_dir (str): path to the directory containing image data
        mask_dir (str): path to the directory containing the ecm masks
        channels (list): list of channels to extract crops from
        crop_size (int): size of the bounding box around each cell or tile
        fovs (list): list of fovs to process
        cell_table (pd.DataFrame): cell table, if None will tile the image


    Returns:
        pd.DataFrame: dataframe of summed crops around cells
    """
    # list to hold dataframes
    crop_df_list = []

    for fov in fovs:
        # load the image data
        img_data = load_utils.load_imgs_from_tree(channel_dir,
                                                  fovs=[fov],
                                                  img_sub_folder='',
                                                  channels=channels) #shape is 1 x 2048 x 2048 x channels
        ecm_mask = io.imread(os.path.join(mask_dir, fov, 'total_ecm.tiff'))

        # combine the image data and the ecm mask into numpy array 
        img_data = np.concatenate((img_data[0].values, ecm_mask[..., None]), axis=-1) #add ecm mask as channel

        # set logic based on whether or not a cell table is provided
        if cell_table is not None:
            cell_table_fov = cell_table[cell_table.fov == fov]
            cell_table_fov = cell_table_fov.reset_index(drop=True)

            # generate the coordinates for cropping each cell
            crop_coords_df = generate_cell_crop_coords(cell_table_fov=cell_table_fov,
                                                       crop_size=crop_size,
                                                       img_shape=img_data.shape[:-1])
        else:
            # generate the coordinates for cropping each tile
            crop_coords_df = generate_tiled_crop_coords(crop_size=crop_size,
                                                        img_shape=img_data.shape[:-1])

        # extract summed counts around each cell
        crop_sums = extract_crop_sums(img_data=img_data, crop_size=crop_size,
                                      crop_coords_df=crop_coords_df)

        # create a dataframe of the summed counts
        crop_sums_df = pd.DataFrame(crop_sums,
                                    columns=channels + ['ecm_mask'])

        # combine the crop_sums_df with the crop_coords_df
        crop_sums_df = pd.concat([crop_coords_df, crop_sums_df], axis=1)
        crop_sums_df['fov'] = fov

        # add the dataframe to the list
        crop_df_list.append(crop_sums_df)

    return pd.concat(crop_df_list, ignore_index=True)


# normalize by ecm area
def normalize_by_ecm_area(crop_sums, crop_size, channels):
    """Normalize the summed pixel values by the area of the ecm mask

    Args:
        crop_sums (pd.DataFrame): dataframe of crop sums
        crop_size (int): size of the crop
        channels (list): list of channels to normalize

    Returns:
        pd.DataFrame: normalized dataframe
    """

    crop_sums['ecm_fraction'] = (crop_sums['ecm_mask'] + 1) / (crop_size ** 2)
    crop_sums.loc[:, channels] = crop_sums.loc[:, channels].div(crop_sums['ecm_fraction'], axis=0)

    return crop_sums