from collections.abc import Iterable
import itertools
from typing import Literal
from matplotlib import ticker
from matplotlib.axes import Axes
from matplotlib.colors import Normalize
from matplotlib.figure import Figure
import numpy as np
import xarray as xr
import os
import pandas as pd
from skimage import morphology
from scipy.ndimage import gaussian_filter
from skimage.segmentation import find_boundaries
from skimage.measure import label
import skimage.io as io

from alpineer.misc_utils import verify_in_list
from alpineer import io_utils, load_utils
from ark.segmentation import marker_quantification

from scipy.stats import spearmanr, ttest_ind, ttest_rel, wilcoxon, mannwhitneyu, pearsonr

import matplotlib.pyplot as plt
import seaborn as sns


def find_conserved_features(paired_df, sample_name_1, sample_name_2, min_samples=20):
    """Given a df with paired features, find features that are conserved between two samples.

    Args:
        paired_df (pd.DataFrame): df with paired features
        sample_name_1 (str): name of sample 1
        sample_name_2 (str): name of sample 2
        min_samples (int): minimum number of samples to calculate correlation
    """

    # set up placeholders
    p_vals = []
    cors = []
    names = []

    # loop over each feature in the df
    for feature_name in paired_df.feature_name_unique.unique():

        # if either sample is missing a feature, skip that comparison
        values = paired_df[(paired_df.feature_name_unique == feature_name)].copy()
        values.dropna(inplace=True)

        # remove rows where both values are 0
        zero_mask = (values[sample_name_1] == 0) & (values[sample_name_2] == 0)
        values = values[~zero_mask]

        if len(values) > min_samples:
            cor, p_val = spearmanr(values[sample_name_1], values[sample_name_2])
            p_vals.append(p_val)
            cors.append(cor)
            names.append(feature_name)
        else:
            p_vals.append(np.nan)
            cors.append(np.nan)
            names.append(feature_name)

    ranked_features = pd.DataFrame({'feature_name_unique': names, 'p_val': p_vals, 'cor': cors})
    ranked_features['log_pval'] = -np.log10(ranked_features.p_val)

    # get ranking of each row by log_pval
    ranked_features['pval_rank'] = ranked_features.log_pval.rank(ascending=False)
    ranked_features['cor_rank'] = ranked_features.cor.rank(ascending=False)
    ranked_features['combined_rank'] = (ranked_features.pval_rank.values + ranked_features.cor_rank.values) / 2

    # generate consistency score
    max_rank = len(~ranked_features.cor.isna())
    normalized_rank = ranked_features.combined_rank / max_rank
    ranked_features['consistency_score'] = 1 - normalized_rank

    return ranked_features


def cluster_df_helper(cell_table, cluster_col_name, result_name, normalize=False):
    """Helper function which creates a df when no subsetting is required

        Args:
            cell_table: the dataframe containing information on each cell
            cluster_col_name: the column name in cell_table that contains the cluster information
            result_name: the name of this statistic in summarized information df
            normalize: whether to report the total or normalized counts in the result

        Returns:
            pd.DataFrame: long format dataframe containing the summarized data"""

    # group each fov by the supplied cluster column, then count and normalize
    grouped = cell_table.groupby(['fov'])
    counts = grouped[cluster_col_name].value_counts(normalize=normalize)
    counts = counts.unstack(level=cluster_col_name, fill_value=0).stack()

    # standardize the column names
    counts = counts.reset_index()
    counts['metric'] = result_name
    counts = counts.rename(columns={cluster_col_name: 'cell_type', 0: 'value'})

    return counts


def create_long_df_by_cluster(cell_table, cluster_col_name, result_name, subset_col=None,
                              normalize=False):
    """Summarize cell counts by cluster, with the option to subset by an additional feature

    Args:
        cell_table (pd.DataFrame): the dataframe containing information on each cell
        cluster_col_name (str): the column name in cell_table that contains the cluster information
        result_name (str): the name of this statistic in the returned df
        subset_col (str): the column name in cell_table to subset by
        normalize (bool): whether to report the total or normalized counts in the result

    Returns:
        pd.DataFrame: long format dataframe containing the summarized data"""

    # first generate df without subsetting
    long_df_all = cluster_df_helper(cell_table, cluster_col_name, result_name, normalize)
    long_df_all['subset'] = 'all'

    # if a subset column is specified, create df stratified by subset
    if subset_col is not None:
        verify_in_list(subset_col=subset_col, cell_table_columns=cell_table.columns)

        # group each fov by fov and cluster
        grouped = cell_table.groupby(['fov', subset_col])
        counts = grouped[cluster_col_name].value_counts(normalize=normalize)

        # unstack and restack to make sure that missing cell populations are filled with zeros
        counts = counts.unstack(level=cluster_col_name, fill_value=0).stack()

        # standardize the column names
        counts = counts.reset_index()
        counts['metric'] = result_name
        counts = counts.rename(columns={cluster_col_name: 'cell_type', subset_col: 'subset',
                                        0: 'value'})

        # combine the two dataframes
        long_df_all = pd.concat([long_df_all, counts], axis=0, ignore_index=True)

    return long_df_all


def functional_df_helper(func_table, cluster_col_name, drop_cols, result_name, normalize=False):
    """Function to summarize functional marker data by cell type

    Args:
        func_table (pd.DataFrame): cell table containing functional markers
        cluster_col_name (str): name of the column in func_table that contains the cluster information
        drop_cols (list): list of columns to drop from func_table
        result_name (str): name of the statistic in the summarized information df
        normalize (bool): whether to report the total or normalized counts in the result

    Returns:
        pd.DataFrame: long format dataframe containing the summarized data"""

    verify_in_list(cell_type_col=cluster_col_name, cell_table_columns=func_table.columns)
    verify_in_list(drop_cols=drop_cols, cell_table_columns=func_table.columns)

    # drop columns from table
    func_table_small = func_table.loc[:, ~func_table.columns.isin(drop_cols)]

    # group by specified columns
    grouped_table = func_table_small.groupby(['fov', cluster_col_name])
    if normalize:
        transformed = grouped_table.agg(np.mean)
    else:
        transformed = grouped_table.agg(np.sum)
    transformed.reset_index(inplace=True)

    # reshape to long df
    long_df = pd.melt(transformed, id_vars=['fov', cluster_col_name], var_name='functional_marker')
    long_df['metric'] = result_name
    long_df = long_df.rename(columns={cluster_col_name: 'cell_type'})

    return long_df


def create_long_df_by_functional(func_table, cluster_col_name, drop_cols, result_name,
                                        subset_col=None, normalize=False):
    """Summarize functional marker positivity by cell type, with the option to subset by an additional feature

    Args:
        func_table (pd.DataFrame): the dataframe containing information on each cell
        cluster_col_name (str): the column name in cell_table that contains the cluster information
        drop_cols (list): list of columns to drop from cell_table
        result_name (str): the name of this statistic in the returned df
        subset_col (str): the column name in cell_table to subset by
        normalize (bool): whether to report the total or normalized counts in the result

    Returns:
        pd.DataFrame: long format dataframe containing the summarized data"""

    # first generate df without subsetting
    drop_cols_all = drop_cols.copy()
    if subset_col is not None:
        drop_cols_all = drop_cols + [subset_col]

    long_df_all = functional_df_helper(func_table, cluster_col_name, drop_cols_all, result_name,
                                       normalize)
    long_df_all['subset'] = 'all'

    # if a subset column is specified, create df stratified by subset
    if subset_col is not None:
        verify_in_list(subset_col=subset_col, cell_table_columns=func_table.columns)

        # drop columns from table
        func_table_small = func_table.loc[:, ~func_table.columns.isin(drop_cols)]

        # group by specified columns
        grouped_table = func_table_small.groupby(['fov', subset_col, cluster_col_name])
        if normalize:
            transformed = grouped_table.agg(np.mean)
        else:
            transformed = grouped_table.agg(np.sum)
        transformed.reset_index(inplace=True)

        # reshape to long df
        long_df = pd.melt(transformed, id_vars=['fov', subset_col, cluster_col_name],
                            var_name='functional_marker')
        long_df['metric'] = result_name
        long_df = long_df.rename(columns={cluster_col_name: 'cell_type', subset_col: 'subset'})

        # combine the two dataframes
        long_df_all = pd.concat([long_df_all, long_df], axis=0, ignore_index=True)

    return long_df_all


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


def create_cell_mask(seg_mask, cell_table, fov_name, cell_types, sigma=10, smooth_thresh=0.3,
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
    cell_subset = cell_subset[cell_subset['cell_cluster_broad'].isin(cell_types)]
    cell_labels = cell_subset['label'].values

    # create mask for cell type
    cell_mask = np.isin(seg_mask, cell_labels)

    # postprocess mask
    cell_mask = create_channel_mask(img=cell_mask, intensity_thresh=smooth_thresh,
                                    sigma=sigma, min_mask_size=min_mask_size,
                                    max_hole_size=max_hole_size)

    return cell_mask


def create_cancer_boundary(img, seg_mask, sigma=10, min_mask_size=3500, max_hole_size=1000,
                           border_size=50, channel_thresh=0.0015):
    """Generate masks representing different tumor regions"""
    img_smoothed = gaussian_filter(img, sigma=sigma)
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
                                                  channels=channels)
        ecm_mask = io.imread(os.path.join(mask_dir, fov, 'total_ecm.tiff'))

        # combine the image data and the ecm mask into numpy array
        img_data = np.concatenate((img_data[0].values, ecm_mask[..., None]), axis=-1)

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


def compare_timepoints(feature_df, timepoint_1_name, timepoint_1_list, timepoint_2_name,
                       timepoint_2_list, paired=None, feature_suff='mean'):
    """Compute change in a feature across two timepoints.

    Args:
        feature_df (pd.DataFrame): dataframe containing features
        timepoint_1_name (str): overall name to give for the first timepoint
        timepoint_1_list (list): list of specific timepoints to include in the first group
        timepoint_2_name (str): overall name to give for the second timepoint
        timepoint_2_list (list): list of specific timepoints to include in the second group
        paired (str): column name to use for paired samples
        feature_suff (str): suffix to add to feature name
    """
    # get unique features
    features = feature_df.feature_name_unique.unique()

    feature_names = []
    timepoint_1_means = []
    timepoint_1_norm_means = []
    timepoint_2_means = []
    timepoint_2_norm_means = []
    log_pvals = []

    analysis_df = feature_df.loc[(feature_df.Timepoint.isin(timepoint_1_list + timepoint_2_list)), :]

    # subset to only include paired samples
    if paired is not None:
        analysis_df = analysis_df.loc[analysis_df[paired], :]
        if len(timepoint_1_list) != 1 or len(timepoint_2_list) != 1:
            raise ValueError('Paired samples only works with one timepoint per group.')

    # loop through each feature separately
    for feature_name_unique in features:
        values = analysis_df.loc[(analysis_df.feature_name_unique == feature_name_unique), :]

        # only keep samples with both timepoints
        if paired is not None:
            values_norm = values.pivot(index='Patient_ID', columns='Timepoint',
                                       values='normalized_mean')
            values_raw = values.pivot(index='Patient_ID', columns='Timepoint', values='raw_mean')
            values_norm = values_norm.dropna()
            values_raw = values_raw.dropna()

            # if there are no paired samples, set to nan
            if values_raw.shape[1] != 2 or len(values_raw) == 0:
                tp_1_vals, tp_1_norm_vals = np.array(np.nan), np.array(np.nan)
                tp_2_vals, tp_2_norm_vals = np.array(np.nan), np.array(np.nan)
            # get the columns corresponding to each timepoint
            else:
                tp_1_vals = values_raw[timepoint_1_list[0]].values
                tp_1_norm_vals = values_norm[timepoint_1_list[0]].values
                tp_2_vals = values_raw[timepoint_2_list[0]].values
                tp_2_norm_vals = values_norm[timepoint_2_list[0]].values

        # for unpaired, just subset to the timepoints of interest
        else:
            tp_1_vals = values.loc[
                values.Timepoint.isin(timepoint_1_list), 'raw_' + feature_suff].values
            tp_1_norm_vals = values.loc[
                values.Timepoint.isin(timepoint_1_list), 'normalized_' + feature_suff].values
            tp_2_vals = values.loc[
                values.Timepoint.isin(timepoint_2_list), 'raw_' + feature_suff].values
            tp_2_norm_vals = values.loc[
                values.Timepoint.isin(timepoint_2_list), 'normalized_' + feature_suff].values

            # if either timepoint is missing, set to nan
            if len(tp_1_vals) == 0 or len(tp_2_vals) == 0:
                tp_1_vals, tp_1_norm_vals = np.array(np.nan), np.array(np.nan)
                tp_2_vals, tp_2_norm_vals = np.array(np.nan), np.array(np.nan)

        timepoint_1_means.append(tp_1_vals.mean())
        timepoint_1_norm_means.append(tp_1_norm_vals.mean())
        timepoint_2_means.append(tp_2_vals.mean())
        timepoint_2_norm_means.append(tp_2_norm_vals.mean())

        # compute t-test for difference between timepoints
        if paired is not None:
            if np.all(tp_1_norm_vals - tp_2_norm_vals == 0):
                t, p = np.nan, np.nan
            else:
                t, p = wilcoxon(tp_1_norm_vals, tp_2_norm_vals)
        else:
            t, p = mannwhitneyu(tp_1_norm_vals, tp_2_norm_vals)

        log_pvals.append(-np.log10(p))

    # construct final df
    means_df = pd.DataFrame({timepoint_1_name + '_mean': timepoint_1_means,
                             timepoint_2_name + '_mean': timepoint_2_means,
                             timepoint_1_name + '_norm_mean': timepoint_1_norm_means,
                             timepoint_2_name + '_norm_mean': timepoint_2_norm_means,
                             'log_pval': log_pvals}, index=features)
    # calculate difference between timepoint 2 and timepoint 1
    means_df['mean_diff'] = means_df[timepoint_2_name + '_norm_mean'].values - means_df[timepoint_1_name + '_norm_mean'].values
    means_df = means_df.reset_index().rename(columns={'index': 'feature_name_unique'})

    return means_df


def compare_populations(feature_df, pop_col, pop_1, pop_2, timepoints, feature_suff='mean',
                        method='ttest'):
    """Compare difference in a feature between two populations.

    Args:
        feature_df (pd.DataFrame): dataframe containing features
        pop_col (str): column name containing population information
        pop_1 (str): name of the first population
        pop_2 (str): name of the second population
        timepoints (list): list of timepoints to include
        feature_suff (str): suffix to add to feature name
        method (str): method to use for comparing populations
    """
    # get unique features
    features = feature_df.feature_name_unique.unique()

    feature_names = []
    pop_1_means = []
    pop_1_norm_means = []
    pop_1_norm_meds = []
    pop_2_means = []
    pop_2_norm_means = []
    pop_2_norm_meds = []
    log_pvals = []

    analysis_df = feature_df.loc[(feature_df.Timepoint.isin(timepoints)), :]

    for feature_name in features:
        values = analysis_df.loc[(analysis_df.feature_name_unique == feature_name), :]
        pop_1_vals = values.loc[values[pop_col] == pop_1, 'raw_' + feature_suff].values
        pop_1_norm_vals = values.loc[values[pop_col] == pop_1, 'normalized_' + feature_suff].values
        pop_2_vals = values.loc[values[pop_col] == pop_2, 'raw_' + feature_suff].values
        pop_2_norm_vals = values.loc[values[pop_col] == pop_2, 'normalized_' + feature_suff].values

        # if insufficient number of samples, ignore that comparison
        if len(pop_1_vals) < 3 or len(pop_2_vals) < 3:
            pop_1_vals, pop_1_norm_vals = np.array(np.nan), np.array(np.nan)
            pop_2_vals, pop_2_norm_vals = np.array(np.nan), np.array(np.nan)

        pop_1_means.append(pop_1_vals.mean())
        pop_1_norm_means.append(pop_1_norm_vals.mean())
        pop_1_norm_meds.append(np.median(pop_1_norm_vals))
        pop_2_means.append(pop_2_vals.mean())
        pop_2_norm_means.append(pop_2_norm_vals.mean())
        pop_2_norm_meds.append(np.median(pop_2_norm_vals))

        # compute difference between timepoints
        if method == 'ttest':
            t, p = ttest_ind(pop_1_norm_vals, pop_2_norm_vals)
        else:
            t, p = mannwhitneyu(pop_1_norm_vals, pop_2_norm_vals)

        log_pvals.append(-np.log10(p))

    means_df = pd.DataFrame({pop_1 + '_mean': pop_1_means,
                             pop_2 + '_mean': pop_2_means,
                             pop_1 + '_norm_mean': pop_1_norm_means,
                             pop_2 + '_norm_mean': pop_2_norm_means,
                             pop_1 + '_norm_med': pop_1_norm_meds,
                             pop_2 + '_norm_med': pop_2_norm_meds,
                             'log_pval': log_pvals}, index=features)
    # calculate difference
    means_df['mean_diff'] = means_df[pop_2 + '_norm_mean'].values - means_df[pop_1 + '_norm_mean'].values
    means_df['med_diff'] = means_df[pop_2 + '_norm_med'].values - means_df[pop_1 + '_norm_med'].values
    means_df = means_df.reset_index().rename(columns={'index': 'feature_name_unique'})

    return means_df


def compare_continuous(feature_df, variable_col, min_samples=20, feature_suff='mean', method='spearman'):
    # set up placeholders
    p_vals = []
    cors = []
    names = []

    # loop over each feature in the df
    for feature_name in feature_df.feature_name_unique.unique():

        # if either sample is missing a feature, skip that comparison
        values = feature_df[(feature_df.feature_name_unique == feature_name)].copy()
        values.dropna(inplace=True)

        if len(values) > min_samples:
            if method == 'spearman':
                cor, p_val = spearmanr(values[variable_col], values['normalized_' + feature_suff])
            elif method == 'pearson':
                cor, p_val = pearsonr(values[variable_col], values['normalized_' + feature_suff])
            p_vals.append(p_val)
            cors.append(cor)
            names.append(feature_name)
        else:
            p_vals.append(np.nan)
            cors.append(np.nan)
            names.append(feature_name)

    ranked_features = pd.DataFrame({'feature_name_unique': names, 'p_val': p_vals, 'cor': cors})
    ranked_features['log_pval'] = -np.log10(ranked_features.p_val)

    # get ranking of each row by log_pval
    ranked_features['pval_rank'] = ranked_features.log_pval.rank(ascending=False)
    ranked_features['cor_rank'] = ranked_features.cor.abs().rank(ascending=False)
    ranked_features['combined_rank'] = (ranked_features.pval_rank.values + ranked_features.cor_rank.values) / 2

    # generate consistency score
    max_rank = len(~ranked_features.cor.isna())
    normalized_rank = ranked_features.combined_rank / max_rank
    ranked_features['feature_score'] = 1 - normalized_rank

    ranked_features = ranked_features.sort_values('feature_score', ascending=False)

    return ranked_features




def summarize_timepoint_enrichment(input_df, feature_df, timepoints, output_dir, pval_thresh=2,
                                   diff_thresh=0.3, plot_type='strip', sort_by='mean_diff'):
    """Generate a summary of the timepoint enrichment results

    Args:
        input_df (pd.DataFrame): dataframe containing timepoint enrichment results
        feature_df (pd.DataFrame): dataframe containing feature information
        timepoints (list): list of timepoints to include
        output_dir (str): path to output directory
        pval_thresh (float): threshold for p-value
        diff_thresh (float): threshold for difference between timepoints
    """

    input_df_filtered = input_df.loc[(input_df.log_pval > pval_thresh) & (np.abs(input_df[sort_by]) > diff_thresh), :]

    input_df_filtered = input_df_filtered.sort_values(sort_by, ascending=False)

    # plot the results
    for idx, feature in enumerate(input_df_filtered.feature_name_unique):
        feature_subset = feature_df.loc[(feature_df.feature_name_unique == feature), :]
        feature_subset = feature_subset.loc[(feature_subset.Timepoint.isin(timepoints)), :]

        g = sns.catplot(data=feature_subset, x='Timepoint', y='raw_mean', kind=plot_type, color='grey')
        g.fig.suptitle(feature)
        g.savefig(os.path.join(output_dir, 'Evolution_{}_{}.png'.format(idx, feature)))
        plt.close()

    sns.catplot(data=input_df_filtered, x=sort_by, y='feature_name_unique', kind='bar', color='grey')
    plt.savefig(os.path.join(output_dir, 'Timepoint_summary.png'))
    plt.close()


def summarize_population_enrichment(input_df, feature_df, timepoints, pop_col, output_dir, pval_thresh=2, diff_thresh=0.3,
                                    sort_by='mean_diff', plot_type='strip'):
    """Generate a summary of the population enrichment results

    Args:
        input_df (pd.DataFrame): dataframe containing population enrichment results
        feature_df (pd.DataFrame): dataframe containing feature information
        timepoints (list): list of timepoints to include
        pop_col (str): column name containing population information
        output_dir (str): path to output directory
        pval_thresh (float): threshold for p-value
        diff_thresh (float): threshold for difference between timepoints
    """

    input_df_filtered = input_df.loc[(input_df.log_pval > pval_thresh) & (np.abs(input_df[sort_by]) > diff_thresh), :]

    input_df_filtered = input_df_filtered.sort_values(sort_by, ascending=False)

    # plot the results
    for idx, feature in enumerate(input_df_filtered.feature_name_unique):
        feature_subset = feature_df.loc[(feature_df.feature_name_unique == feature), :]
        feature_subset = feature_subset.loc[(feature_subset.Timepoint.isin(timepoints)), :]
        if len(feature_subset[pop_col].unique()) != 2:
            continue
        g = sns.catplot(data=feature_subset, x=pop_col, y='raw_mean', kind=plot_type)
        g.fig.suptitle(feature)
        g.savefig(os.path.join(output_dir, 'Evolution_{}_{}.png'.format(idx, feature)))
        plt.close()

    if len(input_df_filtered) == 0:
        return
    sns.catplot(data=input_df_filtered, x=sort_by, y='feature_name_unique', kind='bar', color='grey')
    plt.savefig(os.path.join(output_dir, 'Evolution_summary.png'))
    plt.close()


def summarize_continuous_enrichment(input_df, feature_df, variable_col, timepoint, output_dir, min_score=0.85):

    # generate plot for best ranked features
    plot_features = input_df.loc[input_df.feature_score >= min_score, :]

    # sort by combined rank
    plot_features.sort_values(by='feature_score', inplace=True, ascending=False)

    for i in range(len(plot_features)):
        feature_name = plot_features.iloc[i].feature_name_unique
        feature_score = plot_features.iloc[i].feature_score
        values = feature_df[(feature_df.feature_name_unique == feature_name)].copy()
        values = values.loc[values.Timepoint == timepoint, :]
        values.dropna(inplace=True)

        # add leading 0 for plot
        str_pos = str(i)
        if i < 10:
            str_pos = '0' + str_pos

        # plot
        sns.scatterplot(data=values, x='raw_mean', y=variable_col)
        correlation, p_val = spearmanr(values.raw_mean, values[variable_col])
        plt.title(f'{feature_name} (score={feature_score:.2f}, r={correlation:.2f}, p={p_val:.2f})')
        plt.savefig(os.path.join(output_dir, f'{str_pos}_{feature_name}.png'))
        plt.close()


def compute_feature_enrichment(feature_df, inclusion_col, analysis_col):
    """Compute the enrichment of a set of included features from all features

    Args:
        feature_df (pd.DataFrame): DataFrame with features and inclusion column
        inclusion_col (str): Column with boolean values indicating inclusion
        analysis_col (str): Column to group by for analysis

    Returns:
        feature_props (pd.DataFrame): DataFrame with feature proportions
    """
    # aggregate based on analysis column
    feature_props = feature_df[[inclusion_col, analysis_col]].groupby([inclusion_col, analysis_col]).size().reset_index()
    feature_props = feature_props.pivot(index=analysis_col, columns=inclusion_col, values=0)

    # identify the fraction of features that are included
    selected_frac = np.sum(feature_df[inclusion_col]) / len(feature_df)

    # compute absolute and relative proportions
    feature_props['prop'] = feature_props[True] / (feature_props[True] + feature_props[False])
    feature_props['log2_ratio'] = np.log2(feature_props['prop'] / selected_frac)

    feature_props.reset_index(inplace=True)
    feature_props.sort_values(by='log2_ratio', inplace=True, ascending=False)

    return feature_props


class QuantileNormalization(Normalize):
    def __init__(self,
                 vmin: float = None,
                 vmax: float = None,
                 q: tuple[float, float] = (0.01, 0.99),
                 clip: bool = False,
                 eps: float = 1e-20,
    ) -> None:
        """Normalizes the input data by the qth quantile.

        Args
        ----------
        vmin : float, optional
            If vmin is not given it is initilaized from the minimum value of the
            array, by default None
        vmax : float, optional
            If vmax is not given it is initilaized from the maximum value of the
            array, by default None
        q : tuple[float, float], optional
            A tuple of quatiles where the smallest element is the minimum quantile
            and the largest element is the maximum percentile, by default (0.01, 0.99). Must
            be between 0 and 1 inclusive.
        clip: bool, optional
            If True, the normalized values are clipped to the range [0, 1], by default False
        eps: float, optional
            Small value to add to the denominator to avoid division by zero, by default 1e-20
        """
        super().__init__(vmin, vmax)
        if isinstance(q, tuple):
            if len(q) != 2:
                raise ValueError("q must be a tuple of length 2")
            if not all(0 <= i <= 1 for i in q):
                raise ValueError("q's elements must be between 0 and 1 inclusive")
        else:
            raise ValueError("q must be a tuple")
        
        self.qmin = min(q)
        self.qmax = max(q)
        self.clip = clip
        self.eps = eps

    def __call__(self, value):
        val_qmin, val_qmax = np.quantile(value, [self.qmin, self.qmax])
        
        norm = (value - val_qmin) / (val_qmax - val_qmin + self.eps)
        if self.clip:
            norm = np.clip(norm, 0, 1)
            
        return norm


def _remove_x_axis_ticks(ax: plt.Axes) -> None:
    ax.xaxis.set_major_locator(ticker.NullLocator())
    ax.xaxis.set_major_formatter(ticker.NullFormatter())
    ax.xaxis.set_minor_locator(ticker.NullLocator())
    ax.xaxis.set_minor_formatter(ticker.NullFormatter())


def _remove_y_axis_ticks(ax: plt.Axes) -> None:
    ax.yaxis.set_major_locator(ticker.NullLocator())
    ax.yaxis.set_major_formatter(ticker.NullFormatter())
    ax.yaxis.set_minor_locator(ticker.NullLocator())
    ax.yaxis.set_minor_formatter(ticker.NullFormatter())


def _set_locator_formatter(ax: plt.Axes, axis: Literal["x", "y", "xy", "yx"]) -> None:
    match axis:
        case "x":
            _remove_x_axis_ticks(ax)
        case "y":
            _remove_y_axis_ticks(ax)
        case "xy" | "yx":
            _remove_x_axis_ticks(ax)
            _remove_y_axis_ticks(ax)
        case _:
            raise ValueError("axis must be 'x', 'y' or 'xy' or 'yx'")


def remove_ticks(f: Figure | Axes | Iterable[Axes], axis: Literal["x", "y", "xy", "yx"]) -> None:
    """Removes ticks from the axis of a figure or axis object. If a figure is passed,
    the function will remove the axis-ticks of all the figure's axes.

    Args
    ----------
    f : Figure | Axes | Iterable[Axes]
        The figure or axis object to remove the ticks from.
    axis : Literal["x", "y", "xy", "yx"]
        The axis to remove the ticks from. If "xy" or "yx" is passed, the function will remove
        the ticks from both axes.

    Raises
    ------
    ValueError
        If f is not a Figure or Axis object.
    """
    match f:
        case Figure():
            axes = f.axes
            map(lambda a: _set_locator_formatter(a, axis), axes)
        case plt.Axes():
            _set_locator_formatter(f, axis)
        case Iterable() | list() | np.ndarray():
            assert all(isinstance(a, plt.Axes) for a in f), "f must be an iterable of Axes objects"
            map(lambda a: _set_locator_formatter(a, axis), f)
        case _:
            raise ValueError("f must be a Figure, an Axes object, or a lis tof Axes objects")


def mask_erosion_ufunc(
    x: xr.DataArray,
    connectivity: int = 2,
    mode: Literal["thick", "inner", "outer", "subpixel"] = "thick",
):
    """_summary_

    Parameters
    ----------
    x : xr.DataArray
        The input label image
    connectivity : int, optional
        The connectivity used to find boundaries, by default 2
    mode : Literal["thick", "inner", "outer", "subpixel"], optional
        How to mark the boundaries , by default "thick"

    Returns
    -------
    NDArray
        The mask of the segmentation with eroded boundaries.
    """
    edges = find_boundaries(
        label_img=x, connectivity=connectivity, mode=mode, background=0
    )
    seg_mask = np.where(edges == 0, x, 0)
    return seg_mask