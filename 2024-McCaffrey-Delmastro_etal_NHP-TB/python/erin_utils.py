import pandas as pd
import numpy as np
from skimage.segmentation import find_boundaries
from matplotlib import colors
import skimage.io as io
from matplotlib.colors import ListedColormap, BoundaryNorm
from skimage import morphology
from scipy.ndimage import gaussian_filter
from skimage.measure import label
from matplotlib import cm
from tmi import io_utils, data_utils
import matplotlib.pyplot as plt
import os


def crop_mask(mask, fov, im_size):
    """Takes in a necrosis mask and a fov.

    Returns a cropped necrosis mask.

    Args:
        mask (ndarray):
            Mask for the tiled image.
        fov (str):
            Current fov name.
        im_size (double):
            Size of each tile.

    Returns:
        cropped_mask (ndarray):
            The portion of the mask for the given fov.
    """

    # get x and y positions
    row_pos = int(fov.split('row')[1].split('_')[0])
    col_pos = int(fov.split('col')[1])

    # get range for crop
    rStart = im_size * row_pos
    rEnd = (im_size * row_pos) + im_size

    cStart = im_size * col_pos
    cEnd = (im_size * col_pos) + im_size

    # crop mask
    cropped_mask = mask[rStart:rEnd, cStart:cEnd]

    return cropped_mask


def get_mean_correction_signal(data_dir, segmentation_dir, sample_fovs, correction_channels_nec,
                               correction_channels_cell, necrosis_mask, cellular_mask, seg_suffix):
    """Takes a set of fovs and returns a vector of the mean nuclear correction signal in each.

    Returns a list of values corresponding the mean signal per fov.

    Args:
        data_dir (str):
            Location of correction channel data.
        segmentation_dir (str):
            Location of the segmentation masks.
        sample_fovs (list):
            A list of the fovs to be corrected.
        correction_channels (list):
            The image channels to be averaged.
        necrosis_mask (ndarray):
            Mask for the necrosis zone of the tiled image.
        cellular_mask (ndarray):
            Mask for the cellular zone of the tiled image.
        seg_suffix (str):
            The suffix at the end of the nuclear segmentation mask.

    Returns:
        mean_cell_values (list):
            A list of values corresponding to the mean correction channel signal in each fov for cellular zones.
        mean_necrotic_values (list):
            A list of values corresponding to the mean correction channel signal for necrotic zones.
    """

    # create empty list
    mean_cell_values = []
    mean_necrosis_values = []

    for fov in sample_fovs:

        correction_data_nec = []
        correction_data_cell = []

        # read in channel data
        for i in range(len(correction_channels_nec)):
            channel_data_nec = io.imread(os.path.join(data_dir, fov, correction_channels_nec[i]))
            correction_data_nec.append(channel_data_nec)
        correction_data_nec = sum(correction_data_nec)

        for i in range(len(correction_channels_cell)):
            channel_data_cell = io.imread(os.path.join(data_dir, fov, correction_channels_cell[i]))
            correction_data_cell.append(channel_data_cell)
        correction_data_cell = sum(correction_data_cell)

        # read in segmentation data, convert to binary mask
        seg = io.imread(os.path.join(segmentation_dir, fov + seg_suffix))[0, :, :]
        seg_mask = seg > 0

        # transform the correction data and seg mask by the necrosis mask

        # get the necrosis mask and data
        necrosis_crop = crop_mask(necrosis_mask, fov, im_size=1024)
        correction_data_nec = correction_data_nec * necrosis_crop

        # get the cell_mask and data
        cell_mask = crop_mask(cellular_mask, fov, im_size=1024)
        correction_data_cell = correction_data_cell * cell_mask
        seg_mask_cell = seg_mask * cell_mask

        # subset by segmentation mask
        channel_data_nec = correction_data_nec
        channel_data_cell = correction_data_cell * seg_mask_cell

        # get mean values for necrotic and non-necrotic
        # get a vector of the non-zero values
        channel_values_nec = channel_data_nec[np.nonzero(channel_data_nec)]
        channel_values_cell = channel_data_cell[np.nonzero(channel_data_cell)]

        # get the mean value and append to list
        mean_channel_value_nec = np.mean(channel_values_nec)
        mean_necrosis_values.append(mean_channel_value_nec)

        mean_channel_value_cell = np.mean(channel_values_cell)
        mean_cell_values.append(mean_channel_value_cell)

    return mean_cell_values, mean_necrosis_values


# from ark data_utils
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

    default_label = max(labels_dict.values()) + 1

    # cast to int16 to allow for Photoshop loading
    relabeled_img = np.vectorize(
        lambda x: labels_dict.get(x, default_label) if x != 0 else 0)(img).astype("int16")

    return relabeled_img


# Noah's cell overlay function with modifications
def create_cell_overlay(cell_table,
                        seg_folder,
                        fovs,
                        cluster_col,
                        plot_dir,
                        cluster_vals,
                        color_values):
    """Takes in the information on segmentation and phenotype to produce a CPM.

    Saves the segmentation mask colored by the variable encoded in the indicated cluster col.

    Args:
        cell_table (pandas.array):
            pandas array for annotated single cell data.
        seg_folder (str):
            Location of the segmentation masks.
        fovs (list):
            List of strings indicating the fovs to be plotted.
        cluster_col (str):
            Indicates the column in cell_table that provides the annotation information.
        plot_dir (str):
            Location to save the plotted overlays.
        cluster_vals (str):
            Unique labels in cluster_col.
        color_values (str):
            Unique values to map to cluster_vals.
    """

    # create cmap
    new_cmap = ListedColormap(color_values)
    # bounds = np.append(cluster_vals, max(cluster_vals) + 1)
    # norm = BoundaryNorm(bounds, ncolors=len(color_values))

    for image in fovs:
        seg_mask = io.imread(os.path.join(seg_folder, image + '_labels.tiff'))

        edges = find_boundaries(seg_mask, mode='inner')
        seg_mask = np.where(edges == 0, seg_mask, 0)

        # convert string entries in pandas df to unique integers
        cell_subset_plot = cell_table[cell_table['sample'] == image]
        labels_dict = dict(zip(cell_subset_plot['tiled_label'], cell_subset_plot[cluster_col]))

        # relabel the array
        relabeled_img_array = relabel_segmentation(seg_mask, labels_dict)

        # adjust color map if necessary
        if np.max(relabeled_img_array) < len(color_values) - 1:
            max_cluster = np.max(relabeled_img_array)
            color_values_adj = color_values[0:max_cluster+1]
            sample_cmap = ListedColormap(color_values_adj)
        else:
            sample_cmap = new_cmap

        # im = plt.imshow(relabeled_img_array, cmap=new_cmap, norm=norm)
        # tick_names = ['Empty'] + cluster_vals
        # cbar = plt.colorbar(im, ticks=np.arange(len(tick_names)))
        # cbar.set_ticks(cbar.ax.get_yticks())
        # cbar.ax.set_yticklabels(tick_names)
        # # plt.savefig(os.path.join(plot_dir, save_names[idx]), dpi=300)
        # plt.close()

        plt.imsave(os.path.join(plot_dir, image + '_cell_overlay.png'), relabeled_img_array, cmap=sample_cmap)


def make_corruption_mask(sample_key, im_x, im_y):
    """
    Produces a binary mask to map where there are corrupted_fovs in a sample.

    Args:
        sample_key (pd.array):
            Pandas array providing information of the location of corrupted fovs.
        im_x (int):
            Number of pixels in the x-plane of image.
        im_y (int):
            Number of pixels in the y-plane of image.

    Returns:
        corruption_mask (nd.array):
        Binary mask where value of 1 indicates corruption.
    """
    # get the positions from the mapping
    pos_names = sample_key['position'].values.tolist()

    # get rows
    row_names = [pos for pos_name in pos_names for pos in pos_name.split('row')[1].split('_')[0]]
    row_values = list(map(int, row_names))
    row_num = max(row_values) + 1

    # get cols
    col_names = [pos for pos_name in pos_names for pos in pos_name.split('col')[1]]
    col_values = list(map(int, col_names))
    col_num = max(col_values) + 1

    # get step sizes
    step_x = int(im_x / row_num)
    step_y = int(im_y / col_num)

    # create empty numpy array to make mask
    corruption_mask = np.zeros([im_x, im_y])

    for i in range(0, row_num):
        # define row range
        rStart = step_x * i
        rEnd = (step_x * i) + (step_x)
        for j in range(0, col_num):
            # define col range
            cStart = step_y * j
            cEnd = (step_y * j) + (step_y)
            # get position to reference fov
            position = ('row' + str(i) + "_col" + str(j))
            # print('Current position: ' + position)
            # check if corrupted
            if sample_key[sample_key['position'] == position]['corrupted'].values[0] == 1:
                im = np.ones((1024, 1024))
            else:
                im = np.zeros((1024, 1024))
            # populate mosaic
            corruption_mask[rStart:rEnd, cStart:cEnd] = im[:, :]

    return corruption_mask


def adjust_mask(sample_key, mask):
    """
    Adjusts mask to account for areas of corrupted fovs.

    Args:
        sample (str):
            Sample to be charted.
        sample_key (pd.array):
            Pandas array providing information of the location of corrupted fovs.
        mask (nd.array):
            Unadjusted mask
    Returns:
        adjusted_mask (nd.array):
        Mask where regions of corruption are zero-ed out.
    """
    # get the x and y size of the image
    im_x = mask.shape[0]
    im_y = mask.shape[1]

    # get corruption mask
    corruption_mask = make_corruption_mask(sample_key, im_x, im_y)

    # invert corruption mask, multiply by input mask
    inverted_mask = (~corruption_mask.astype(bool)).astype(int)

    # adjust mask
    adjusted_mask = mask * inverted_mask

    return adjusted_mask


def create_channel_mask(img, intensity_thresh, sigma, min_mask_size=0, max_hole_size=100000):
    """Generates a binary mask from a single channel image (from Noah)

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
