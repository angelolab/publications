import numpy as np
import skimage.io as io
import matplotlib.pyplot as plt
import os
import skimage

import numpy as np
import pandas as pd

import skimage.measure
from skimage import morphology
from scipy.ndimage import gaussian_filter

from alpineer.io_utils import list_folders
from python_files import utils

# This script creates image masks defining the tumor compartments and slide background to be
# be used in subsequent feature extraction pipeline

# set up paths
base_dir = '/Volumes/Shared/Noah Greenwald/TONIC_Cohort'
channel_dir = os.path.join(base_dir, 'image_data/samples/')
seg_dir = os.path.join(base_dir,'segmentation_data/deepcell_output')
mask_dir = os.path.join(base_dir, 'intermediate_files/mask_dir/')
analysis_dir = os.path.join(base_dir,'analysis_files')


cell_table_clusters = pd.read_csv(os.path.join(analysis_dir, 'cell_table_clusters.csv'))
folders = list_folders(channel_dir)

# create directories to hold masks
intermediate_dir = os.path.join(mask_dir, 'intermediate_masks')
if not os.path.exists(intermediate_dir):
    os.makedirs(intermediate_dir)

individual_dir = os.path.join(mask_dir, 'individual_masks')
if not os.path.exists(individual_dir):
    os.makedirs(individual_dir)

# loop over each FOV and generate the appropriate masks
for folder in folders:

    ecad = io.imread(os.path.join(channel_dir, folder, 'ECAD.tiff'))

    intermediate_folder = os.path.join(intermediate_dir, folder)
    if not os.path.exists(intermediate_folder):
        os.mkdir(intermediate_folder)

    # generate cancer/stroma mask by combining segmentation mask with ECAD channel
    seg_label = io.imread(os.path.join(seg_dir, folder + '_whole_cell.tiff'))[0]
    seg_mask = utils.create_cell_mask(seg_label, cell_table_clusters, folder, ['Cancer'])
    cancer_mask = utils.create_cancer_boundary(ecad, seg_mask, min_mask_size=7000)
    cancer_mask = cancer_mask.astype(np.uint8)
    io.imsave(os.path.join(intermediate_folder, 'cancer_mask.png'), cancer_mask,
              check_contrast=False)

    # create mask for TLS
    tls_mask = utils.create_cell_mask(seg_label, cell_table_clusters, folder, ['B', 'T'], sigma=4)
    tls_label_mask = skimage.measure.label(tls_mask)
    tls_label_mask = morphology.remove_small_objects(tls_label_mask, min_size=25000)
    tls_label_mask = morphology.remove_small_holes(tls_label_mask, area_threshold=7000)
    tls_label_mask = skimage.measure.label(tls_label_mask)

    # get location of all T cells in image
    cell_subset = cell_table_clusters[cell_table_clusters['fov'] == folder]
    cell_subset = cell_subset[cell_subset['cell_cluster_broad'] == 'T']
    cell_labels = cell_subset['label'].values
    t_mask = np.isin(seg_label, cell_labels)

    # get location of all B cells in image
    cell_subset = cell_table_clusters[cell_table_clusters['fov'] == folder]
    cell_subset = cell_subset[cell_subset['cell_cluster_broad'] == 'B']
    cell_labels = cell_subset['label'].values
    b_mask = np.isin(seg_label, cell_labels)

    # create mask to hold T only aggregates
    tagg_label_mask = np.zeros(tls_label_mask.shape)

    # figure out which TLS objects contain both B and T cells
    for i in range(1, tls_label_mask.max() + 1):
        tls_object = tls_label_mask == i

        # get proportion of pixels that are B cells
        b_pixels = np.sum(b_mask[tls_object])
        t_pixels = np.sum(t_mask[tls_object])
        b_prop = b_pixels / (b_pixels + t_pixels)

        if b_prop < 0.2:
            # not a tls, remove from tls mask and add to tagg
            tls_label_mask[tls_object] = 0
            tagg_label_mask[tls_object] = 1

    # convert to binary
    tls_label_mask = tls_label_mask > 0
    io.imsave(os.path.join(intermediate_folder, 'tls_mask.png'), tls_label_mask.astype(np.uint8),
              check_contrast=False)

    io.imsave(os.path.join(intermediate_folder, 'tagg_mask.png'), tagg_label_mask.astype(np.uint8),
                check_contrast=False)

    # create mask for slide background
    gold = io.imread(os.path.join(channel_dir, folder, 'Au.tiff'))

    gold_mask = utils.create_channel_mask(img=gold, sigma=2, intensity_thresh=350,
                                          min_mask_size=5000, max_hole_size=1000)

    # erode edges of gold mask so that it doesn't encroach on other masks
    for _ in range(15):
        gold_mask = morphology.binary_erosion(gold_mask)

    # any cell can't be in the gold mask
    gold_mask[seg_label > 0] = 0
    io.imsave(os.path.join(intermediate_folder, 'gold_mask.png'), gold_mask.astype(np.uint8),
                check_contrast=False)


# remove any overlapping pixels from different masks, then save individually
for folder in folders:
    # read in generated masks
    intermediate_folder = os.path.join(intermediate_dir, folder)
    cancer_mask = io.imread(os.path.join(intermediate_folder, 'cancer_mask.png'))
    gold_mask = io.imread(os.path.join(intermediate_folder, 'gold_mask.png'))
    tls_mask = io.imread(os.path.join(intermediate_folder, 'tls_mask.png'))
    tagg_mask = io.imread(os.path.join(intermediate_folder, 'tagg_mask.png'))

    # create a single unified mask; TLS and background override tumor compartments
    cancer_mask[gold_mask == 1] = 0
    cancer_mask[tls_mask == 1] = 5
    cancer_mask[tagg_mask == 1] = 6

    # save individual masks
    processed_folder = os.path.join(individual_dir, folder)
    if not os.path.exists(processed_folder):
        os.mkdir(processed_folder)

    for idx, name in zip(range(0, 7), ['empty_slide', 'stroma_core', 'stroma_border',
                                       'cancer_border', 'cancer_core', 'tls', 'tagg']):
        channel_img = cancer_mask == idx
        io.imsave(os.path.join(processed_folder, name + '.tiff'), channel_img.astype(np.uint8),
                  check_contrast=False)

# compute the area of each mask
area_df = utils.calculate_mask_areas(mask_dir=individual_dir, fovs=folders)

# combine tls and tagg masks into single immune_agg compartment
for fov in np.unique(area_df.fov):
    fov_df = area_df[area_df.fov == fov]
    tls_tagg_sum = fov_df[fov_df.compartment == 'tls'].area.values[0] + fov_df[fov_df.compartment == 'tagg'].area.values[0]
    area_df = pd.concat([pd.DataFrame([['immune_agg', tls_tagg_sum, fov]], columns=area_df.columns), area_df], ignore_index=True)
area_df.to_csv(os.path.join(mask_dir, 'fov_annotation_mask_area.csv'), index=False)

# create combined images for visualization
for folder in folders[:20]:
    cluster_overlay = io.imread(os.path.join(base_dir, 'overlay_dir/cell_cluster_overlay', folder + '.png'))
    compartment_overlay = io.imread(os.path.join(base_dir, 'overlay_dir/compartment_overlay', folder + '.png'))
    gold_chan = io.imread(os.path.join(channel_dir, folder, 'Au.tiff'))
    border_mask = io.imread(os.path.join(intermediate_dir, folder, 'cancer_mask.png'))
    gold_mask = io.imread(os.path.join(intermediate_dir, folder, 'gold_mask.png'))
    tls_mask = io.imread(os.path.join(intermediate_dir, folder, 'tls_mask.png'))
    tagg_mask = io.imread(os.path.join(intermediate_dir, folder, 'tagg_mask.png'))

    # create a single unified mask; TLS and background override tumor compartments
    border_mask[gold_mask == 1] = 0
    border_mask[tls_mask == 1] = 5
    border_mask[tagg_mask == 1] = 6

    # make top row shorter than bottom row
    fig, ax = plt.subplots(2, 2, figsize=(15, 10), gridspec_kw={'height_ratios': [1, 2]})
    ax[1, 0].imshow(cluster_overlay)
    ax[1, 0].axis('off')
    ax[0, 0].imshow(gold_chan)
    ax[0, 1].imshow(border_mask)
    ax[1, 1].imshow(compartment_overlay)
    ax[1, 1].axis('off')

    plt.tight_layout()
    plt.savefig(os.path.join(base_dir, 'overlay_dir/combined_mask_overlay', folder + '.png'))
    plt.close()


# assign cells to the correct compartment
all_assignment_table = pd.DataFrame()
for i in range(0, 1400, 200):
    assignment_table = utils.assign_cells_to_mask(seg_dir=seg_dir, mask_dir=individual_dir, fovs=folders[i:i+200])
    # assignment_table.to_csv(os.path.join(mask_dir, 'annotation_files', 'cell_annotation_mask_{}'.format(i)), index=False)
    # assignment_table = pd.read_csv(os.path.join(mask_dir, 'annotation_files', 'cell_annotation_mask_{}'.format(i)))
    all_assignment_table = pd.concat([all_assignment_table, assignment_table])

# replace tls and tagg assignments with immune_agg
all_assignment_table = all_assignment_table.replace({'tls': 'immune_agg', 'tagg': 'immune_agg'})

all_assignment_table.to_csv(os.path.join(mask_dir, 'cell_annotation_mask.csv'), index=False)
