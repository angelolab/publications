import numpy as np
import skimage.io as io
import os
from tmi import io_utils
import matplotlib.pyplot as plt

# define directories
base_dir = '/Volumes/T7 Shield/MIBI_data/NHP_TB_Panel1/Panel1'
mask_dir = os.path.join(base_dir, 'masks')

# define samples
# samples = io_utils.list_folders(os.path.join(base_dir, 'normalized/per_sample'))
samples = ['sample49', 'sample57']

# create directory for masks
output_dir = os.path.join(mask_dir, 'cellular_masks')
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

for sample in samples:

    print('Working on ' + sample)

    # read in cellular and necrosis masks
    au_mask = io.imread(os.path.join(mask_dir, 'au_masks', sample + '_gold_mask.png')).astype(np.uint8)
    necrosis_mask = io.imread(os.path.join(mask_dir, 'necrosis_masks', sample, sample + '_Panel1_Necrosis.tif')).astype(np.uint8)

    # rescale 0-1
    necrosis_mask[np.where(necrosis_mask > 1)] = 1

    # combine
    combo_mask = (au_mask + necrosis_mask)
    combo_mask[np.where(combo_mask > 1)] = 1
    combo_mask = combo_mask.astype(bool)

    # invert
    cellular_mask = np.logical_not(combo_mask)

    # save
    io.imsave(os.path.join(output_dir, sample + '_cellular_mask.png'), cellular_mask.astype(np.uint8),
              check_contrast=False)