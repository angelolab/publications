import numpy as np
import skimage.io as io
import os
import erin_utils
import pandas as pd
from skimage import morphology
from tmi import io_utils
import matplotlib.pyplot as plt
from skimage.measure import label, regionprops_table

# define directories
base_dir = '/Volumes/T7 Shield/MIBI_data/NHP_TB_Cohort/Panel2'
# im_dir = os.path.join(base_dir, 'tile_corrected/tiled')
# im_dir = '/Volumes/T7 Shield/MIBI_data/NHP_TB_Cohort/Panel2/GIS_inputs/glycozone/unpadded'
im_dir = base_dir

output_dir = os.path.join(base_dir, 'masks')
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# define samples
samples = io_utils.list_files('/Volumes/T7 Shield/MIBI_data/NHP_TB_Cohort/Panel2/masks/polygons', substrs=['poly_padded'])
cohort_samples = [sub.replace('_poly_padded.tif', '') for sub in samples]

for sample in cohort_samples:

    print('Working on ' + sample)

    # create mask for slide background
    im_IDO1 = io.imread(os.path.join(im_dir, 'tile_corrected/tiled', sample, 'IDO1.tiff'))
    # im_GLUT1 = io.imread(os.path.join(im_dir, 'GIS_inputs/glycozone/unpadded', sample + '_GLUT1-IDOneg.tiff'))

    disk = morphology.disk(10)  # creates a disk of 1 with radius = 10

    # IDO1 mask
    counts, bins = np.histogram(im_IDO1, bins=100)
    plt.stairs(counts, bins)

    mask_IDO1 = erin_utils.create_channel_mask(img=im_IDO1, intensity_thresh=0.00005, sigma=5,
                                          min_mask_size=50000, max_hole_size=1000)
    mask_IDO1 = morphology.binary_closing(mask_IDO1, disk)
    plt.imshow(mask_IDO1)

    # evaluate size and remove objects to get clean mask
    # IDO1_label_img = label(mask_IDO1)
    # plt.imshow(IDO1_label_img, 'jet')
    # plt.colorbar()
    # IDO1_props = pd.DataFrame(regionprops_table(IDO1_label_img, properties=('label', 'area')))
    # IDO1_props
    # IDO1_remove = [1,2]
    # for remove_label in IDO1_remove:
    #     IDO1_label_img[np.where(IDO1_label_img == remove_label)] = 0
    # plt.imshow(IDO1_label_img, 'jet')
    # IDO1_binarize = IDO1_label_img > 0
    # plt.imshow(IDO1_binarize)

    io.imsave(os.path.join(output_dir, 'IDO_masks', sample + '_IDO_mask.png'), mask_IDO1, check_contrast=False)

    # # GLUT1 mask
    # counts, bins = np.histogram(im_GLUT1, bins=100)
    # plt.stairs(counts, bins)
    #
    # mask_GLUT1 = erin_utils.create_channel_mask(img=im_GLUT1, intensity_thresh=0.00003, sigma=5,
    #                                             min_mask_size=30000, max_hole_size=1000)
    # plt.imshow(mask_GLUT1)
    #
    # mask_GLUT1 = morphology.binary_closing(mask_GLUT1, disk)
    # plt.imshow(mask_GLUT1)
    #
    # for i in range(1):
    #     mask_GLUT1 = morphology.binary_erosion(mask_GLUT1)
    # plt.imshow(mask_GLUT1)
    #
    # mask_GLUT1 = morphology.remove_small_objects(mask_GLUT1, min_size=50000)
    # plt.imshow(mask_GLUT1)
    #
    # # evaluate size and remove objects to get clean mask
    # GLUT1_label_img = label(mask_GLUT1)
    # plt.imshow(GLUT1_label_img, 'jet')
    # plt.colorbar()
    # GLUT1_props = pd.DataFrame(regionprops_table(GLUT1_label_img, properties=('label', 'area')))
    # GLUT1_props
    # GLUT1_remove = [6,8,9,10,11,12,14,15,16]
    # for remove_label in GLUT1_remove:
    #     GLUT1_label_img[np.where(GLUT1_label_img == remove_label)] = 0
    # plt.imshow(GLUT1_label_img, 'jet')
    # GLUT1_binarize = GLUT1_label_img > 0
    # plt.imshow(GLUT1_binarize)
    #
    # io.imsave(os.path.join(output_dir, 'glyco_masks', sample + '_glyco_mask.png'), GLUT1_binarize, check_contrast=False)

    output_dir = '/Users/erinmccaffrey/Library/CloudStorage/GoogleDrive-erinmcc@stanford.edu/My ' \
                 'Drive/Grad_School/AngeloLab/MIBIProjects/NHP_TB/Manuscript/Figures/Supplemental_Figures/FigureS4/content'
    io.imsave(os.path.join(output_dir, 'sample53_mask1.png'), mask_IDO1, check_contrast=False)

