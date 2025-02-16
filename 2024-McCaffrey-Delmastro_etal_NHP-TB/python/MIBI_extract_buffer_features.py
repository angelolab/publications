import matplotlib.pyplot as plt
import os
import numpy as np
import pandas as pd
import skimage.io as io
from tmi import io_utils
from natsort import natsorted
import matplotlib.pyplot as plt
from skimage.measure import label, regionprops_table


# define directories
base_dir = '/Volumes/T7 Shield/MIBI_data/NHP_TB_Cohort/Panel2'
sc_data_dir = os.path.join(base_dir)
segmentation_dir = os.path.join(base_dir, 'final_segmentation')
buffer_mask_dir = os.path.join(base_dir, 'masks', 'buffers')
cell_mask_dir = os.path.join(base_dir, 'masks', 'cellular_masks')
im_dir = os.path.join(base_dir, 'tile_corrected', 'tiled')

# dead in the sc data
sc_data = pd.read_csv(os.path.join(sc_data_dir, 'cohort_cell_table.csv'))

# define the samples to process
samples = io_utils.list_files(buffer_mask_dir, substrs=['sample'])
cohort_samples = [sub.replace('_buffers.png', '') for sub in samples]
cohort_samples = natsorted(cohort_samples)

# define cell phenotypes
cell_phenos = np.unique(sc_data['pheno_corrected']).tolist()

# define channels
channel_tiffs = io_utils.list_files(os.path.join(im_dir,cohort_samples[0]), substrs=['.tiff'])
remove = ['ASCT2.tiff', 'Au.tiff', 'HIF1a.tiff', 'HH3.tiff', 'Noodle.tiff', 'Na.tiff', 'PD-L1.tiff', 'Ta.tiff']
channel_tiffs_filtered = list(filter(lambda i: i not in remove, channel_tiffs))
channels = [sub.replace('.tiff', '') for sub in channel_tiffs_filtered]
# sc_channels = [sub.replace('.tiff', '_sc') for sub in channel_tiffs_filtered]

# create empty output
output_columns = ['sample', 'total_buffer_num', 'buffer_num', 'buffer_area', 'buffer_cell_num'] + cell_phenos + channels
buffer_data = pd.DataFrame(columns=output_columns)

# create counter for easy row indexing
row_idx = 0

# iterate through samples
for sample in cohort_samples:

    print("Working on: " + sample)

    # read in the necessary masks
    buffer_mask = io.imread(os.path.join(buffer_mask_dir, sample + '_buffers.png'))
    cellular_mask = io.imread(os.path.join(cell_mask_dir, sample + '_cellular_mask.png'))
    seg_mask = io.imread(os.path.join(segmentation_dir, sample + '_labels.tiff'))

    # generate region props for segmentation
    label_props = pd.DataFrame(regionprops_table(seg_mask.astype(int), properties=('label', 'area')))

    # process buffer mask to only include cellular area
    buffer_mask_cellular = buffer_mask * (cellular_mask.astype(int))

    # determine total number of buffers
    buffer_num = np.max(buffer_mask)

    # process buffers
    for buffer in range(1, buffer_num+1):

        print("Buffer " + str(buffer) + " out of " + str(buffer_num))

        # create df row
        buffer_row = {'sample': sample, 'total_buffer_num': buffer_num, 'buffer_num': buffer}
        # buffer_data = pd.concat([buffer_data, buffer_row], ignore_index=True)
        buffer_data = buffer_data.append(buffer_row, ignore_index=True)

        # get a mask of just the buffer
        buffer_sub_mask = (buffer_mask == buffer).astype(int)

        # get the buffer area and append
        buffer_area = np.sum(buffer_sub_mask)
        buffer_data.at[row_idx, 'buffer_area'] = buffer_area

        # cell data #

        # get the cells that are in that buffer
        cell_buffer_mask = seg_mask * buffer_sub_mask

        # generate region props on the buffer-overlapped labels
        buffer_label_props = pd.DataFrame(regionprops_table(cell_buffer_mask.astype(int), properties=('label', 'area')))

        # subset the segmentation props on the overlapped labels
        buffer_labels = np.unique(cell_buffer_mask).astype(int).tolist()
        label_props_sub = label_props[label_props['label'].isin(buffer_labels)]

        # combine and get percentage overlap
        label_props_sub = label_props_sub.merge(buffer_label_props, on='label')
        label_props_sub['overlap_freq'] = label_props_sub['area_y'] / label_props_sub['area_x']

        # extract labels for those with at least 50% overlap
        overlap_labels = label_props_sub.loc[label_props_sub['overlap_freq'] >= 0.5, 'label'].tolist()

        # get the cell phenotypes and their frequency
        overlap_phenos = sc_data.loc[(sc_data['sample'] == sample) & (sc_data['tiled_label'].isin(overlap_labels)), 'pheno_corrected']
        overlap_pheno_counts = overlap_phenos.value_counts()

        # append information to df
        buffer_data.at[row_idx, 'buffer_cell_num'] = len(overlap_labels)
        for pheno, count in overlap_pheno_counts.items():
            buffer_data.at[row_idx, pheno] = count

        # channels #
        for channel in channels:

            # read in the channel data
            channel_data = io.imread(os.path.join(im_dir, sample, channel + '.tiff'))

            # subset to get the masked portion
            channel_data_buffer = channel_data * buffer_sub_mask

            # append to df
            buffer_data.at[row_idx, channel] = np.sum(channel_data_buffer)

        # update the counter
        row_idx = row_idx + 1

# save=
buffer_data.to_csv(os.path.join(base_dir, 'buffer_data.csv'))