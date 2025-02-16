"""This script takes the buffer zone and metabolic masks for each FOV and determines the
proportion of buffer area that overlaps with the zone. This information is extracetd as
a csv."""

import os
import pandas as pd
from tmi import io_utils
from natsort import natsorted
import skimage.io as io
import matplotlib.pyplot as plt
import numpy as np

# Define directories
base_dir = '/Volumes/T7 Shield/MIBI_data/NHP_TB_Cohort/Panel2'
mask_dir = os.path.join(base_dir, 'masks')
buffer_mask_dir = os.path.join(mask_dir, 'buffers')
cell_mask_dir = os.path.join(mask_dir, 'cellular_masks')
glyco_mask_dir = os.path.join(mask_dir, 'glyco_masks')
IDO_mask_dir = os.path.join(mask_dir, 'IDO_masks')

# Create an empty output
output_columns = ['sample', 'total_buffer_num', 'buffer_num', 'buffer_area', 'buffer_glyco_area', 'buffer_IDO1_area']
buffer_data = pd.DataFrame(columns=output_columns)

# get the sample inform
samples = io_utils.list_files(buffer_mask_dir, substrs=['sample'])
cohort_samples = [sub.replace('_buffers.png', '') for sub in samples]
cohort_samples = natsorted(cohort_samples)

# create counter for easy row indexing
row_idx = 0

# Iterate
for sample in cohort_samples:

    print('Working on: ', sample)

    # read in the buffer data
    buffer_mask = io.imread(os.path.join(buffer_mask_dir, sample + '_buffers.png'))

    # get the cellular mask
    cellular_mask = io.imread(os.path.join(cell_mask_dir, sample + '_cellular_mask.png'))

    # read in the metabolic masks
    glyco_mask = io.imread(os.path.join(glyco_mask_dir, sample + '_glyco_mask.png'))
    IDO_mask = io.imread(os.path.join(IDO_mask_dir, sample + '_IDO_mask.png'))

    # process buffer mask to only include cellular area
    buffer_mask_cellular = buffer_mask * (cellular_mask.astype(int))

    # determine total number of buffers
    buffer_num = np.max(buffer_mask)

    # process buffers
    for buffer in range(1, buffer_num+1):

        print("Buffer " + str(buffer) + " out of " + str(buffer_num))

        # create df row
        buffer_row = {'sample': sample, 'total_buffer_num': buffer_num, 'buffer_num': buffer}
        buffer_data = buffer_data.append(buffer_row, ignore_index=True)

        # get a mask of just the buffer
        buffer_sub_mask = (buffer_mask == buffer).astype(int)

        # get the buffer area and append
        buffer_area = np.sum(buffer_sub_mask)
        buffer_data.at[row_idx, 'buffer_area'] = buffer_area

        # get the metabolic masks for the buffer
        glyco_buffer_mask = glyco_mask.astype(bool) * buffer_sub_mask
        IDO_buffer_mask = IDO_mask.astype(bool) * buffer_sub_mask

        # append the metabolic areas
        glyco_buffer_area = np.sum(glyco_buffer_mask)
        buffer_data.at[row_idx, 'buffer_glyco_area'] = glyco_buffer_area
        IDO_buffer_area = np.sum(IDO_buffer_mask)
        buffer_data.at[row_idx, 'buffer_IDO1_area'] = IDO_buffer_area

        # update the counter
        row_idx = row_idx + 1

buffer_data.to_csv(os.path.join(base_dir, 'metabolic_buffer_data.csv'))
