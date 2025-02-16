"""
MIBI_assign_cells_to_metabolic_zones.py

Overview: This script takes in the segmentation amd metabolic zone masks for each sample.
By combining these masks together each cell is assigned as belonging to the hypoxic or IDO
zones. For cells belonging to both zones, the % area of overlap between a cell and the zones
was used a tiebreaker. Results as exported as a .csv file.

"""

# library imports
import os
import pandas as pd
import skimage.io as io
import numpy as np

# Define directories
base_dir = '/Volumes/T7 Shield/MIBI_data/NHP_TB_Cohort/Panel2'
seg_dir = os.path.join(base_dir, 'final_segmentation')
mask_dir = os.path.join(base_dir, 'masks')
sc_data_dir = os.path.join(base_dir)

# Read in the single cell data
sc_data = pd.read_csv(os.path.join(sc_data_dir, 'cohort_cell_table.csv'))

# Define the samples
samples = sc_data['sample'].unique()

# Create column for hypoxic and IDO1 zone area overlaps and assignments
sc_data['glyco_zone'] = 0
sc_data['IDO1_zone'] = 0
sc_data['glyco_area'] = 0
sc_data['IDO1_area'] = 0

# Iterate over all samples
for sample in samples:

    print('Working on: ', sample)

    # read in the segmentation data
    seg_mask = io.imread(os.path.join(seg_dir, sample + '_labels.tiff'))

    # read in the masks
    glyco_mask = io.imread(os.path.join(mask_dir, 'glyco_masks', sample + '_glyco_mask.png'))
    IDO_mask = io.imread(os.path.join(mask_dir, 'IDO_masks', sample + '_IDO_mask.png'))

    # append the area values
    glyco_area = np.sum(glyco_mask.astype(bool))
    IDO_area = np.sum(IDO_mask.astype(bool))
    sc_data.loc[(sc_data['sample'] == sample), 'glyco_area'] = glyco_area
    sc_data.loc[(sc_data['sample'] == sample), 'IDO1_area'] = IDO_area

    # get the labels associated with that fov
    labels = sc_data.loc[sc_data['sample'] == sample, 'tiled_label'].tolist()

    # get the overlap between the cell label mask the zone masks
    labels_glyco_mask = seg_mask * glyco_mask.astype(bool)
    labels_IDO_mask = seg_mask * IDO_mask.astype(bool)

    # extract the labels
    glyco_labels = np.unique(labels_glyco_mask).tolist()
    IDO_labels = np.unique(labels_IDO_mask).tolist()

    # update the assignments
    sc_data.loc[(sc_data['sample'] == sample) & (sc_data['tiled_label'].isin(glyco_labels)), 'glyco_zone'] = 1
    sc_data.loc[(sc_data['sample'] == sample) & (sc_data['tiled_label'].isin(IDO_labels)), 'IDO1_zone'] = 1

    # break the ties
    common_labels = list(set(glyco_labels).intersection(set(IDO_labels)))

    counter = 1

    for common_label in common_labels:

        print("label: ", counter, "out of: ", len(common_labels))
        counter += 1

        # get mask for label
        label_mask = (seg_mask == common_label).astype(int)

        # get the overlap between the cell label mask the zone masks
        label_glyco_mask = label_mask * glyco_mask.astype(bool)
        label_IDO_mask = label_mask * IDO_mask.astype(bool)

        # determine the % area of the cell in the zones
        label_area = np.sum(label_mask)
        label_glyco_area = np.sum(label_glyco_mask)
        label_IDO_area = np.sum(label_IDO_mask)

        perc_glyco = label_glyco_area / label_area
        perc_IDO = label_IDO_area / label_area

        # update labels
        if perc_glyco < perc_IDO:
            sc_data.loc[(sc_data['sample'] == sample) & (sc_data['tiled_label'] == common_label), 'glyco_zone'] = 0
        if perc_IDO < perc_glyco:
            sc_data.loc[(sc_data['sample'] == sample) & (sc_data['tiled_label'] == common_label), 'IDO1_zone'] = 0

# trim off everything except the sample, tiled_label, pheno_corrected, and zone data
export_data = sc_data[['sample', 'tiled_label', 'pheno_corrected', 'glyco_zone', 'IDO1_zone', 'glyco_area', 'IDO1_area']]
export_data.to_csv(os.path.join(sc_data_dir, 'cell_cohort_data_metabolic_zones.csv'))