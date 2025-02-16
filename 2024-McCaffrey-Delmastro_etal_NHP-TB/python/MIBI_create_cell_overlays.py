"""Produces segmentation mask colored by the cell phenotype"""

import erin_utils as eu
import os
import pandas as pd
from tmi import io_utils

# Define directories
base_dir = '/Volumes/T7 Shield/MIBI_data/NHP_TB_Cohort/Panel2'
sc_data_dir = os.path.join(base_dir)
segmentation_dir = os.path.join(base_dir, 'final_segmentation')

# Create output directory
# overlay_dir = os.path.join(base_dir)
overlay_dir = "/Users/erinmccaffrey/Library/CloudStorage/GoogleDrive-erinmcc@stanford.edu/My " \
              "Drive/Grad_School/AngeloLab/MIBIProjects/NHP_TB/Manuscript/Figures/Figure4"
if not os.path.exists(overlay_dir):
    os.mkdir(overlay_dir)

# read in sc data
sc_data = pd.read_csv(os.path.join(sc_data_dir, 'cohort_cell_table.csv'))

# define column with cluster
cluster_col = 'cell_pheno_num'

# read in color key and extract values
color_key = pd.read_csv(os.path.join(base_dir, 'keys/cell_color_key_custom.csv'))
color_key = color_key.sort_values(by='Num')
cluster_vals = color_key['Num'].values.tolist()
colors = color_key['Hex'].values.tolist()

# define fovs
# samples = sc_data['sample'].unique()
# get list of samples to process
# samples = io_utils.list_folders(os.path.join(base_dir, 'exports'), substrs=['Series'])
samples = ['sample33','sample40']

# create overlays
eu.create_cell_overlay(cell_table=sc_data,
                       seg_folder=segmentation_dir,
                       fovs=samples,
                       cluster_col=cluster_col,
                       plot_dir=overlay_dir,
                       cluster_vals=cluster_vals,
                       color_values=colors)