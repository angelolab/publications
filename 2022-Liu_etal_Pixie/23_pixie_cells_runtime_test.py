"""
Runtime test for Pixie
Author: Candace Liu
Date: 1/12/23

"""

import shutil
import time
import random

from datetime import datetime as dt
import os
import subprocess

import feather
import json
from matplotlib import rc_file_defaults
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
import scipy.stats as stats
import seaborn as sns

from ark.analysis import visualize
from ark.phenotyping import pixel_cluster_utils, cell_cluster_utils
from ark.utils import data_utils, example_dataset, io_utils, load_utils, plot_utils
from ark.utils.metacluster_remap_gui import colormap_helper, MetaClusterData, MetaClusterGui, metaclusterdata_from_files

random.seed(10)
dataset_sizes = [10,50,100,165]

# pixel clustering
base_dir = "cimac_data"
tiff_dir = os.path.join(base_dir, "input_data", "single_channel_inputs")
img_sub_folder = None
segmentation_dir = "deepcell_output"
seg_suffix = "_feature_0.tif"
all_fovs = io_utils.list_folders(tiff_dir)
multiprocess = False
batch_size = 5

pixel_cluster_prefix = "example"
pixel_output_dir = os.path.join("pixie", "%s_pixel_output_dir" % pixel_cluster_prefix)

pixel_data_dir = os.path.join(pixel_output_dir, 'pixel_mat_data')
pixel_subset_dir = os.path.join(pixel_output_dir, 'pixel_mat_subset')
norm_vals_name = os.path.join(pixel_output_dir, 'channel_norm_post_rowsum.feather') 
 
channels = ["beta-tubulin","CD11c","CD20","CD31",
            "CD3","CD45","CD4","CD56","CD68",
            "CD8","dsDNA","PANCK","PAX5","Vimentin"]
blur_factor = 2

pixel_som_weights_name = os.path.join(pixel_output_dir, 'pixel_som_weights.feather')
pc_chan_avg_som_cluster_name = os.path.join(pixel_output_dir, 'pixel_channel_avg_som_cluster.csv')
pc_chan_avg_meta_cluster_name = os.path.join(pixel_output_dir, 'pixel_channel_avg_meta_cluster.csv')
pixel_meta_cluster_remap_name = os.path.join(pixel_output_dir, 'pixel_meta_cluster_mapping.csv')

max_k = 20
cap = 3

# cell clustering
cell_table_path = os.path.join(base_dir, "single_cell_output", "cell_table_size_normalized.csv")

cell_cluster_prefix = "example"
cell_output_dir = '%s_cell_output_dir' % cell_cluster_prefix

cell_som_weights_name = os.path.join("pixie", cell_output_dir, 'cell_som_weights.feather')
cluster_counts_name = os.path.join("pixie", cell_output_dir, 'cluster_counts.feather')
cluster_counts_size_norm_name = os.path.join("pixie", cell_output_dir, 'cluster_counts_size_norm.feather')
weighted_cell_channel_name = os.path.join("pixie", cell_output_dir, 'weighted_cell_channel.feather')
cell_som_cluster_count_avg_name = os.path.join("pixie", cell_output_dir, 'cell_som_cluster_count_avg.csv')
cell_meta_cluster_count_avg_name = os.path.join("pixie", cell_output_dir, 'cell_meta_cluster_count_avg.csv')
cell_som_cluster_channel_avg_name = os.path.join("pixie", cell_output_dir, 'cell_som_cluster_channel_avg.csv')
cell_meta_cluster_channel_avg_name = os.path.join("pixie", cell_output_dir, 'cell_meta_cluster_channel_avg.csv')
cell_meta_cluster_remap_name = os.path.join("pixie", cell_output_dir, 'cell_meta_cluster_mapping.csv')

pixel_cluster_col = 'pixel_meta_cluster_rename'

if pixel_cluster_col == 'pixel_som_cluster':
    pc_chan_avg_name = pc_chan_avg_som_cluster_name
elif pixel_cluster_col == 'pixel_meta_cluster_rename':
    pc_chan_avg_name = pc_chan_avg_meta_cluster_name


runtime_dict= {}
for size in dataset_sizes:
    print(size)

    fovs = random.sample(all_fovs, k=size)

    if not os.path.exists(os.path.join(base_dir, pixel_output_dir)):
        os.makedirs(os.path.join(base_dir, pixel_output_dir))

    if not os.path.exists(os.path.join(base_dir, "pixie", cell_output_dir)):
        os.mkdir(os.path.join(base_dir, "pixie", cell_output_dir))

    pixel_cluster_utils.create_pixel_matrix(
        fovs,
        channels,
        base_dir,
        tiff_dir,
        os.path.join(base_dir, segmentation_dir),
        img_sub_folder=img_sub_folder,
        seg_suffix=seg_suffix,
        pixel_output_dir=pixel_output_dir,
        data_dir=pixel_data_dir,
        subset_dir=pixel_subset_dir,
        norm_vals_name=norm_vals_name,
        blur_factor=blur_factor,
        subset_proportion=0.01,
        multiprocess=multiprocess,
        batch_size=batch_size
    )
    
    pixel_pysom = pixel_cluster_utils.train_pixel_som(
        fovs,
        channels,
        base_dir,
        subset_dir=pixel_subset_dir,
        norm_vals_name=norm_vals_name,
        som_weights_name=pixel_som_weights_name,
        num_passes=1
    )
    
    pixel_cluster_utils.cluster_pixels(
        fovs,
        channels,
        base_dir,
        pixel_pysom,
        data_dir=pixel_data_dir,
        multiprocess=multiprocess,
        batch_size=batch_size
    )
    
    pixel_cluster_utils.generate_som_avg_files(
        fovs,
        channels,
        base_dir,
        pixel_pysom,
        data_dir=pixel_data_dir,
        pc_chan_avg_som_cluster_name=pc_chan_avg_som_cluster_name
    )
    
    pixel_cc = pixel_cluster_utils.pixel_consensus_cluster(
        fovs,
        channels,
        base_dir,
        max_k=max_k,
        cap=cap,
        data_dir=pixel_data_dir,
        pc_chan_avg_som_cluster_name=pc_chan_avg_som_cluster_name,
        multiprocess=multiprocess,
        batch_size=batch_size
    )
    
    pixel_cluster_utils.generate_meta_avg_files(
        fovs,
        channels,
        base_dir,
        pixel_cc,
        data_dir=pixel_data_dir,
        pc_chan_avg_som_cluster_name=pc_chan_avg_som_cluster_name,
        pc_chan_avg_meta_cluster_name=pc_chan_avg_meta_cluster_name
    )

    for fov in fovs:
        one_fov_dat = feather.read_dataframe(os.path.join(base_dir, pixel_data_dir,fov+'.feather'))
        one_fov_dat = one_fov_dat.rename(columns={"pixel_meta_cluster":"pixel_meta_cluster_rename"})
        feather.write_dataframe(one_fov_dat, os.path.join(base_dir, pixel_data_dir,fov+'.feather'), compression='uncompressed')

    rename_dat = pd.read_csv(os.path.join(base_dir,pc_chan_avg_meta_cluster_name))
    rename_dat = rename_dat.rename(columns={"pixel_meta_cluster":"pixel_meta_cluster_rename"})
    rename_dat.to_csv(os.path.join(base_dir,pc_chan_avg_meta_cluster_name), index=False)

    start_time = time.time()

    cell_pysom = cell_cluster_utils.train_cell_som(
        fovs,
        channels,
        base_dir,
        pixel_data_dir=pixel_data_dir,
        cell_table_path=cell_table_path,
        cluster_counts_name=cluster_counts_name,
        cluster_counts_size_norm_name=cluster_counts_size_norm_name,
        pixel_cluster_col=pixel_cluster_col,
        pc_chan_avg_name=pc_chan_avg_name,
        som_weights_name=cell_som_weights_name,
        weighted_cell_channel_name=weighted_cell_channel_name,
        num_passes=10
    )
    
    cell_cluster_utils.cluster_cells(
        base_dir,
        cell_pysom,
        pixel_cluster_col_prefix=pixel_cluster_col,
        cell_som_cluster_count_avg_name=cell_som_cluster_count_avg_name
    )
    
    _ = cell_cluster_utils.cell_consensus_cluster(
        fovs=fovs,
        channels=channels,
        base_dir=base_dir,
        pixel_cluster_col=pixel_cluster_col,
        max_k=max_k,
        cap=cap,
        cluster_counts_size_norm_name=cluster_counts_size_norm_name,
        cell_som_cluster_count_avg_name=cell_som_cluster_count_avg_name,
        cell_meta_cluster_count_avg_name=cell_meta_cluster_count_avg_name,
        weighted_cell_channel_name=weighted_cell_channel_name,
        cell_som_cluster_channel_avg_name=cell_som_cluster_channel_avg_name,
        cell_meta_cluster_channel_avg_name=cell_meta_cluster_channel_avg_name
    )

    end_time = time.time()
    elapsed = end_time - start_time

    num_cell = len(cell_pysom.cell_data.index)
    print(num_cell)

    runtime_dict[num_cell] = elapsed
    shutil.rmtree(os.path.join(base_dir, pixel_output_dir))
    shutil.rmtree(os.path.join(base_dir, "pixie", cell_output_dir))


runtime_df = pd.DataFrame({'size':list(runtime_dict.keys()), 'runtime':list(runtime_dict.values())})
runtime_df.to_csv("pixie_runtimes_cell.csv", index=False)

