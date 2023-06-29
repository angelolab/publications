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

from ark.analysis import visualize
from ark.phenotyping import pixel_cluster_utils
from ark.utils import data_utils, io_utils, load_utils, plot_utils
from ark.utils.metacluster_remap_gui import colormap_helper, MetaClusterData, MetaClusterGui, metaclusterdata_from_files

random.seed(10)
dataset_sizes = [1,10,30,50]

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

runtime_dict= {}
for size in dataset_sizes:
    print(size)

    fovs = random.sample(all_fovs, k=size)

    start_time = time.time()

    if not os.path.exists(os.path.join(base_dir, pixel_output_dir)):
        os.makedirs(os.path.join(base_dir, pixel_output_dir))

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
        subset_proportion=1,
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

    end_time = time.time()
    elapsed = end_time - start_time

    num_px = len(pixel_pysom.train_data.index)
    print(num_px)

    runtime_dict[num_px] = elapsed
    shutil.rmtree(os.path.join(base_dir, pixel_output_dir))


runtime_df = pd.DataFrame({'size':list(runtime_dict.keys()), 'runtime':list(runtime_dict.values())})
runtime_df.to_csv("pixie_runtimes_pixels.csv", index=False)

