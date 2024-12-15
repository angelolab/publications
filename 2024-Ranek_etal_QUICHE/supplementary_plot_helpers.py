import re
import os
from os import PathLike
import json
import math
import natsort as ns
import numpy as np
import pandas as pd
import pathlib
import shutil
import datetime
import copy
import itertools
import warnings
from collections.abc import Iterable
from scipy.ndimage import gaussian_filter
from scipy import stats
from scipy.stats import ranksums
from statsmodels.stats.multitest import multipletests
import seaborn as sns
import skimage.io as io
from skimage.measure import label
from skimage import morphology
from skimage.io import imread
from skimage.segmentation import find_boundaries
import xarray as xr
from PIL import Image, ImageDraw, ImageFont
from typing import Dict, List, Optional, Tuple, TypedDict, Union, Literal, Iterable
import geopandas as gpd
from ark.utils import plot_utils, data_utils
from alpineer.io_utils import list_folders, list_files, remove_file_extensions, validate_paths
from alpineer.load_utils import load_imgs_from_tree, load_imgs_from_dir
from alpineer.misc_utils import verify_in_list
from alpineer import image_utils, io_utils, load_utils, misc_utils
from matplotlib.patches import Ellipse
from matplotlib.gridspec import GridSpec
from matplotlib import ticker
from matplotlib.axes import Axes
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter, MaxNLocator
import matplotlib.image as mpimg
import matplotlib.gridspec as gridspec
from matplotlib.colors import ListedColormap, Normalize
from dataclasses import dataclass, field
from pandas.core.groupby import DataFrameGroupBy
from requests.exceptions import HTTPError
from tqdm.auto import tqdm

# # NOTE: the following loop computes the mean for every marker across all images. precomputed to save time and memory :)
# BASE_DIR = "/Volumes/Shared/Noah Greenwald/TNBC_Cohorts/SPAIN"
# IMAGE_DIR = os.path.join(BASE_DIR, "image_data", "samples")

# exclude_chans = ["Au", "CD11c_nuc_exclude", "CK17_smoothed", "ECAD_smoothed", "FOXP3_nuc_include",
#                  "LAG3", "Noodle", "chan_39", "chan_45", "chan_48", "chan_115", "chan_141", "Calprotectin_old"]

# # plot two sample FOVs
# test_samples_fov = io_utils.list_folders(IMAGE_DIR)[0]
# channels = sorted(io_utils.remove_file_extensions(
#     io_utils.list_files(os.path.join(IMAGE_DIR, test_samples_fov), substrs=".tiff")
# ))
# for ec in exclude_chans:
#     if ec in channels:
#         channels.remove(ec)
# MEAN_PANEL_NORM = {}
# for chan in channels:
#     MEAN_PANEL_NORM[chan] = np.mean(load_imgs_from_tree(IMAGE_DIR, fovs=None, channels=[chan], max_image_size=2048))

# tiled regions
REGION_PARAM_FIELDS = [
    "region_name",
    "region_start_row",
    "region_start_col",
    "fov_num_row",
    "fov_num_col",
    "row_fov_size",
    "col_fov_size",
    "region_rand",
]

MICRON_TO_STAGE_X_MULTIPLIER = 0.001001
MICRON_TO_STAGE_X_OFFSET = 0.3116
MICRON_TO_STAGE_Y_MULTIPLIER = 0.001018
MICRON_TO_STAGE_Y_OFFSET = 0.6294

STAGE_LEFT_BOUNDARY = -1.03
STAGE_RIGHT_BOUNDARY = 23.19
STAGE_TOP_BOUNDARY = 56.16
STAGE_BOTTOM_BOUNDARY = 1.57
OPTICAL_LEFT_BOUNDARY = 386.0
OPTICAL_RIGHT_BOUNDARY = 730.1
OPTICAL_TOP_BOUNDARY = 302.1
OPTICAL_BOTTOM_BOUNDARY = 1085.7

OPTICAL_BOUNDARY_TOL = 100
STAGE_BOUNDARY_TOL = 1

COREG_PARAM_BASELINE = {
    "STAGE_TO_OPTICAL_X_MULTIPLIER": 15.05,
    "STAGE_TO_OPTICAL_X_OFFSET": 26.15,
    "STAGE_TO_OPTICAL_Y_MULTIPLIER": -15.03,
    "STAGE_TO_OPTICAL_Y_OFFSET": -76.16,
}


# QC channels to ignore
QC_CHANNEL_IGNORE = ["Au", "Fe", "Na", "Ta", "Noodle"]

# QC metric .csv suffix and column naming
QC_SUFFIXES = ["nonzero_mean_stats", "total_intensity_stats", "percentile_99_9_stats"]
QC_COLUMNS = ["Non-zero mean intensity", "Total intensity", "99.9% intensity value"]

MEAN_PANEL_NORM = {"CD11c": 0.006193546577196089, "CD14": 0.01984555177226498, "CD163": 0.026490620436259955,
    "CD20": 0.012355682807796918, "CD3": 0.006915745669193154, "CD31": 0.018580328706651567,
    "CD38": 0.014254272705785212, "CD4": 0.011660068838085572, "CD45": 0.015016060634967094,
    "CD45RB": 0.008236598627789901, "CD45RO": 0.01470480803636466, "CD56": 0.0039886591958356934,
    "CD57": 0.012048721429121926, "CD68": 0.011606977635707979, "CD69": 0.008835169089640722,
    "CD8": 0.01140980883839861, "CK17": 0.015449040598057523, "Calprotectin": 0.00495033742848854,
    "ChyTr": 0.027970794698707765, "Collagen1": 0.022180374726308422, "ECAD": 0.02324031755306159,
    "FAP": 0.021780513481618562, "FOXP3": 0.00494151681211686, "Fe": 0.34932304124394165,
    "Fibronectin": 0.02718638734057556, "GLUT1": 0.019362882625847296,
    "H3K27me3": 0.07062930678187326, "H3K9ac": 0.07087346982563525, "HLA1": 0.022028920388760115,
    "HLADR": 0.014832535896920995, "IDO": 0.00431968466707603, "Ki67": 0.030366892417654723,
    "PD1": 0.003349747752931683, "PDL1": 0.007616826308262865, "SMA": 0.2710457265857868,
    "TBET": 0.008260657932221848, "TCF1": 0.006155141651624279, "TIM3": 0.006329398943399673,
    "Vim": 0.06671803387741954}

MARKER_INFO = {
    "Ki67": {
        "populations": ["Cancer", "Mast"],
        "threshold": 0.002,
        "x_range": (0, 0.012),
        "x_ticks": np.array([0, 0.004, 0.008, 0.012]),
        "x_tick_labels": np.array([0, 0.004, 0.008, 0.012]),
    },
    "CD38": {
        "populations": ["Endothelium", "Cancer_EMT"],
        "threshold": 0.004,
        "x_range": (0, 0.02),
        "x_ticks": np.array([0, 0.005, 0.01, 0.015, 0.02]),
        "x_tick_labels": np.array([0, 0.005, 0.01, 0.015, 0.02]),
    },
    "CD45RB": {
        "populations": ["CD4T", "Stroma"],
        "threshold": 0.001,
        "x_range": (0, 0.015),
        "x_ticks": np.array([0, 0.005, 0.010, 0.015]),
        "x_tick_labels": np.array([0, 0.005, 0.010, 0.015])
    },
    "CD45RO": {
        "populations": ["CD4T", "Fibroblast"],
        "threshold": 0.002,
        "x_range": (0, 0.02),
        "x_ticks": np.array([0, 0.005, 0.01, 0.015, 0.02]),
        "x_tick_labels": np.array([0, 0.005, 0.01, 0.015, 0.02])
    },
    "CD57": {
        "populations": ["CD8T", "B"],
        "threshold": 0.002,
        "x_range": (0, 0.006),
        "x_ticks": np.array([0, 0.002, 0.004, 0.006]),
        "x_tick_labels": np.array([0, 0.002, 0.004, 0.006])
    },
    "CD69": {
        "populations": ["Treg", "Cancer"],
        "threshold": 0.002,
        "x_range": (0, 0.008),
        "x_ticks": np.array([0, 0.002, 0.004, 0.006, 0.008]),
        "x_tick_labels": np.array([0, 0.002, 0.004, 0.006, 0.008])
    },
    "GLUT1": {
        "populations": ["Cancer_EMT", "M2_Mac"],
        "threshold": 0.002,
        "x_range": (0, 0.02),
        "x_ticks": np.array([0, 0.005, 0.01, 0.015, 0.02]),
        "x_tick_labels": np.array([0, 0.005, 0.01, 0.015, 0.02])
    },
    "IDO": {
        "populations": ["APC", "M1_Mac"],
        "threshold": 0.001,
        "x_range": (0, 0.003),
        "x_ticks": np.array([0, 0.001, 0.002, 0.003]),
        "x_tick_labels": np.array([0, 0.001, 0.002, 0.003])
    },
    "PD1": {
        "populations": ["CD8T", "Stroma"],
        "threshold": 0.0005,
        "x_range": (0, 0.002),
        "x_ticks": np.array([0, 0.0005, 0.001, 0.0015, 0.002]),
        "x_tick_labels": np.array([0, 0.0005, 0.001, 0.0015, 0.002])
    },
    "PDL1": {
        "populations": ["Cancer", "Stroma"],
        "threshold": 0.001,
        "x_range": (0, 0.003),
        "x_ticks": np.array([0, 0.001, 0.002, 0.003]),
        "x_tick_labels": np.array([0, 0.001, 0.002, 0.003]),
    },
    "HLA1": {
        "populations": ["APC", "Stroma"],
        "threshold": 0.001,
        "x_range": (0, 0.025),
        "x_ticks": np.array([0, 0.0125, 0.025]),
        "x_tick_labels": np.array([0, 0.0125, 0.025])
    },
    "HLADR": {
        "populations": ["APC", "Neutrophil"],
        "threshold": 0.001,
        "x_range": (0, 0.025),
        "x_ticks": np.array([0, 0.0125, 0.025]),
        "x_tick_labels": np.array([0, 0.0125, 0.025])
    },
    "TBET": {
        "populations": ["NK", "B"],
        "threshold": 0.0015,
        "x_range": (0, 0.0045),
        "x_ticks": np.array([0, 0.0015, 0.003, 0.0045]),
        "x_tick_labels": np.array([0, 0.0015, 0.003, 0.0045])
    },
    "TCF1": {
        "populations": ["CD4T", "M1_Mac"],
        "threshold": 0.001,
        "x_range": (0, 0.003),
        "x_ticks": np.array([0, 0.001, 0.002, 0.003]),
        "x_tick_labels": np.array([0, 0.001, 0.002, 0.003])
    },
    "TIM3": {
        "populations": ["Monocyte", "Endothelium"],
        "threshold": 0.001,
        "x_range": (0, 0.004),
        "x_ticks": np.array([0, 0.001, 0.002, 0.003, 0.004]),
        "x_tick_labels": np.array([0, 0.001, 0.002, 0.003, 0.004])
    },
    "Vim": {
        "populations": ["Endothelium", "B"],
        "threshold": 0.002,
        "x_range": (0, 0.06),
        "x_ticks": np.array([0, 0.02, 0.04, 0.06]),
        "x_tick_labels": np.array([0, 0.02, 0.04, 0.06])
    },
    "Fe": {
        "populations": ["Fibroblast", "Cancer"],
        "threshold": 0.1,
        "x_range": (0, 0.3),
        "x_ticks": np.array([0, 0.1, 0.2, 0.3]),
        "x_tick_labels": np.array([0, 0.1, 0.2, 0.3]),
    }
}

def validate_panel(
    data_dir: Union[str, pathlib.Path], fov: str, save_dir: Union[str, pathlib.Path], 
    channels: Optional[List[str]] = None, img_sub_folder: str = "", padding: int = 10,
    num_rows: Optional[int] = None, font_size: int = 200
):
    """Given a FOV in an image folder, stitch and annotate each channel

    Args:
        data_dir (Union[str, pathlib.Path]):
            The directory containing the image data for each FOV
        fov (str):
            The name of the FOV to stitch
        save_dir (Union[str, pathlib.Path]):
            The directory to save the stitched image
        channels (Optional[List[str]]):
            The list of channels to tile. If None, uses all.
        img_sub_folder (str):
            The sub folder name inside each FOV directory, set to "" if None
        padding (int):
            Amount of padding to add around each channel in the stitched image
        num_rows (int):
            The number of rows, if None uses the rounded sqrt of total num of images
        font_size (int):
            The font size to use for annotations
    """
    # verify the FOV is valid
    all_fovs: List[str] = list_folders(data_dir)
    verify_in_list(
        specified_fov=fov,
        valid_fovs=all_fovs
    )

    # verify save_dir is valid before defining the save path
    validate_paths([save_dir])
    stitched_img_path: pathlib.Path = pathlib.Path(save_dir) / f"{fov}_channels_stitched.tiff"

    # validate the channels provided, or set to all if None
    all_channels = remove_file_extensions(list_files(os.path.join(data_dir, fov), substrs=".tiff"))
    if not channels:
        channels = all_channels
    verify_in_list(
        specified_channels=channels,
        valid_channels=all_channels
    )

    # sort the channels to ensure they get tiled in alphabetical order, regardless of case
    channels = sorted(channels, key=str.lower)

    # load the data and get the channel names and image dimensions
    image_data: xr.DataArray = load_imgs_from_tree(
        data_dir=data_dir, fovs=[fov], channels=channels, img_sub_folder=img_sub_folder
    )[0, ...]

    # normalize each channel by the mean 99.9% value across the cohort, for clearer visualization
    for chan in image_data.channels.values:
        norm_val = MEAN_PANEL_NORM[chan]
        image_data.loc[..., chan] = image_data.loc[..., chan] / norm_val

    # ensure the channels dimension is the 0th for annotation purposes
    image_data = image_data.transpose("channels", "rows", "cols")

    # generate the stitched image and save
    panel_tiled: Image = stitch_and_annotate_padded_img(
        image_data, padding=padding, num_rows=num_rows, font_size=font_size, annotate=False
    )

    panel_tiled.save(stitched_img_path)

def stitch_and_annotate_padded_img(image_data: xr.DataArray, padding: int = 25,
                                   num_rows: Optional[int] = None, font_size: int = 100,
                                   annotate: bool = False, step: int = 1):
    """Stitch an image with (c, x, y) dimensions. If specified, annotate each image with labels
    contained in the cth dimension.

    Args:
        image_data (xr.DataArray):
            The image data to tile, should be 3D
        padding (int):
            Amount of padding to add around each channel in the stitched image
        num_rows (int):
            The number of rows, if None uses the rounded sqrt of total num of images
        font_size (int):
            The font size to use for annotations
        annotate (bool):
            Whether to annotate the images with labels in dimension c
        step (int):
            The step size to use before adding an image to the tile

    Returns:
        Image:
            The PIL image instance which contains the stitched (with padding) and annotated image
    """
    # param validation
    if padding < 0:
        raise ValueError("padding must be a non-negative integer")
    if font_size <= 0:
        raise ValueError("font_size must be a positive integer")

    images_to_select = np.arange(0, image_data.shape[0], step=step)
    image_data = image_data[images_to_select, ...]

    # define the number of rows and columns
    if num_rows:
        num_cols = math.ceil(image_data.shape[0] / num_rows)
    else:
        num_cols: int = math.isqrt(image_data.shape[0])
        num_rows: int = math.ceil(image_data.shape[0] / num_cols)
    row_len: int = image_data.shape[1]
    col_len: int = image_data.shape[2]
    # create the blank image, start with a fully white slate
    stitched_image: np.ndarray = np.zeros(
        (
            num_rows * row_len + (num_rows - 1) * padding,
            num_cols * col_len + (num_cols - 1) * padding
        )
    )
    stitched_image.fill(255)

    # retrieve the annotation labels
    annotation_labels = list(image_data.coords[image_data.dims[0]].values)

    # stitch the channels
    img_idx: int = 0
    for row in range(num_rows):
        for col in range(num_cols):
            stitched_image[
                (row * row_len + padding * row) : ((row + 1) * row_len + padding * row),
                (col * col_len + padding * col) : ((col + 1) * col_len + padding * col)
            ] = image_data[img_idx, ...]
            img_idx += 1
            if img_idx == len(annotation_labels):
                break

    # define a draw instance for annotating the channel name
    stitched_image_im: Image = Image.fromarray(stitched_image)

    # annotate with labels in c-axis if arg set
    if annotate:
        imdraw: ImageDraw = ImageDraw.Draw(stitched_image_im)
        imfont: ImageFont = ImageFont.truetype("Arial Unicode.ttf", font_size)

        img_idx = 0
        for row in range(num_rows):
            for col in range(num_cols):
                imdraw.text(
                    (col * col_len + padding * col, row * row_len + padding * row),
                    annotation_labels[img_idx],
                    font=imfont,
                    fill=255
                )
                img_idx += 1
                if img_idx == len(annotation_labels):
                    break

    return stitched_image_im

class QuantileNormalization(Normalize):
    def __init__(self,
                 vmin: float = None,
                 vmax: float = None,
                 q: tuple[float, float] = (0.01, 0.99),
                 clip: bool = False,
                 eps: float = 1e-20,
    ) -> None:
        """Normalizes the input data by the qth quantile.

        Args
        ----------
        vmin : float, optional
            If vmin is not given it is initilaized from the minimum value of the
            array, by default None
        vmax : float, optional
            If vmax is not given it is initilaized from the maximum value of the
            array, by default None
        q : tuple[float, float], optional
            A tuple of quatiles where the smallest element is the minimum quantile
            and the largest element is the maximum percentile, by default (0.01, 0.99). Must
            be between 0 and 1 inclusive.
        clip: bool, optional
            If True, the normalized values are clipped to the range [0, 1], by default False
        eps: float, optional
            Small value to add to the denominator to avoid division by zero, by default 1e-20
        """
        super().__init__(vmin, vmax)
        if isinstance(q, tuple):
            if len(q) != 2:
                raise ValueError("q must be a tuple of length 2")
            if not all(0 <= i <= 1 for i in q):
                raise ValueError("q's elements must be between 0 and 1 inclusive")
        else:
            raise ValueError("q must be a tuple")
        
        self.qmin = min(q)
        self.qmax = max(q)
        self.clip = clip
        self.eps = eps

    def __call__(self, value):
        val_qmin, val_qmax = np.quantile(value, [self.qmin, self.qmax])
        
        norm = (value - val_qmin) / (val_qmax - val_qmin + self.eps)
        if self.clip:
            norm = np.clip(norm, 0, 1)
            
        return norm
    
def stitch_before_after_rosetta(
    pre_rosetta_dir: Union[str, pathlib.Path], post_rosetta_dir: Union[str, pathlib.Path],
    save_dir: Union[str, pathlib.Path],
    run_name: str, fov_indices: Optional[List[int]], target_channel: str, source_channel: str = "Noodle",
    pre_rosetta_subdir: str = "", post_rosetta_subdir: str = "",
    img_size_scale: float = 0.5, percent_norm: Optional[float] = 99.999,
    padding: int = 25, font_size: int = 175, step: int = 1,
    save_separate: bool = False
):
    """Generates two stitched images: before and after Rosetta

    pre_rosetta_dir (Union[str, pathlib.Path]):
        The directory containing the run data before Rosetta
    post_rosetta_dir (Union[str, pathlib.Path]):
        The directory containing the run data after Rosetta
    save_dir (Union[str, pathlib.Path]):
        The directory to save both the pre- and post-Rosetta tiled images
    run_name (str):
        The name of the run to tile, should be present in both pre_rosetta_dir and post_rosetta_dir
    fov_indices (Optional[List[int]]):
        The list of indices to select. If None, use all.
    target_channel (str):
        The name of the channel to tile inside run_name
    source_channel (str):
        The name of the source channel that was subtracted from target_channel
    pre_rosetta_subdir (str):
        If applicable, the name of the subdirectory inside each FOV folder of the pre-Rosetta data
    post_rosetta_subdir (str):
        If applicable, the name of the subdirectory inside each FOV folder of the post-Rosetta data
    percent_norm (int):
        Percentile normalization param to enable easy visualization, if None then skips this step
    img_size_scale (float):
        Amount to scale down image. Set to None for no scaling
    padding (int):
        Amount of padding to add around each channel in the stitched image
    font_size (int):
        The font size to use for annotations
    step (int):
        The step size to use before adding an image to the tile
    save_separate (bool):
        If set, then save each FOV separately, otherwise save full tile
    """
    # verify that the run_name specified appears in both pre and post norm folders
    all_pre_rosetta_runs: List[str] = list_folders(pre_rosetta_dir)
    all_post_rosetta_runs: List[str] = list_folders(post_rosetta_dir)
    verify_in_list(
        specified_run=run_name,
        all_pre_norm_runs=all_pre_rosetta_runs
    )
    verify_in_list(
        specified_run=run_name,
        all_post_norm_runs=all_post_rosetta_runs
    )

    # verify save_dir is valid before defining the save paths
    validate_paths([save_dir])
    rosetta_stitched_path: pathlib.Path = \
        pathlib.Path(save_dir) / f"{run_name}_{target_channel}_pre_post_Rosetta.tiff"

    # define full paths to pre- and post-Rosetta data
    pre_rosetta_run_path: pathlib.Path = pathlib.Path(pre_rosetta_dir) / run_name
    post_rosetta_run_path: pathlib.Path = pathlib.Path(post_rosetta_dir) / run_name

    # get all the FOVs in natsorted order
    # NOTE: assumed that the FOVs are the same pre and post, since the run names are the same
    all_fovs: List[str] = ns.natsorted(list_folders(pre_rosetta_run_path))
    all_fovs = [item for item in all_fovs if 'stitched' not in item]
    # load Noodle, pre-, and post-Rosetta data in acquisition order, drop channel axis as it's 1-D
    # ensure the pre-Rosetta Noodle is loaded
    noodle_data: xr.DataArray = load_imgs_from_tree(
        data_dir=pre_rosetta_run_path, fovs=all_fovs, channels=[source_channel],
        img_sub_folder=pre_rosetta_subdir, max_image_size=2048
    )[..., 0]
    pre_rosetta_data: xr.DataArray = load_imgs_from_tree(
        data_dir=pre_rosetta_run_path, fovs=all_fovs, channels=[target_channel],
        img_sub_folder=pre_rosetta_subdir, max_image_size=2048
    )[..., 0]
    post_rosetta_data: xr.DataArray = load_imgs_from_tree(
        data_dir=post_rosetta_run_path, fovs=all_fovs, channels=[target_channel],
        img_sub_folder=post_rosetta_subdir, max_image_size=2048
    )[..., 0]

    # divide pre-Rosetta by 200 to ensure same scale
    pre_rosetta_data = pre_rosetta_data / 200

    # reassign coordinate with FOV names that don't contain "-scan-1" or additional dashes
    fovs_condensed: np.ndarray = np.array([f"FOV{af.split('-')[1]}" for af in all_fovs])
    noodle_data = noodle_data.assign_coords({"fovs": fovs_condensed})
    pre_rosetta_data = pre_rosetta_data.assign_coords({"fovs": fovs_condensed})
    post_rosetta_data = post_rosetta_data.assign_coords({"fovs": fovs_condensed})

    # the top should be original, middle Noodle, bottom Rosetta-ed
    # NOTE: leave out Noodle row from dimensions for now for proper rescaling and percentile norm
    stitched_pre_post_rosetta: np.ndarray = np.zeros(
        (2048 * 2, 2048 * len(fovs_condensed))
    )
    for fov_i, fov_name in enumerate(fovs_condensed):
        # add the rescaled pre- and post-Rosetta images first
        stitched_pre_post_rosetta[
            0:2048, (2048 * fov_i):(2048 * (fov_i + 1))
        ] = pre_rosetta_data[fov_i, ...].values
        stitched_pre_post_rosetta[
            2048:4096, (2048 * fov_i):(2048 * (fov_i + 1))
        ] = post_rosetta_data[fov_i, ...].values

    # define the Noodle row
    stitched_noodle: np.ndarray = np.zeros((2048, 2048 * len(fovs_condensed)))
    for fov_i, fov_name in enumerate(fovs_condensed):
        stitched_noodle[
            :, (2048 * fov_i):(2048 * (fov_i + 1))
        ] = noodle_data[fov_i, ...].values

    # run percent normalization on Noodle data if specified
    if percent_norm:
        source_percentile: float = np.percentile(stitched_noodle, percent_norm)
        non_source_percentile: float = np.percentile(stitched_pre_post_rosetta, percent_norm)
        perc_ratio: float = source_percentile / non_source_percentile
        stitched_noodle = stitched_noodle / perc_ratio

    # combine the Noodle data with the stitched data, swap so that Noodle is in the middle
    stitched_pre_post_rosetta = np.vstack(
        (stitched_pre_post_rosetta[:2048, :], stitched_noodle, stitched_pre_post_rosetta[2048:, :])
    )

    # subset on just the FOV indices selected
    # NOTE: because of how percent norm works, better to run on all FOVs first to ensure brightness
    # as opposed to subsetting first, which often leads to dimmer images
    if fov_indices:
        indices_select = []
        for fi in fov_indices:
            indices_select.extend(list(np.arange(2048 * fi, 2048 * (fi + 1))))
        stitched_pre_post_rosetta = stitched_pre_post_rosetta[:, indices_select]

    if save_separate:
        # save each individual image separately
        for fov_i, fov_num in enumerate(fov_indices):
            stitched_rosetta_pil: Image = Image.fromarray(
                stitched_pre_post_rosetta[:, (2048 * fov_i):(2048 * (fov_i + 1))]
            )
            rosetta_stitched_path: pathlib.Path = \
                pathlib.Path(save_dir) / f"{target_channel}_image_{fov_num}.tiff"
            stitched_rosetta_pil.save(rosetta_stitched_path)

    else:
        # save the full stitched image
        stitched_rosetta_pil: Image = Image.fromarray(np.round(stitched_pre_post_rosetta, 3))
        stitched_rosetta_pil.save(rosetta_stitched_path)

def stitch_before_after_norm(
    pre_norm_dir: Union[str, pathlib.Path], post_norm_dir: Union[str, pathlib.Path],
    save_dir: Union[str, pathlib.Path], run_name: str,
    fov_indices: Optional[List[int]], channel: str, pre_norm_subdir: str = "",
    post_norm_subdir: str = "", padding: int = 25, font_size: int = 100, step: int = 1, num_rows = 7):
    """Generates two stitched images: before and after normalization

    pre_norm_dir (Union[str, pathlib.Path]):
        The directory containing the run data before normalization
    post_norm_dir (Union[str, pathlib.Path]):
        The directory containing the run data after normalization
    save_dir (Union[str, pathlib.Path]):
        The directory to save both the pre- and post-norm tiled images
    run_name (str):
        The name of the run to tile, should be present in both pre_norm_dir and post_norm_dir
    fov_indices (Optional[List[int]]):
        The list of indices to select. If None, use all.
    channel (str):
        The name of the channel to tile inside run_name
    pre_norm_subdir (str):
        If applicable, the name of the subdirectory inside each FOV folder of the pre-norm data
    post_norm_subdir (str):
        If applicable, the name of the subdirectory inside each FOV folder of the post-norm data
    padding (int):
        Amount of padding to add around each channel in the stitched image
    font_size (int):
        The font size to use for annotations
    step (int):
        The step size to use before adding an image to the tile
    """
    # verify that the run_name specified appears in both pre and post norm folders
    all_pre_norm_runs: List[str] = list_folders(pre_norm_dir)
    all_post_norm_runs: List[str] = list_folders(post_norm_dir)
    verify_in_list(
        specified_run=run_name,
        all_pre_norm_runs=all_pre_norm_runs
    )
    verify_in_list(
        specified_run=run_name,
        all_post_norm_runs=all_post_norm_runs
    )

    # verify save_dir is valid before defining the save paths
    validate_paths([save_dir])
    pre_norm_stitched_path: pathlib.Path = \
        pathlib.Path(save_dir) / f"{run_name}_{channel}_pre_norm_stitched.tiff"
    post_norm_stitched_path: pathlib.Path = \
        pathlib.Path(save_dir) / f"{run_name}_{channel}_post_norm_stitched.tiff"

    pre_norm_run_path: pathlib.Path = pathlib.Path(pre_norm_dir) / run_name
    post_norm_run_path: pathlib.Path = pathlib.Path(post_norm_dir) / run_name

    # get all the FOVs in natsorted order
    # NOTE: assumed that the FOVs are the same pre and post, since the run names are the same
    all_fovs: List[str] = ns.natsorted(list_folders(pre_norm_run_path))

    # load pre- and post-norm data in acquisition order, drop channel axis as it's 1-D
    pre_norm_data: xr.DataArray = load_imgs_from_tree(
        data_dir=pre_norm_run_path, fovs=all_fovs, channels=[channel],
        img_sub_folder=pre_norm_subdir, max_image_size=2048
    )[..., 0]
    post_norm_data: xr.DataArray = load_imgs_from_tree(
        data_dir=post_norm_run_path, fovs=all_fovs, channels=[channel],
        img_sub_folder=post_norm_subdir, max_image_size=2048
    )[..., 0]

    if fov_indices:
        pre_norm_data = pre_norm_data[fov_indices, ...]
        post_norm_data = post_norm_data[fov_indices, ...]

    # reassign coordinate with FOV names that don't contain "-scan-1" or additional dashes
    fovs_condensed: np.ndarray = np.array([f"FOV{af.split('-')[1]}" for af in all_fovs])
    if fov_indices:
        fovs_condensed = fovs_condensed[fov_indices]

    pre_norm_data = pre_norm_data.assign_coords({"fovs": fovs_condensed})
    post_norm_data = post_norm_data.assign_coords({"fovs": fovs_condensed})

    # generate and save the pre- and post-norm tiled images
    pre_norm_tiled: Image = stitch_and_annotate_padded_img(pre_norm_data, num_rows, padding, font_size, step=step)
    post_norm_tiled: Image = stitch_and_annotate_padded_img(post_norm_data, num_rows, padding, font_size, step=step)

    pre_norm_tiled.save(pre_norm_stitched_path)
    post_norm_tiled.save(post_norm_stitched_path)

class MembraneMarkersSegmentationPlot:
    def __init__(
        self,
        fov: str,
        image_data: PathLike,
        segmentation_dir: PathLike,
        membrane_channels,
        overlay_channels,
        q: tuple[float, float] = (0.05, 0.95),
        clip: bool = False,
        figsize: Tuple[int, int] = (12, 4),
        layout: Literal["constrained", "tight"] = None,
        image_type: Literal["png", "pdf", "svg"] = "pdf",
    ):
        """Creates a figure with two subplots, one for each membrane marker used for segmentation,
        and one for the overlay of the membrane and nuclear markers.

        Args
        ----------
        fov : str
            The name of the FOV to be plotted
        image_data : PathLike
            The directory containing the image data.
        segmentation_dir : PathLike
            The directory containing the segmentation data.
        membrane_channels : List[str]
            The names of the membrane markers to be plotted.
        overlay_channels : str | List[str]
            The overlay channels to be plotted, can be either "nuclear_channel",
            "membrane_channel", or both.
        overlay_cmap: str, optional
            The colormap to use for the overlay, by default "viridis_r"
        q : tuple[float, float], optional
            A tuple of quatiles where the smallest element is the minimum quantile
            and the largest element is the maximum percentile, by default (0.05, 0.95). Must
            be between 0 and 1 inclusive.
        clip : bool, optional
            If True, the normalized values are clipped to the range [0, 1], by default False
        figsize : Tuple[int, int], optional
            The size of the figure, by default (8, 4)
        layout : Literal["constrained", "tight"], optional
            The layout engine, defaults to None, by default None
        image_type : Literal["png", "pdf", "svg"], optional
            The file type to save the plot as, by default "pdf"
        """
        self.fov_name = fov
        self.membrane_channels = membrane_channels
        self.overlay_channels = overlay_channels
        self.figsize = figsize
        self.layout = layout
        self.n_chans = len(set(membrane_channels))
        self.q = q
        self.image_data = image_data
        self.seg_dir = segmentation_dir
        self.image_type = image_type
        self.clip = clip

        self.fig = plt.figure(figsize=figsize, layout=layout)
        self.subfigs = self.fig.subfigures(
            nrows=1, ncols=2, wspace=0.05, width_ratios=[1, 1]
        )

    def make_plot(self, save_dir: PathLike):
        """Plots the membrane markers and overlay and saves the figure to the specified directory.

        Args
        ----------
        save_dir : PathLike
            The directory to save the figure to.
        """
        self.fov_xr = load_imgs_from_tree(
            data_dir=self.image_data,
            fovs=[self.fov_name],
            channels=self.membrane_channels,
        )

        self.fov_overlay = plot_utils.create_overlay(
            fov=self.fov_name,
            segmentation_dir=self.seg_dir +  "/deepcell_output",
            data_dir=self.seg_dir +  "/deepcell_input",
            img_overlay_chans=self.overlay_channels,
            seg_overlay_comp="whole_cell",
        )

        self.fov_cell_segmentation = load_imgs_from_dir(
            data_dir=self.seg_dir +  "/deepcell_output",
            files=[f"{self.fov_name}_whole_cell.tiff"],
        )

        self.fig.suptitle(
            t=f"{self.fov_name} Membrane Markers and Segmentation", fontsize=8
        )

        self._plot_mem_markers()
        self._plot_overlay_segmentation()
        self.fig.savefig(
            fname=save_dir / f"{self.fov_name}_membrane_markers_overlay.{self.image_type}",
        )
        plt.close(self.fig)

    def _plot_mem_markers(self):
        self.subfigs[0].suptitle("Membrane Markers", fontsize=6)

        markers_subplots = self.subfigs[0].subplots(
            nrows=2,
            ncols=int(np.ceil((self.n_chans + 1) / 2)),
            sharex=True,
            sharey=True,
            gridspec_kw={"wspace": 0.05, "hspace": 0.05},
        )

        channel_axes = markers_subplots.flat[: self.n_chans]

        self.subfigs[0].add_subplot

        for ax, channel in zip(channel_axes, ns.natsorted(self.membrane_channels)):
            chan_data = self.fov_xr.sel({"channels": channel}).squeeze()

            ax.imshow(
                X=chan_data,
                cmap="gray",
                norm=QuantileNormalization(q=self.q, clip=self.clip),
                interpolation="none",
                aspect="equal",
            )
            ax.set_title(channel, fontsize=6)
            remove_ticks(ax, "xy")

        ax_sum = markers_subplots.flat[self.n_chans]

        ax_sum.imshow(
            X=self.fov_xr.sum("channels").squeeze(),
            cmap="gray",
            norm=QuantileNormalization(q=self.q, clip=self.clip),
            interpolation="none",
            aspect="equal",
        )

        ax_sum.set_title("Sum", fontsize=6)
        remove_ticks(ax_sum, "xy")

        # Clean up and remove the empty subplots
        for ax in markers_subplots.flat[self.n_chans + 1 :]:
            ax.remove()

    def _plot_overlay_segmentation(self)-> None:
        cell_seg_ax, overlay_ax = self.subfigs[1].subplots(
            nrows=1, ncols=2, sharex=True, sharey=True
        )
        overlay_ax.set_title("Nuclear and Membrane Overlay", fontsize=6)
        overlay_ax.imshow(
            X=self.fov_overlay,
            norm=QuantileNormalization(q=self.q, clip=self.clip),
            interpolation="none",
            aspect="equal",
        )

        cell_seg_ax.set_title("Cell Segmentation", fontsize=6)
        cell_seg_ax.imshow(
            X=xr.apply_ufunc(
                mask_erosion_ufunc,
                self.fov_cell_segmentation.squeeze(),
                input_core_dims=[["rows", "cols"]],
                output_core_dims=[["rows", "cols"]],
                kwargs={"connectivity": 1, "mode": "thick"},
            ).pipe(lambda x: x.where(cond=x < 1, other=0.5)),
            cmap="grey",
            interpolation="none",
            aspect="equal",
            vmin=0,
            vmax=1,
        )

        remove_ticks(overlay_ax, "xy")
        remove_ticks(cell_seg_ax, "xy")

class SegmentationOverlayPlot:
    def __init__(
        self,
        fov: str,
        segmentation_dir: PathLike,
        overlay_channels = ["nuclear_channel", "membrane_channel"],
        q: tuple[float, float] = (0.05, 0.95),
        clip: bool = False,
        figsize: tuple[float, float] =(8, 4),
        layout: Literal["constrained", "tight"] = "constrained",
        image_type: Literal["png", "pdf", "svg"] = "svg",
    ) -> None:
        """Creates a figure with two subplots, one for the cell segmentation and one for the overlay

        Parameters
        ----------
        fov : str
            The name of the FOV to be plotted
        segmentation_dir : PathLike
            The directory containing the segmentation data.
        overlay_channels : str | List[str]
            The overlay channels to be plotted, can be either/both "nuclear_channel" or "membrane_channel",
            defaults to ["nuclear_channel", "membrane_channel"].
        q : tuple[float, float], optional
            A tuple of quatiles where the smallest element is the minimum quantile
            and the largest element is the maximum percentile, by default (0.05, 0.95). Must
            be between 0 and 1 inclusive.
        clip : bool, optional
            If True, the normalized values are clipped to the range [0, 1], by default False
        figsize : Tuple[int, int], optional
            The size of the figure, by default (8, 4)
        layout : Literal["constrained", "tight"], optional
            The layout engine, defaults to None, by default None
        image_type : Literal["png", "pdf", "svg"], optional
            The file type to save the plot as, by default "pdf"
        """
        self.fov_name = fov
        self.seg_dir = segmentation_dir
        self.overlay_channels = overlay_channels
        self.q = q
        self.clip = clip
        self.figsize = figsize
        self.layout = layout
        self.image_type = image_type

        self.fig = plt.figure(figsize=self.figsize, layout=layout)

    def make_plot(self, save_dir: PathLike) -> None:
        self.fov_overlay = plot_utils.create_overlay(
            fov=self.fov_name,
            segmentation_dir=self.seg_dir + "/deepcell_output",
            data_dir=self.seg_dir + "/deepcell_input",
            img_overlay_chans=self.overlay_channels,
            seg_overlay_comp="whole_cell",
        )
        self.fov_cell_segmentation = load_imgs_from_dir(
            data_dir=self.seg_dir + "/deepcell_output",
            files=[f"{self.fov_name}_whole_cell.tiff"],
        )
        self.fig.suptitle(
            t=f"{self.fov_name} Cell Segmentation and Overlay", fontsize=8
        )
        self._plot_overlay_segmentation()
        self.fig.savefig(
            save_dir / f"{self.fov_name}_segmentation_overlay.{self.image_type}"
        )

        plt.close(self.fig)

    def _plot_overlay_segmentation(self):
        cell_seg_ax, overlay_ax = self.fig.subplots(
            nrows=1, ncols=2, sharex=True, sharey=True
        )
        overlay_ax.set_title("Nuclear and Membrane Overlay", fontsize=6)
        overlay_ax.imshow(
            X=self.fov_overlay,
            norm=QuantileNormalization(q=self.q, clip=self.clip),
            interpolation="none",
            aspect="equal",
        )

        cell_seg_ax.set_title("Cell Segmentation", fontsize=6)
        cell_seg_ax.imshow(
            X=xr.apply_ufunc(
                mask_erosion_ufunc,
                self.fov_cell_segmentation.squeeze(),
                input_core_dims=[["rows", "cols"]],
                output_core_dims=[["rows", "cols"]],
                kwargs={"connectivity": 2, "mode": "thick"},
            ).pipe(lambda x: x.where(cond=x < 1, other=0.5)),
            cmap="grey",
            interpolation="none",
            aspect="equal",
            vmin=0,
            vmax=1,
        )

        remove_ticks(overlay_ax, "xy")
        remove_ticks(cell_seg_ax, "xy")

class CorePlot:
    def __init__(
        self,
        fov: str,
        hne_path: PathLike,
        seg_dir: PathLike,
        overlay_channels: list[str] = ["nuclear_channel", "membrane_channel"],
        figsize: tuple[float, float] = (13, 4),
        layout: Literal["constrained", "tight"] = "constrained",
        image_type: Literal["png", "pdf", "svg"] = "pdf",
    ):
        """Generates a figure with three subplots: one for the HnE core, one for the HnE FOV crop,
        and one for the overlay of the nuclear and membrane channels.

        Parameters
        ----------
        fov : str
            The name of the FOV to be plotted
        hne_path : PathLike
            The directory containing the fovs with their HnE OME-TIFFs.
        seg_dir : PathLike
            The directory containing the segmentation data.
        overlay_channels : str | List[str]
            The overlay channels to be plotted, can be either/both "nuclear_channel" or "membrane_channel",
            defaults to ["nuclear_channel", "membrane_channel"].
        figsize : Tuple[int, int], optional
            The size of the figure, by default (8, 4)
        layout : Literal["constrained", "tight"], optional
            The layout engine, defaults to None, by default None
        image_type : Literal["png", "pdf", "svg"], optional
            The file type to save the plot as, by default "pdf"
        """
        self.fov_name = fov
        self.hne_path = hne_path
        self.seg_dir = seg_dir
        self.overlay_channels = overlay_channels
        self.figsize = figsize
        self.layout = layout
        self.image_type = image_type

        self.fig = plt.figure(figsize=self.figsize, layout=self.layout)

        self.axes = self.fig.subplots(nrows=1, ncols=3, width_ratios=[1, 1, 1])

    def make_plot(self, save_dir: PathLike):
        self.hne_core = imread(
            self.hne_path / self.fov_name / "core.ome.tiff",
            plugin="tifffile",
            is_ome=True,
        )

        self.hne_fov = imread(
            self.hne_path / self.fov_name / "fov.ome.tiff",
            plugin="tifffile",
            is_ome=True,
        )
        self.fov_loc = gpd.read_file(self.hne_path / self.fov_name / "loc.geojson")
        self.fov_overlay = plot_utils.create_overlay(
            fov=self.fov_name,
            segmentation_dir=self.seg_dir + "/deepcell_output",
            data_dir=self.seg_dir + "/deepcell_input",
            img_overlay_chans=self.overlay_channels,
            seg_overlay_comp="whole_cell",
        )

        self.fig.suptitle(
            t=f"{self.fov_name} HnE Core and Cell Segmentation and Overlay", fontsize=8
        )

        self._plot_core()
        self._plot_fov_overlay()

        self.fig.savefig(
            save_dir / f"{self.fov_name}_hne_core_overlay.{self.image_type}"
        )
        plt.close(self.fig)

    def _plot_core(self):
        hne_core_ax = self.axes[0]
        hne_core_ax.set_title(label="HnE Core", fontsize=6)
        hne_core_ax.imshow(X=self.hne_core, aspect="equal", interpolation="none")
        self.fov_loc.buffer(0.1, cap_style=1, join_style=1, resolution=32).plot(
            ax=hne_core_ax,
            facecolor="none",
            edgecolor="black",
            linewidth=1,
            aspect="equal",
        )

        remove_ticks(hne_core_ax, axis="xy")

    def _plot_fov_overlay(self):
        hne_fov_ax = self.axes[1]
        hne_fov_ax.set_title(label="HnE FOV Crop", fontsize=6)
        hne_fov_ax.imshow(X=self.hne_fov, aspect="equal", interpolation="none")

        overlay_ax = self.axes[2]
        overlay_ax.set_title("Nuclear and Membrane Channel Overlay", fontsize=6)
        overlay_ax.imshow(
            X=self.fov_overlay,
            interpolation="none",
            aspect="equal",
        )

        remove_ticks(hne_fov_ax, axis="xy")
        remove_ticks(overlay_ax, axis="xy")

def _set_locator_formatter(ax: Axes, axis: str) -> None:
    """Helper function to remove ticks and formatters for the specified axis."""
    if "x" in axis:
        ax.xaxis.set_major_locator(plt.NullLocator())
        ax.xaxis.set_minor_locator(plt.NullLocator())
        ax.xaxis.set_major_formatter(plt.NullFormatter())
        ax.xaxis.set_minor_formatter(plt.NullFormatter())
    if "y" in axis:
        ax.yaxis.set_major_locator(plt.NullLocator())
        ax.yaxis.set_minor_locator(plt.NullLocator())
        ax.yaxis.set_major_formatter(plt.NullFormatter())
        ax.yaxis.set_minor_formatter(plt.NullFormatter())

def remove_ticks(
    f: Union[Figure, Axes, Iterable[Axes]],
    axis: Literal["x", "y", "xy", "yx"],
) -> None:
    """
    Removes ticks from the axis of a figure or axis object. If a figure is passed,
    the function will remove the axis-ticks of all the figure's axes.

    Args
    -----
    f : Figure | Axes | Iterable[Axes]
        The figure or axis object to remove the ticks from.
    axis : Literal["x", "y", "xy", "yx"]
        The axis to remove the ticks from. If "xy" or "yx" is passed, the function will remove
        the ticks from both axes.

    Raises
    ------
    ValueError
        If f is not a Figure, Axes object, or an iterable of Axes objects.
    """
    if isinstance(f, Figure):
        axes = f.axes
    elif isinstance(f, Axes):
        axes = [f]
    elif isinstance(f, Iterable):
        axes = list(f)  # Ensure iterable is expanded
        if not all(isinstance(a, Axes) for a in axes):
            raise ValueError("All elements of the iterable must be Axes objects.")
    else:
        raise ValueError("f must be a Figure, an Axes object, or an iterable of Axes objects.")

    for ax in axes:
        _set_locator_formatter(ax, axis)

def _remove_x_axis_ticks(ax: plt.Axes) -> None:
    ax.xaxis.set_major_locator(ticker.NullLocator())
    ax.xaxis.set_major_formatter(ticker.NullFormatter())
    ax.xaxis.set_minor_locator(ticker.NullLocator())
    ax.xaxis.set_minor_formatter(ticker.NullFormatter())


def _remove_y_axis_ticks(ax: plt.Axes) -> None:
    ax.yaxis.set_major_locator(ticker.NullLocator())
    ax.yaxis.set_major_formatter(ticker.NullFormatter())
    ax.yaxis.set_minor_locator(ticker.NullLocator())
    ax.yaxis.set_minor_formatter(ticker.NullFormatter())

def _remove_x_axis_ticks(ax: plt.Axes) -> None:
    """Helper function to remove x-axis ticks and labels."""
    ax.xaxis.set_major_locator(plt.NullLocator())
    ax.xaxis.set_minor_locator(plt.NullLocator())
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    ax.xaxis.set_minor_formatter(plt.NullFormatter())

def _remove_y_axis_ticks(ax: plt.Axes) -> None:
    """Helper function to remove y-axis ticks and labels."""
    ax.yaxis.set_major_locator(plt.NullLocator())
    ax.yaxis.set_minor_locator(plt.NullLocator())
    ax.yaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_minor_formatter(plt.NullFormatter())

def _set_locator_formatter(ax: plt.Axes, axis: Literal["x", "y", "xy", "yx"]) -> None:
    """Removes ticks from the specified axis or axes of the provided Axes object."""
    if axis == "x":
        _remove_x_axis_ticks(ax)
    elif axis == "y":
        _remove_y_axis_ticks(ax)
    elif axis in {"xy", "yx"}:
        _remove_x_axis_ticks(ax)
        _remove_y_axis_ticks(ax)
    else:
        raise ValueError("axis must be 'x', 'y', 'xy', or 'yx'")
        
def mask_erosion_ufunc(
    x: xr.DataArray,
    connectivity: int = 2,
    mode: Literal["thick", "inner", "outer", "subpixel"] = "thick",
):
    """_summary_

    Parameters
    ----------
    x : xr.DataArray
        The input label image
    connectivity : int, optional
        The connectivity used to find boundaries, by default 2
    mode : Literal["thick", "inner", "outer", "subpixel"], optional
        How to mark the boundaries , by default "thick"

    Returns
    -------
    NDArray
        The mask of the segmentation with eroded boundaries.
    """
    edges = find_boundaries(
        label_img=x, connectivity=connectivity, mode=mode, background=0
    )
    seg_mask = np.where(edges == 0, x, 0)
    return seg_mask

# functional marker thresholding helpers
class MarkerDict(TypedDict):
    populations: List[str]
    threshold: float
    x_range: Optional[Tuple[float, float]]
    x_ticks: Optional[np.ndarray]
    x_tick_labels: Optional[np.ndarray]

def functional_marker_thresholding(
    cell_table: pd.DataFrame, save_dir: Union[str, pathlib.Path],
    marker_info: Dict[str, MarkerDict], pop_col: str = "cell_cluster",
    figsize: Optional[Tuple[float, float]] = None
):
    """For a set of markers, visualize their distribution across the entire cohort, plus just 
    against the specified populations.

    Args:
        cell_table (pd.DataFrame):
            Cell table with clustered cell populations
        save_dir (Union[str, pathlib.Path]):
            The directory to save the marker distribution histograms
        marker_info (str):
            For each marker, define the populations, threshold, x-range, and x-tick locations
            NOTE: assumes that each marker is being visualized against the same number of
            populations
        pop_col (str):
            Column containing the names of the cell populations
        fig_size (Optional[Tuple[float, float]]):
            The figure size to use for the image.
            If None use default sizing (18.6, 6.6 * len(populations))
    """
    # verify save_dir is valid
    validate_paths([save_dir])

    # verify figsize is valid if set
    if figsize and (len(figsize) != 2 or figsize[0] <= 0 or figsize[1] <= 0):
        raise ValueError(
            "Invalid figsize: it must be in the form (size_x, size_y), size_x > 0, size_y > 0"
        )

    # define the subplots
    markers = list(marker_info.keys())
    figsize = figsize if figsize else (18.6, 6.6 * len(populations))
    fig, axs = plt.subplots(
        len(marker_info),
        len(marker_info[markers[0]]["populations"]) + 1,
        figsize=figsize
    )

    # retrieve all the markers and populations in the cell table (done for validation)
    all_markers: np.ndarray = cell_table.columns.values
    all_populations: np.ndarray = cell_table[pop_col].unique()

    # set axs_row and axs_col as counters to position the titles correctly
    axs_row: int = 0
    axs_col: int = 0

    # iterate over each marker
    for marker in markers:
        # retrieve all populations associated with the marker
        populations: List[str] = marker_info[marker]["populations"]

        # Verify that the marker and all populations specified are valid
        verify_in_list(
            specified_marker=marker,
            cell_table_columns=all_markers
        )

        verify_in_list(
            specified_populations=populations,
            cell_table_populations=all_populations
        )

        # limit x_range to 99.9% of the marker in question if x_range not specified
        x_range = marker_info[marker].get(
            "x_range", (0, np.quantile(cell_table[marker].values, 0.999))
        )

        # retrieve the x ticks and x tick labels
        x_ticks = marker_info[marker].get("x_ticks", None)
        x_tick_labels = marker_info[marker].get("x_tick_labels", None)

        # the first subplot should always be the distribution of the marker against all populations
        threshold: float = marker_info[marker]["threshold"]
        axs[axs_row][0].hist(
            cell_table[marker].values,
            50,
            density=True,
            facecolor='g',
            alpha=0.75,
            range=x_range
        )
        axs[axs_row][0].set_title(
            "{} in all populations".format(marker),
            fontsize=28
        )
        axs[axs_row][0].axvline(x=threshold)

        if isinstance(x_ticks, np.ndarray):
            axs[axs_row][0].set_xticks(x_ticks)

        if isinstance(x_tick_labels, np.ndarray):
            axs[axs_row][0].set_xticklabels(x_tick_labels, fontsize=24)

        axs[axs_row][0].tick_params(axis="y", labelsize=24)

        # add additional subplots to the figure based on the specified populations
        for i, pop in zip(np.arange(1, len(populations) + 1), populations):
            cell_table_marker_sub: pd.DataFrame = cell_table.loc[
                cell_table[pop_col] == pop, marker
            ].values
            axs[axs_row][i].hist(
                cell_table_marker_sub,
                50,
                density=True,
                facecolor='g',
                alpha=0.75,
                range=x_range
            )
            axs[axs_row][i].set_title(
                "{} in {}".format(marker, pop),
                fontsize=28
            )
            axs[axs_row][i].axvline(x=threshold)

            if isinstance(x_ticks, np.ndarray):
                axs[axs_row][i].set_xticks(x_ticks)

            if isinstance(x_tick_labels, np.ndarray):
                axs[axs_row][i].set_xticklabels(x_tick_labels, fontsize=24)

            axs[axs_row][i].tick_params(axis="y", labelsize=24)

        # update axs_row to the next column
        axs_row += 1

    plt.tight_layout()

    # save the figure to save_dir
    fig.savefig(pathlib.Path(save_dir) / f"functional_marker_thresholds.png", dpi=300)

def compute_nonzero_mean_intensity(image_data):
    """Compute the nonzero mean of a specific fov/chan pair.

    Args:
        image_data (numpy.ndarray):
            the image data for a specific fov/chan pair

    Returns:
        float:
            The nonzero mean intensity of the fov/chan pair (`np.nan` if channel contains all 0s)
    """
    # take just the non-zero pixels
    image_data_nonzero = image_data[image_data != 0]

    # take the mean of the non-zero pixels and assign to (fov, channel) in array
    # unless there are no non-zero pixels, in which case default to 0
    if len(image_data_nonzero) > 0:
        nonzero_mean_intensity = image_data_nonzero.mean()
    else:
        nonzero_mean_intensity = 0

    return nonzero_mean_intensity

def compute_total_intensity(image_data):
    """Compute the sum of all pixels of a specific fov/chan pair.

    Args:
        image_data (numpy.ndarray):
            the image data for a specific fov/chan pair

    Returns:
        float:
            The total intensity of the fov/chan pair (`np.nan` if channel contains all 0s)
    """
    return np.sum(image_data)

def compute_99_9_intensity(image_data):
    """Compute the 99.9% pixel intensity value of a specific fov/chan pair.

    Args:
        image_data (numpy.ndarray):
            the image data for a specific fov/chan pair

    Returns:
        float:
            The 99.9% pixel intensity value of a specific fov/chan pair
    """
    return np.percentile(image_data, q=99.9)

def sort_bin_file_fovs(fovs, suffix_ignore=None):
    """Sort a list of fovs in a bin file by fov and scan number.

    fovs (list):
        a list of fovs prefixed with `'fov-m-scan-n'`
    suffix_ignore (str):
        removes this at the end of each fov name, needed if sorting fov-level QC `.csv` files

    Returns:
        list:
            fov name list sorted by ascending fov number, the ascending scan number
    """
    # set suffix_ignore to the empty string if None
    if suffix_ignore is None:
        suffix_ignore = ""

    # TODO: if anyone can do this using a walrus operator I'd appreciate it!
    return sorted(
        fovs,
        key=lambda f: (
            int(f.replace(suffix_ignore, "").split("-")[1]),
            int(f.replace(suffix_ignore, "").split("-")[3]),
        ),
    )

def compute_qc_metrics(
    extracted_imgs_path, fov_name, gaussian_blur=False, blur_factor=1, save_csv=None
):
    """Compute the QC metric matrices for the image data provided.

    Args:
        extracted_imgs_path (str):
            the directory where extracted images are stored
        fov_name (str):
            the name of the FOV to extract from `bin_file_path`, needs to correspond with JSON name
        gaussian_blur (bool):
            whether or not to add Gaussian blurring
        blur_factor (int):
            the sigma (standard deviation) to use for Gaussian blurring
            set to 0 to use raw inputs without Gaussian blurring
            ignored if `gaussian_blur` set to `False`
        save_csv (str):
            path to save csvs of the qc metrics to

    Returns:
        None
    """
    # path validation checks
    if not os.path.exists(extracted_imgs_path):
        raise FileNotFoundError("extracted_imgs_path %s does not exist" % extracted_imgs_path)

    # retrieve the image data from extracted tiff files
    # the image coords should be: ['fov', 'type', 'x', 'y', 'channel']
    image_data = load_utils.load_imgs_from_tree(extracted_imgs_path, fovs=[fov_name])
    image_data = format_img_data(image_data)

    metric_csvs = compute_qc_metrics_direct(image_data, fov_name, gaussian_blur, blur_factor)
    if save_csv:
        for metric_name, data in metric_csvs.items():
            data.to_csv(os.path.join(save_csv, metric_name), index=False)

def compute_qc_metrics_direct(image_data, fov_name, gaussian_blur=False, blur_factor=1):
    """Compute the QC metric matrices for the image data provided.

    Args:
        image_data (xr.DataArray):
            image data in 'extract_bin_files' output format
        fov_name (str):
            the name of the FOV to extract from `bin_file_path`, needs to correspond with JSON name
        gaussian_blur (bool):
            whether or not to add Gaussian blurring
        blur_factor (int):
            the sigma (standard deviation) to use for Gaussian blurring
            set to 0 to use raw inputs without Gaussian blurring
            ignored if `gaussian_blur` set to `False`

    """
    # there's only 1 FOV and 1 type ('pulse'), so subset on that
    image_data = image_data.loc[fov_name, "pulse", :, :, :]

    # define the list of channels to use
    chans = image_data.channel.values

    # define numpy arrays for all the metrics to extract, more efficient indexing than pandas
    blank_arr = np.zeros(image_data.shape[2], dtype="float32")
    nonzero_mean_intensity = copy.deepcopy(blank_arr)
    total_intensity = copy.deepcopy(blank_arr)
    intensity_99_9 = copy.deepcopy(blank_arr)

    # it's faster to loop through the individual channels rather than broadcasting
    for i, chan in enumerate(chans):
        # subset on the channel, cast to float32 to prevent truncation
        image_data_np = image_data.loc[:, :, chan].values.astype(np.float32)

        # STEP 1: gaussian blur (if specified)
        if gaussian_blur:
            image_data_np = gaussian_filter(
                image_data_np, sigma=blur_factor, mode="nearest", truncate=2.0
            )

        # STEP 2: extract non-zero mean intensity
        nonzero_mean_intensity[i] = compute_nonzero_mean_intensity(image_data_np)

        # STEP 3: extract total intensity
        total_intensity[i] = compute_total_intensity(image_data_np)

        # STEP 4: take 99.9% value of the data and assign
        intensity_99_9[i] = compute_99_9_intensity(image_data_np)

    # define the list of numpy arrays for looping
    metric_data = [nonzero_mean_intensity, total_intensity, intensity_99_9]

    metric_csvs = {}

    for ms, md, mc in zip(QC_SUFFIXES, metric_data, QC_COLUMNS):
        # define the dataframe for this metric
        metric_df = pd.DataFrame(columns=["fov", "channel", mc], dtype=object)

        # assign the metric data
        metric_df[mc] = md

        # assign the fov and channel names
        metric_df["fov"] = fov_name
        metric_df["channel"] = chans

        metric_csvs[f"{fov_name}_{ms}.csv"] = metric_df

    return metric_csvs

def combine_qc_metrics(qc_metrics_dir, warn_overwrite=True):
    """Aggregates the QC results of each FOV into one `.csv`.

    Args:
        qc_metrics_dir (str):
            the name of the folder containing the QC metric files
        warn_overwrite (bool):
            whether to warn if existing combined CSV found for each metric
    """
    # path validation check
    if not os.path.exists(qc_metrics_dir):
        raise FileNotFoundError("qc_metrics_dir %s does not exist" % qc_metrics_dir)

    for ms in QC_SUFFIXES:
        # define an aggregated metric DataFrame
        metric_df = pd.DataFrame()

        # list all the files corresponding to this metric
        metric_files = io_utils.list_files(qc_metrics_dir, substrs=ms + ".csv")

        # don't consider any existing combined .csv files, just the fov-level .csv files
        metric_files = [mf for mf in metric_files if "combined" not in mf]

        # sort the files to ensure consistency
        metric_files = sort_bin_file_fovs(metric_files, suffix_ignore="_%s.csv" % ms)

        # iterate over each metric file and append the data to metric_df
        for mf in metric_files:
            metric_df = pd.concat([metric_df, pd.read_csv(os.path.join(qc_metrics_dir, mf))])

        # write the aggregated metric data
        # NOTE: if this combined metric file already exists, it will be overwritten
        if os.path.exists(os.path.join(qc_metrics_dir, "combined_%s.csv" % ms)) and warn_overwrite:
            warnings.warn(
                "Removing previously generated combined %s file in %s" % (ms, qc_metrics_dir)
            )
        metric_df.to_csv(os.path.join(qc_metrics_dir, "combined_%s.csv" % ms), index=False)

def format_img_data(img_data):
    """Formats the image array from load_imgs_from_tree to be same structure as the array returned
    by extract_bin_files. Works for one FOV data at a time.

    Args:
        img_data (str): current image data array as produced by load function
    Returns:
         xarray.DataArray: image data array with shape [fov, type, x, y, channel]
    """
    # add type dimension
    img_data = img_data.assign_coords(type="pulse")
    img_data = img_data.expand_dims("type", 1)

    # edit dimension names
    img_data = img_data.rename({"fovs": "fov", "rows": "x", "cols": "y", "channels": "channel"})

    return img_data

def qc_filtering(qc_metrics: List[str]) -> Tuple[List[str], List[str]]:
    """Filters the QC columns and suffixes based on the user specified QC metrics,
    then sorts the suffixes w.r.t the columns. Refer to `settings.py` for the
    Column and Suffixes available.

    Args:
        qc_metrics (List[str]): A list of QC metrics to use. Options include:

                - `"Non-zero mean intensity"`
                - `"Total intensity"`
                - `"99.9% intensity value"`


    Returns:
        Tuple[List[str], List[str]]: Returns the QC Columns and the QC Suffixes
    """
    # Filter out unused QC columns and suffixes
    if qc_metrics is not None:
        selected_qcs: List[bool] = [qcm in qc_metrics for qcm in QC_COLUMNS]
        qc_cols = list(itertools.compress(QC_COLUMNS, selected_qcs))
        qc_suffixes = list(itertools.compress(QC_SUFFIXES, selected_qcs))
    else:
        qc_cols: List[str] = QC_COLUMNS
        qc_suffixes: List[str] = QC_SUFFIXES

    return qc_cols, qc_suffixes

def _channel_filtering(
    df: pd.DataFrame, channel_include: List[str] = None, channel_exclude: List[str] = None
) -> pd.DataFrame:
    """Filters the DataFrame based on the included and excluded channels. In addition
    the default ignored channels; Au, Fe, Na, Ta, Noodle, are removed.

    Args:
        df (pd.DataFrame): The DataFrame to filter.
        channel_include (List[str], optional): A list of channels to include. Defaults to None.
        channel_exclude (List[str], optional): A list of channels to exclude. Defaults to None.

    Returns:
        pd.DataFrame: The filtered DataFrame.
    """
    if (
        isinstance(channel_include, list)
        and isinstance(channel_exclude, list)
        and not set(channel_exclude).isdisjoint(set(channel_include))
    ):
        raise ValueError("You cannot include and exclude the same channel.")

    # Filter out the default ignored channels
    df = df[~df["channel"].isin(QC_CHANNEL_IGNORE)]

    # Remove the excluded channels
    # If a channel does not exist, it is ignored
    if channel_exclude is not None:
        df: pd.DataFrame = df[~df["channel"].isin(channel_exclude)]

    # Then filter the excluded channels
    # If a channel does not exist, it is ignored
    if channel_include is not None:
        df: pd.DataFrame = df[df["channel"].isin(channel_include)]
    return df

class QCTMA:
    """Computes the QC metrics for a given list of TMAs of interest and saves TMA specific QC files
    in the `qc_tma_metrics_dir` directory.

    Args:
        cohort_path (Union[str,pathlib.Path]): The directory where the extracted images are
            stored.
        qc_tma_metrics_dir (Union[str, pathlib.path]): The directory where to save the QC TMA
            metrics.
        qc_metrics (List[str]): A list of QC metrics to use. Options include:

                - `"Non-zero mean intensity"`
                - `"Total intensity"`
                - `"99.9% intensity value"`

    Attributes:
        qc_cols (List[str]): A list of the QC columns.
        qc_suffixes (List[str]): A list of the QC suffixes, ordered w.r.t `qc_cols`.
        search_term (re.Pattern): The regex pattern to extract n,m from FOV names of the form RnCm.
        tma_avg_zscores (Dict[str, xr.DataArray]): A dictionary containing the average z-scores for
            each TMA for each QC Metric in `qc_metrics`.
    """

    qc_metrics: Optional[List[str]]
    cohort_path: Union[str, pathlib.Path]
    metrics_dir: Union[str, pathlib.Path]

    # Fields initialized after `__post_init__`
    search_term: re.Pattern = field(init=False)
    qc_cols: List[str] = field(init=False)
    qc_suffixes: List[str] = field(init=False)

    # Set by methods
    tma_avg_zscores: Dict[str, xr.DataArray] = field(init=False)

    def __post_init__(self):
        """Initialize QCTMA."""
        # Input validation: Ensure that the paths exist
        io_utils.validate_paths([self.cohort_path, self.metrics_dir])

        self.qc_cols, self.qc_suffixes = qc_filtering(qc_metrics=self.qc_metrics)

        # Create regex pattern for searching RnCm
        self.search_term: re.Pattern = re.compile(r"R\+?(\d+)C\+?(\d+)")

        # Set the tma_avg_zscores to be an empty dictionary
        self.tma_avg_zscores = {}

    def _get_r_c(self, fov_name: pd.Series) -> pd.Series:
        """Extracts the row and column value from a FOV's name containing RnCm.

        Args:
            fov_name (pd.Series): The FOV's name.

        Returns:
            pd.Series: Returns `n` and `m` as a series of integers.
        """
        r, c = map(int, re.search(self.search_term, fov_name).group(1, 2))
        return pd.Series([r, c])

    def _create_r_c_tma_matrix(
        self, group: DataFrameGroupBy, n_cols: int, n_rows: int, qc_col: str
    ) -> pd.Series:
        """Z-scores all FOVS for a given channel and creates a matrix of size `n_rows` by `n_cols`
        as the TMA grid.

        Args:
            group (DataFrameGroupBy): Each group consists of an individual channel, and all of it's
                associated FOVs.
            n_cols (int): The number of columns in the matrix.
            n_rows (int): The number of rows in the matrix.
            qc_col (str): The column to get the the QC data.

        Returns:
            pd.Series[np.ndarray]: Returns the a series containing the z-score matrix.
        """
        rc_array: np.ndarray = np.full(shape=(n_cols, n_rows), fill_value=np.nan)
        rc_array[group["column"] - 1, group["row"] - 1] = stats.zscore(group[qc_col])

        return pd.Series([rc_array])

    def compute_qc_tma_metrics(self, tmas: List[str]):
        """Calculates the QC metrics for a user specified list of TMAs.

        Args:
            tmas (List[str]): The FOVs with the TMA in the folder name to gather.
        """
        with tqdm(
            total=len(tmas), desc="Computing QC TMA Metrics", unit="TMA", leave=True
        ) as tma_pbar:
            for tma in tmas:
                self._compute_qc_tma_metrics(tma=tma)
                tma_pbar.set_postfix(TMA=tma)
                tma_pbar.update(n=1)

    def _compute_qc_tma_metrics(self, tma: str):
        """Computes the FOV QC metrics for all FOVs in a given TMA.
        If the QC metrics have already been computed, then.

        Args:
            tma (str): The TMA to compute the QC metrics for.
        """
        # cannot use `io_utils.list_folders` because it cannot do a partial "exact match"
        # i.e. if we want to match `tma_1_Rn_Cm` but not `tma_10_Rn_Cm`, `io_utils.list_folders`
        # will return both for `tma_1`
        fovs: List[str] = ns.natsorted(
            seq=(p.name for p in pathlib.Path(self.cohort_path).glob(f"{tma}_*"))
        )

        # Compute the QC metrics
        with tqdm(fovs, desc="Computing QC Metrics", unit="FOV", leave=False) as pbar:
            for fov in pbar:
                # Gather the qc tma files for the current fov if they exist
                pre_computed_metrics = filter(
                    lambda f: "combined" not in f,
                    io_utils.list_files(
                        dir_name=self.metrics_dir,
                        substrs=[f"{fov}_{qc_suffix}.csv" for qc_suffix in self.qc_suffixes],
                    ),
                )

                # only compute if any QC files are missing for the current fov
                if len(list(pre_computed_metrics)) != len(self.qc_cols):
                    compute_qc_metrics(
                        extracted_imgs_path=self.cohort_path,
                        fov_name=fov,
                        save_csv=self.metrics_dir,
                    )
                    pbar.set_postfix(FOV=fov, status="Computing")
                else:
                    pbar.set_postfix(FOV=fov, status="Already Computed")

        # Generate the combined metrics for each TMA
        for qc_suffix in self.qc_suffixes:
            metric_files: List[str] = ns.natsorted(
                (
                    io_utils.list_files(
                        dir_name=self.metrics_dir,
                        substrs=[f"{fov}_{qc_suffix}.csv" for fov in fovs],
                    )
                )
            )

            # Define an aggregated metric DataFrame
            combined_metric_tissue_df: pd.DataFrame = pd.concat(
                (pd.read_csv(os.path.join(self.metrics_dir, mf)) for mf in metric_files)
            )

            combined_metric_tissue_df.to_csv(
                os.path.join(self.metrics_dir, f"{tma}_combined_{qc_suffix}.csv"),
                index=False,
            )

    def qc_tma_metrics_zscore(self, tmas: List[str], channel_exclude: List[str] = None):
        """Creates the average zscore for a given TMA across all FOVs and unexcluded channels.
        By default the following channels are excluded: Au, Fe, Na, Ta, Noodle.

        Args:
            tmas (List[str]): The FOVs withmetet the TMA in the folder name to gather.
            channel_exclude (List[str], optional): An optional list of channels to further filter
                out. Defaults to None.
        """
        max_col, max_row = 0, 0
        with tqdm(total=len(tmas), desc="Computing QC TMA Metric Z-scores", unit="TMA") as pbar:
            for tma in tmas:
                self.tma_avg_zscores[tma] = self._compute_qc_tma_metrics_zscore(
                    tma, channel_exclude=channel_exclude
                )
                max_col = max(self.tma_avg_zscores[tma].shape[1], max_col)
                max_row = max(self.tma_avg_zscores[tma].shape[2], max_row)

                pbar.set_postfix(TMA=tma)
                pbar.update()

        # also average z-scores and store
        all_tmas = np.full(
            shape=(len(tmas), len(self.qc_metrics), max_col, max_row), fill_value=np.nan
        )
        for i, tma in enumerate(tmas):
            col, row = self.tma_avg_zscores[tma].shape[1], self.tma_avg_zscores[tma].shape[2]
            all_tmas[i, :, :col, :row] = self.tma_avg_zscores[tma]

        self.tma_avg_zscores["cross_TMA_averages"] = xr.DataArray(
            data=np.stack(np.nanmean(all_tmas, axis=0)),
            coords=[self.qc_cols, np.arange(max_col), np.arange(max_row)],
            dims=["qc_col", "cols", "rows"],
        )

    def _compute_qc_tma_metrics_zscore(
        self,
        tma: str,
        channel_exclude: List[str] = None,
    ) -> xr.DataArray:
        """Creates the average z-score for a given TMA across all FOVs and unexcluded channels.
        By default the following channels are excluded: Au, Fe, Na, Ta, Noodle.

        Args:
            tma (str): The TMA to compute the average z-score for.
            channel_exclude (List[str], optional): An optional list of channels to further filter
                out. Defaults to None.

        Returns:
            xr.DataArray: An xarray DataArray containing the average z-score for each channel across
                a TMA.
        """
        # Sort the loaded combined csv files based on the filtered `qc_suffixes`
        combined_metric_tmas: List[str] = ns.natsorted(
            io_utils.list_files(self.metrics_dir, substrs=f"{tma}_combined"),
            key=lambda tma_mf: (i for i, qc_s in enumerate(self.qc_suffixes) if qc_s in tma_mf),
        )

        zscore_channels_matrix = []
        n_cols: int = None
        n_rows: int = None

        for cmt, qc_col in zip(combined_metric_tmas, self.qc_cols):
            # Open and filter the default ignored channels, along with the user specified channels
            cmt_df: pd.DataFrame = _channel_filtering(
                df=pd.read_csv(os.path.join(self.metrics_dir, cmt)), channel_exclude=channel_exclude
            )

            cmt_df[["column", "row"]] = cmt_df["fov"].apply(self._get_r_c)

            # Get matrix dimensions
            n_cols = cmt_df["column"].max()
            n_rows = cmt_df["row"].max()

            # Z-score all FOVs per channel, and then create the heatmap matrix
            zscore_channel_tmas: pd.DataFrame = cmt_df.groupby(by="channel", sort=True).apply(
                lambda group: self._create_r_c_tma_matrix(group, n_cols, n_rows, qc_col)
            )
            zscore_channel_matrices: np.ndarray = np.array(
                [c_tma[0] for c_tma in zscore_channel_tmas.values],
            )

            avg_zscore = np.mean(zscore_channel_matrices, axis=0)

            zscore_channels_matrix.append(avg_zscore)

        return xr.DataArray(
            data=np.stack(zscore_channels_matrix),
            coords=[self.qc_cols, np.arange(n_cols), np.arange(n_rows)],
            dims=["qc_col", "cols", "rows"],
        )

class QCControlMetrics:
    """Computes QC Metrics for a set of control sample FOVs across various runs, and saves the QC
    files in the `longitudinal_control_metrics_dir`.

    Args:
        cohort_path (Union[str,pathlib.Path]): The directory where the extracted images are
        stored for the control FOVs.
        longitudinal_control_metrics_dir (Union[str, pathlib.Path]): The directory where to save
        the QC Control metrics.
        qc_metrics (List[str]): A list of QC metrics to use. Options include:

                - `"Non-zero mean intensity"`
                - `"Total intensity"`
                - `"99.9% intensity value"`

    Attributes:
        qc_cols (List[str]): A list of the QC columns.
        qc_suffixes (List[str]): A list of the QC suffixes, ordered w.r.t `qc_cols`.
    """

    qc_metrics: Optional[List[str]]
    cohort_path: Union[str, pathlib.Path]
    metrics_dir: Union[str, pathlib.Path]

    # Fields initialized after `__post_init__`
    qc_cols: List[str] = field(init=False)
    qc_suffixes: List[str] = field(init=False)
    longitudinal_control_metrics: Dict[Tuple[str, str], pd.DataFrame] = field(init=False)

    def __post_init__(self):
        """Initialize QCControlMetrics."""
        # Input validation: Ensure that the paths exist
        io_utils.validate_paths([self.cohort_path, self.metrics_dir])

        self.qc_cols, self.qc_suffixes = qc_filtering(qc_metrics=self.qc_metrics)

        self.longitudinal_control_metrics = {}

    def compute_control_qc_metrics(
        self,
        control_sample_name: str,
        fovs: List[str],
        channel_exclude: List[str] = None,
        channel_include: List[str] = None,
    ) -> None:
        """Computes QC metrics for a set of Control Sample FOVs and saves their QC files in the
        `longitudinal_control_metrics_dir`. Calculates the following metrics for the specified
        control samples:
                - `"Non-zero mean intensity"`
                - `"Total intensity"`
                - `"99.9% intensity value"`.

        Args:
            control_sample_name (str): An identifier for naming the control sample.
            fovs (List[str]): A list of control samples to find QC metrics for.
            channel_exclude (List[str], optional): A list of channels to exclude. Defaults to None.
            channel_include (List[str], optional): A list of channels to include. Defaults to None.


        Raises:
            ValueError: Errors if `tissues` is either None, or a list of size 0.
        """
        if fovs is None or not isinstance(fovs, list):
            raise ValueError("The tissues must be specified as a list of strings")

        with tqdm(
            total=len(fovs),
            desc=f"Computing QC Longitudinal Control metrics - {control_sample_name}",
            unit="FOVs",
        ) as pbar:
            for fov in ns.natsorted(fovs):
                # Gather the qc files for the current fov if they exist
                pre_computed_metrics = filter(
                    lambda f: "combined" not in f,
                    io_utils.list_files(
                        dir_name=self.metrics_dir,
                        substrs=[f"{fov}_{qc_suffix}.csv" for qc_suffix in self.qc_suffixes],
                    ),
                )

                if len(list(pre_computed_metrics)) != len(self.qc_cols):
                    compute_qc_metrics(
                        extracted_imgs_path=self.cohort_path,
                        fov_name=fov,
                        save_csv=self.metrics_dir,
                    )
                    pbar.set_postfix(FOV=fov, status="Computing")
                else:
                    pbar.set_postfix(FOV=fov, status="Already Computed")
                pbar.update()

        # Combine metrics for the set of FOVs into a single file per QC metric
        for qc_col, qc_suffix in zip(self.qc_cols, self.qc_suffixes):
            metric_files = filter(
                lambda f: "combined" not in f,
                io_utils.list_files(
                    dir_name=self.metrics_dir,
                    substrs=[f"{fov}_{qc_suffix}.csv" for fov in fovs],
                ),
            )

            # Define an aggregated metric DataFrame, and filter channels
            combined_lc_df: pd.DataFrame = _channel_filtering(
                df=pd.concat(
                    (pd.read_csv(os.path.join(self.metrics_dir, mf)) for mf in metric_files),
                ),
                channel_include=channel_include,
                channel_exclude=channel_exclude,
            )

            self.longitudinal_control_metrics.update(
                {(control_sample_name, qc_col): combined_lc_df}
            )

            combined_lc_df.to_csv(
                os.path.join(self.metrics_dir, f"{control_sample_name}_combined_{qc_suffix}.csv"),
                index=False,
            )

    def transformed_control_effects_data(
        self, control_sample_name: str, qc_metric: str, to_csv: bool = False
    ) -> pd.DataFrame:
        """Creates a transformed DataFrame for the Longitudinal Control effects data, normalizing by the mean,
        then taking the `log2` of each value.

        Args:
            control_sample_name (str): A control sample to tranform the longitudinal control effects for.
            qc_metric (str): The metric to transform.
            to_csv (bool, optional): Whether to save the transformed data to a csv. Defaults to False.

        Returns:
            pd.DataFrame: The transformed QC Longitudinal Control data.
        """
        misc_utils.verify_in_list(user_metric=qc_metric, qc_metrics=self.qc_cols)

        try:
            df: pd.DataFrame = self.longitudinal_control_metrics[control_sample_name, qc_metric]
        except KeyError:
            # A qc file which isn't stored in the longitudinal_control_metrics dictionary, try to load it
            # in if it exists as a file
            df: pd.DataFrame = pd.read_csv(
                os.path.join(
                    self.metrics_dir,
                    f"{control_sample_name}_combined_{self.qc_suffixes[self.qc_cols.index(qc_metric)]}.csv",
                )
            )
        except FileNotFoundError as e:
            raise FileNotFoundError(
                f"QC Metric Not Found for the Control Sample {control_sample_name}"
            ) from e

        # Apply a log2 transformation to the mean normalized data.
        log2_norm_df: pd.DataFrame = df.pivot(
            index="channel", columns="fov", values=qc_metric
        ).transform(func=lambda row: np.log2(row / row.mean()), axis=1)

        mean_log2_norm_df: pd.DataFrame = (
            log2_norm_df.mean(axis=0)
            .to_frame(name="mean")
            .transpose()
            .sort_values(by="mean", axis=1)
        )

        transformed_df: pd.DataFrame = pd.concat(
            objs=[log2_norm_df, mean_log2_norm_df]
        ).sort_values(by="mean", axis=1, inplace=False)

        transformed_df.rename_axis("channel", axis=0, inplace=True)
        transformed_df.rename_axis("fov", axis=1, inplace=True)

        # Save the pivoted dataframe to a csv
        if to_csv:
            qc_suffix: str = self.qc_suffixes[self.qc_cols.index(qc_metric)]
            transformed_df.to_csv(
                os.path.join(
                    self.metrics_dir,
                    f"{control_sample_name}_transformed_{qc_suffix}.csv",
                ),
                index=True,
            )

        return transformed_df

def visualize_qc_metrics(
    metric_name: str,
    qc_metric_dir: Union[str, pathlib.Path],
    save_dir: Union[str, pathlib.Path],
    channel_filters: Optional[List[str]] = ["chan_"],
    axes_font_size: int = 16,
    wrap: int = 6,
    dpi: int = 300,
    return_plot: bool = False,
) -> Optional[sns.FacetGrid]:
    """Visualize the barplot of a specific QC metric.

    Args:
        metric_name (str):
            The name of the QC metric to plot. Used as the y-axis label. Options include:
            `"Non-zero mean intensity"`, `"Total intensity"`, `"99.9% intensity value"`.
        qc_metric_dir (Union[str, pathlib.Path]):
            The path to the directory containing the `'combined_{qc_metric}.csv'` files
        save_dir (Optional[Union[str, pathlib.Path]], optional):
            The name of the directory to save the plot to. Defaults to None.
        channel_filters (List[str], optional):
            A list of channels to filter out.
        axes_font_size (int, optional):
            The font size of the axes labels. Defaults to 16.
        wrap (int, optional):
            The number of plots to display per row. Defaults to 6.
        dpi (Optional[int], optional):
            The resolution of the image to use for saving. Defaults to None.
        return_plot (bool):
            If `True`, this will return the plot. Defaults to `False`

    Raises:
        ValueError:
            When an invalid metric is provided.
        FileNotFoundError:
            The QC metric directory `qc_metric_dir` does not exist.
        FileNotFoundError:
            The QC metric `combined_csv` file is does not exist in `qc_metric_dir`.

    Returns:
        Optional[sns.FacetGrid]: Returns the Seaborn FacetGrid catplot of the QC metrics.
    """
    # verify the metric provided is valid
    if metric_name not in QC_COLUMNS:
        raise ValueError(
            "Invalid metric %s provided, must be set to 'Non-zero mean intensity', "
            "'Total intensity', or '99.9%% intensity value'" % metric_name
        )

    # verify the path to the QC metric datasets exist
    if not os.path.exists(qc_metric_dir):
        raise FileNotFoundError("qc_metric_dir %s does not exist" % qc_metric_dir)

    # get the file name of the combined QC metric .csv file to use
    qc_metric_index = QC_COLUMNS.index(metric_name)
    qc_metric_suffix = QC_SUFFIXES[qc_metric_index]
    qc_metric_path = os.path.join(qc_metric_dir, "combined_%s.csv" % qc_metric_suffix)

    # ensure the user set the right qc_metric_dir
    if not os.path.exists(qc_metric_path):
        raise FileNotFoundError(
            "Could not locate %s, ensure qc_metric_dir is correct" % qc_metric_path
        )

    # read in the QC metric data
    qc_metric_df = pd.read_csv(qc_metric_path)

    # filter out naturally-occurring elements as well as Noodle
    qc_metric_df = qc_metric_df[~qc_metric_df["channel"].isin(QC_CHANNEL_IGNORE)]

    # filter out any channel in the channel_filters list
    if channel_filters is not None:
        qc_metric_df: pd.DataFrame = qc_metric_df[
            ~qc_metric_df["channel"].str.contains("|".join(channel_filters))
        ]

    # catplot allows for easy facets on a barplot
    qc_fg: sns.FacetGrid = sns.catplot(
        x="fov",
        y=metric_name,
        col="channel",
        col_wrap=wrap,
        data=qc_metric_df,
        kind="bar",
        color="black",
        sharex=True,
        sharey=False)

    qc_fg.set_titles(template="{col_name}")
    qc_fg.figure.supxlabel(t="fov", x=0.5, y=0, ha="center", size=axes_font_size)
    qc_fg.figure.supylabel(t=f"{metric_name}", x=0, y=0.5, va="center", size=axes_font_size)
    qc_fg.set(xticks=[])

    qc_fg.set_axis_labels(x_var="", y_var="")
    qc_fg.set_xticklabels([])

    qc_fg.savefig(os.path.join(save_dir, f"{metric_name}_barplot_stats.png"), dpi=dpi)

    if return_plot:
        return qc_fg

def qc_tmas_metrics_plot(
    qc_tmas: QCTMA,
    tmas: List[str],
    save_figure: bool = False,
    dpi: int = 300,
) -> None:
    """Produces the QC TMA metrics plot for a given set of QC metrics applied to a user specified
    TMA. The figures are saved in `qc_tma_metrics_dir/figures`.

    Args:
        qc_tmas (QCTMA): The class which contains the QC TMA data, filepaths, and methods.
        QC matrix.
        tmas (str): The TMAs to plot the QC metrics for.
        save_figure (bool, optional): If `True`, the figure is saved in a subdirectory in the
        `QCTMA.qc_tma_metrics_dir` directory. Defaults to `False`.
        dpi (int, optional): Dots per inch, the resolution of the image. Defaults to 300.
    """
    if save_figure:
        fig_dir: pathlib.Path = pathlib.Path(qc_tmas.metrics_dir) / "figures"
        fig_dir.mkdir(parents=True, exist_ok=True)

    # also plot averages
    tmas.append("cross_TMA_averages")

    with tqdm(total=len(tmas), desc="Plotting QC TMA Metric Z-scores", unit="TMAs") as pbar:
        for tma in tmas:
            _qc_tma_metrics_plot(qc_tmas, tma, fig_dir=fig_dir, save_figure=save_figure, dpi=dpi)
            pbar.set_postfix(TMA=tma)
            pbar.update(n=1)

def _qc_tma_metrics_plot(
    qc_tmas: QCTMA,
    tma: str,
    fig_dir: pathlib.Path,
    save_figure: bool = False,
    dpi: int = 300,
) -> None:
    """Produces the QC TMA metrics plot for a given set of QC metrics applied to a user specified
    TMA. The figures are saved in `qc_tma_metrics_dir/figures`.

    Args:
        qc_tmas (QCTMA): The class which contains the QC TMA data, filepaths, and methods.
        tma (str): The TMA to plot the metrics for.
        fig_dir (pathlib.Path) : Path of where to save the plots.
        save_figure (bool, optional): If `True`, the figure is saved in a subdirectory in the
        `qc_tma_metrics_dir` directory. Defaults to `False`.
        dpi (int, optional): Dots per inch, the resolution of the image. Defaults to 300.
    """
    for qc_metric, suffix in zip(qc_tmas.qc_cols, qc_tmas.qc_suffixes):
        qc_tma_data: np.ndarray = qc_tmas.tma_avg_zscores[tma].loc[qc_metric].values

        # Set up the Figure for multiple axes
        fig: Figure = plt.figure(dpi=dpi)
        qc_tma_data = qc_tma_data.round(decimals=2)

        # Heatmap
        _norm = Normalize(vmin=-1, vmax=1)
        ax_heatmap = fig.add_subplot()
        sns.heatmap(
            data=qc_tma_data,
            square=True,
            ax=ax_heatmap,
            linewidths=1,
            linecolor="black",
            annot_kws={"size": 7},
            cbar_kws={"shrink": 0.5},
            annot=True,
            cmap=sns.color_palette(palette="vlag", as_cmap=True),
            norm=_norm,
        )
        # Set ticks
        ax_heatmap.set_xticks(
            ticks=ax_heatmap.get_xticks(),
            labels=[f"{i+1}" for i in range(qc_tma_data.shape[1])],
            rotation=0,
        )

        ax_heatmap.set_yticks(
            ticks=ax_heatmap.get_yticks(),
            labels=[f"{i+1}" for i in range(qc_tma_data.shape[0])],
            rotation=0,
        )

        ax_heatmap.set_xlabel("Column")
        ax_heatmap.set_ylabel("Row")
        if tma == "cross_TMA_averages":
            ax_heatmap.set_title("Average Scores Across All TMAs\n")
        else:
            ax_heatmap.set_title(f"{tma} - Average {qc_metric}\n")

        if save_figure:
            fig.savefig(
                fname=fig_dir / f"{tma}_{suffix}.pdf",
                dpi=dpi,
                bbox_inches="tight",
            )
            plt.close(fig)

def longitudinal_control_heatmap(
    qc_control: QCControlMetrics,
    control_sample_name: str,
    save_lc_df: bool = False,
    save_figure: bool = False,
    display_figure: bool = False,
    figsize: tuple[int, int] = (12, 12),
    dpi: int = 300,
) -> None:
    """Generates a heatmap of the QC metrics for the QC Control FOVs.

    Args:
        qc_control (QCControlMetrics): The class which contains the QC LC data, filepaths
        , and methods.
        control_sample_name (List[str]): A list of tissues to plot the QC metrics for.
        save_lc_df: (bool, optional): If `True`, the longitudinal control data is saved as a `csv` in
        `qc_control.metrics_dir`. Defaults to `False`.
        save_figure (bool, optional): If `True`, the figure is saved in a subdirectory in the
        `longitudinal_control_metrics_dir` directory. Defaults to `False`.
        display_figure (bool, optional): If `True`, the figure is displayed. Defaults to `False`.
        figsize: (tuple[int, int], optional): The size of the figure. Defaults to (12, 12).
        dpi (int, optional): Dots per inch, the resolution of the image. Defaults to 300.

    Raises:
        ValueError: Raised when the input tissues are not a list of strings.
    """
    if control_sample_name is None or not isinstance(control_sample_name, str):
        raise ValueError("The control sample name must be string.")
    if save_figure:
        fig_dir: pathlib.Path = pathlib.Path(qc_control.metrics_dir) / "figures"
        fig_dir.mkdir(parents=True, exist_ok=True)

    for qc_col, qc_suffix in zip(qc_control.qc_cols, qc_control.qc_suffixes):
        # Try to read the transformed df if it exists
        try:
            t_df: pd.DataFrame = pd.read_csv(
                os.path.join(
                    qc_control.metrics_dir, f"{control_sample_name}_transformed_{qc_suffix}.csv"
                )
            )
            t_df.rename_axis("fov", axis=1, inplace=True)
        # If it doesn't exist, transform the data and save it.
        except FileNotFoundError:
            t_df: pd.DataFrame = qc_control.transformed_control_effects_data(
                control_sample_name=control_sample_name, qc_metric=qc_col, to_csv=save_lc_df
            )

        # Set up the Figure for multiple axes
        fig: Figure = plt.figure(figsize=figsize, dpi=dpi)
        fig.set_layout_engine(layout="constrained")
        gs = gridspec.GridSpec(nrows=2, ncols=1, figure=fig, height_ratios=[len(t_df.index) - 1, 1])

        # Colorbar Normalization
        _norm = Normalize(vmin=-1, vmax=1)
        _cmap = sns.color_palette("vlag", as_cmap=True)

        fig.suptitle(f"{control_sample_name} - QC: {qc_col}")

        # Annontation kwargs
        annotation_kws = {
            "horizontalalignment": "center",
            "verticalalignment": "center",
            "fontsize": 8,
        }

        # Heatmap
        ax_heatmap: Axes = fig.add_subplot(gs[0, 0])

        sns.heatmap(
            t_df[~t_df.index.isin(["mean"])],
            ax=ax_heatmap,
            linewidths=1,
            linecolor="black",
            cbar_kws={"shrink": 0.5},
            annot=True,
            annot_kws=annotation_kws,
            xticklabels=False,
            norm=_norm,
            cmap=_cmap,
        )

        # cbar title
        ax_heatmap.collections[0].colorbar.ax.set_title(r"$\log_2(QC)$")

        # Axes labels, and ticks
        ax_heatmap.set_yticks(
            ticks=ax_heatmap.get_yticks(),
            labels=ax_heatmap.get_yticklabels(),
            rotation=0,
        )
        ax_heatmap.set_xlabel(None)

        # Averaged values
        ax_avg: Axes = fig.add_subplot(gs[1, 0])

        sns.heatmap(
            data=t_df[t_df.index.isin(["mean"])],
            ax=ax_avg,
            linewidths=1,
            linecolor="black",
            annot=True,
            annot_kws=annotation_kws,
            fmt=".2f",
            cmap=ListedColormap(["white"]),
            cbar=False,
        )
        ax_avg.set_yticks(
            ticks=ax_avg.get_yticks(),
            labels=["Mean"],
            rotation=0,
        )
        ax_avg.set_xticks(
            ticks=ax_avg.get_xticks(),
            labels=ax_avg.get_xticklabels(),
            rotation=45,
            ha="right",
            rotation_mode="anchor",
        )
        ax_heatmap.set_ylabel("Channel")
        ax_avg.set_xlabel("FOV")

        # Save figure
        if save_figure:
            fig.savefig(
                fname=pathlib.Path(qc_control.metrics_dir)
                / "figures"
                / f"{control_sample_name}_heatmap_{qc_suffix}.png",
                dpi=dpi,
                bbox_inches="tight",
            )
        if display_figure:
            fig.show()
        else:
            plt.close(fig)

def plot_unstructured_niche(df, palette, figsize, hue, group, group_id, legend_title, save_directory, filename_save, axes_label = True, ext = '.pdf', legend = True):
    sns.set_style('ticks')
    sns.set_context('paper')
    plt.figure(figsize=figsize)
    g = sns.scatterplot(x = 'X0', y = 'Y0', data = df, linewidth = 0, s = 2, color = '#4D4D4D', label = 'random', alpha = 0.3)
    g = sns.scatterplot(x = 'X0', y = 'Y0', data = df[np.isin(df[group], group_id)], linewidth = 0.1, s = 5, marker = 'o', hue = hue, palette = palette)
    if legend == True:
        plt.legend(title=legend_title, bbox_to_anchor=(0.98, 1), loc='upper left', fontsize = 8)
    else:
        g.legend_.remove()
    
    if axes_label ==True:
        g.set_ylim(-25, 2100)
        g.set_xlim(-25, 2100)
    else:
        g.axis("off")
    g.tick_params(labelsize=10)
    # Removing the spines for a cleaner look
    sns.despine()
    if save_directory is not None:
        plt.savefig(os.path.join(save_directory, filename_save+ext), bbox_inches = 'tight')
        plt.close()
    else:
        plt.show()
        plt.close()

def plot_unstructured_niche_cat(df, figsize, hue, legend_title, save_directory, filename_save, axes_label = True, ext = '.pdf', legend = True):
    sns.set_style('ticks')
    sns.set_style('ticks')
    sns.set_context('paper')
    plt.figure(figsize=figsize)
    g = sns.scatterplot(x = 'X0', y = 'Y0', data = df, linewidth = 0.1, s = 5, marker = 'o', hue = hue, palette = 'Set2')
    if legend== True:
        plt.legend(title=legend_title, bbox_to_anchor=(0.98, 1), loc='upper left', fontsize = 8)
    else:
       g.legend_.remove()
    
    if axes_label ==True:
        g.set_ylim(-25, 2100)
        g.set_xlim(-25, 2100)
    else:
        g.axis("off")
    g.tick_params(labelsize=10)
    # Removing the spines for a cleaner look
    sns.despine()
    if save_directory is not None:
        plt.savefig(os.path.join(save_directory, filename_save+ext), bbox_inches = 'tight')
        plt.close()
    else:
        plt.show()
        plt.close()

def plot_unstructured_niche_score(df, cmap, figsize, hue, group, group_id, legend_title, save_directory, filename_save, vmin, vcenter, vmax, axes_label = True, cbar = True, ext = '.pdf'):
    sns.set_style('ticks')
    sns.set_context('paper')
    from matplotlib.colors import TwoSlopeNorm
    # Create a figure and axis object
    fig, ax = plt.subplots(figsize=figsize)

    # Plot all data points with a base color
    sns.scatterplot(x='X0', y='Y0', data=df, ax=ax, linewidth=0, s=2, color='#4D4D4D', label='random', alpha = 0.3)
    ax.legend_.remove()

    # Set up normalization with a center at 0
    subset = df[np.isin(df[group], group_id)]
    norm = TwoSlopeNorm(vmin=vmin, vmax=vmax, vcenter=vcenter)

    # Plot specific group with continuous hue and centered norm
    scatter = ax.scatter(x=subset['X0'], y=subset['Y0'], c=subset[hue], s=5, marker='o', 
                         cmap=cmap, norm=norm, linewidth=0.1)

    if cbar == True:
        # Adding color bar for continuous hue with centered norm
        cbar = fig.colorbar(scatter, ax=ax, label=legend_title, aspect=40, pad=0.02)
        cbar.set_label(legend_title, fontsize=10)

    if axes_label == True:
        ax.set_ylim(-25, 2100)
        ax.set_xlim(-25, 2100)
    else:
        ax.axis("off")
    
    ax.set_xlabel('u', fontsize=12)
    ax.set_ylabel('v', fontsize=12)
    ax.tick_params(labelsize=10)

    # Removing the spines for a cleaner look
    sns.despine()

    if save_directory is not None:
        plt.savefig(os.path.join(save_directory, filename_save + ext), bbox_inches='tight')
        plt.close()
    else:
        plt.show()
        plt.close()

def plot_niche_score(adata_sub, fov, seg_dir, niche_key, niche, metric = 'spatialFDR', vmin = 0, vmax = 10, fontsize = 12, cmap = 'vlag', background = [0.3, 0.3, 0.3, 1],figsize = (6, 6), save_directory = None, filename_save = None, cbar = True, ext = '.pdf'):
    import matplotlib.cm as cm
    seg_img = io.imread(os.path.join(seg_dir, fov +'_whole_cell.tiff'))
    seg_img = np.squeeze(seg_img)
    seg_img[find_boundaries(seg_img, mode='inner')] = 0

    cell_table_subset = adata_sub.obs.loc[:, ['fov', 'label', niche_key, metric]]
    score = cell_table_subset[~cell_table_subset[metric].isnull()]
    score = cell_table_subset[np.isin(cell_table_subset[niche_key], niche)]
    score_df = pd.DataFrame(score[metric].values, index = score.label.astype('int'))

    # Initialize the score array with a default value (e.g., NaN or 0)
    scores = np.full(seg_img.shape, np.nan, dtype=np.float32)

    # Iterate over unique cell IDs (excluding the background, usually 0)
    for cell_id in tqdm(np.unique(seg_img)[1:]):
        m = seg_img == cell_id  # Create a mask for the current cell_id
        if cell_id in score_df.index:
            scores[m] = score_df.loc[cell_id].values[0]
        else:
            scores[m] = np.nan  # Assign NaN to areas with no value for a grey color later

    # Use the vlag colormap without normalization
    cmap_ = cm.get_cmap(cmap)

    
    norm = Normalize(vmin=vmin, vmax=vmax)

    # Apply the colormap to the normalized scores
    colored_image = cmap_(norm(scores))

    # Handle NaN values to make them grey
    nan_mask = np.isnan(scores)
    colored_image[nan_mask] = background  # Assign grey color for NaN

    # Handle seg_img = 0 entries to make them black
    black_mask = seg_img == 0
    colored_image[black_mask] = [0, 0, 0, 1]  # Assign black color for seg_img == 0

    # Create the figure and axis for plotting
    fig, ax = plt.subplots(figsize=figsize)
    im = ax.imshow(colored_image)
    ax.axis("off")
    if cbar == True:
        # Adding color bar for continuous hue with centered norm
        cbar = fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax, aspect=40, pad=0.02)
        cbar.set_label(metric, fontsize=10)

        cbar.outline.set_edgecolor('black')  # Set the border color to black
        cbar.outline.set_linewidth(1)  # Set the border width
        cbar.ax.tick_params(labelsize=fontsize)
    if filename_save is not None:
        plt.savefig(os.path.join(save_directory, filename_save+ext), bbox_inches = 'tight')
    else:
        plt.show()
    plt.close()

def load_and_plot_image(ax, filepath, title):
    if os.path.exists(filepath):
        img = mpimg.imread(filepath)
        ax.imshow(img)
        ax.axis("off")
        ax.set_title(title, fontsize=10)
    else:
        ax.axis("off")
        ax.set_title(f"Missing: {title}", fontsize=10)

def move_files_to_tmp(file_list, save_directory, tmp_folder):
    """
    Move files from save_directory to a unique temporary folder.

    Args:
    - file_list: List of file paths relative to save_directory.
    - save_directory: Base directory where files are located.
    - tmp_folder: Base folder to move the files into.

    Returns:
    - None
    """
    # Ensure tmp_folder is unique
    if os.path.exists(tmp_folder):
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        tmp_folder = f"{tmp_folder}_{timestamp}"
    
    # Create the unique tmp folder
    os.makedirs(tmp_folder, exist_ok=True)

    # Move each file to the tmp folder
    for file in file_list:
        full_path = os.path.join(save_directory, file)
        if os.path.exists(full_path):
            tmp_path = os.path.join(tmp_folder, os.path.basename(full_path))
            shutil.move(full_path, tmp_path)  # Move the file
            print(f"Moved: {full_path} -> {tmp_path}")
        else:
            print(f"File does not exist and was skipped: {full_path}")
           
def compute_logFC(adata, patient_key, cell_type_key, outcome_key, outcome_dict, id1, id2):
    counts = adata.obs.groupby([patient_key, cell_type_key]).size().unstack().copy()
    norm_counts = counts.div(counts.sum(axis = 1), axis = 0).reset_index()
    norm_counts[outcome_key] = pd.Series(norm_counts[patient_key]).map(outcome_dict)
    norm_counts.drop(columns = patient_key, inplace = True)
    med_counts = norm_counts.groupby([outcome_key]).mean()
    logFC = np.log2(med_counts.loc[id1, :].div(med_counts.loc[id2, :], axis = 0))
    return logFC

def compute_kmeans(adata_niche,
                    unidentified_idx,
                    n_clusters = 2,
                    random_state = 0,
                    key_added = 'kmeans_cluster'):
    from sklearn.cluster import KMeans
    bool_idx = ~np.isin(adata_niche.obs_names, unidentified_idx)
    df = adata_niche[~np.isin(adata_niche.obs_names, unidentified_idx)].to_df()
    cluster_fit = KMeans(n_clusters=n_clusters, random_state=random_state, n_init='auto').fit(df)
    cluster_pred = cluster_fit.predict(df)
    cluster_pred = cluster_pred + 1
    cluster_df = pd.DataFrame(cluster_pred, index = df.index, columns = [key_added])
    if key_added in adata_niche.obs.columns:
        adata_niche.obs.drop(columns = key_added, inplace = True)

    adata_niche.obs = pd.merge(adata_niche.obs, cluster_df, left_index = True, right_index = True, how = 'left')
    adata_niche.obs[key_added][bool_idx] = np.round(adata_niche.obs[key_added][bool_idx])
    adata_niche.obs[key_added][~bool_idx] = np.nan
    adata_niche.obs[key_added] = pd.Categorical(adata_niche.obs[key_added])
    return adata_niche, cluster_fit

def boxplot(mdata,
            feature_key="quiche",
            alpha: float = 0.05,
            niches=None,
            figsize=(6, 12),
            annot_key='quiche_niche',
            xlim=None,
            fontsize=8,
            colors=['#377eb8', '#e41a1c'],
            save_directory='figures',
            colors_dict=None,
            circle_x = -3.8,
            circle_width = 0.23, 
            filename_save=None):
    
    try:
        nhood_adata = mdata[feature_key].T.copy()
    except KeyError:
        raise RuntimeError(
            "mdata should be a MuData object with two slots: feature_key and 'milo'. Run 'milopy.count_nhoods(adata)' first."
        ) from None

    if niches is not None:
        nhood_adata = nhood_adata[np.isin(nhood_adata.obs[annot_key], niches)]

    try:
        nhood_adata.obs[annot_key]
    except KeyError:
        raise RuntimeError(f"{annot_key} not defined") from None

    try:
        nhood_adata.obs["logFC"]
    except KeyError:
        raise RuntimeError("Unable to find 'logFC' in mdata.uns['nhood_adata'].obs. Run 'core.da_nhoods(adata)' first.") from None

    sorted_annos = (
        nhood_adata.obs[[annot_key, "logFC"]].groupby(annot_key).median().sort_values("logFC", ascending=True).index
    )
    
    sns.set_style('ticks')
    anno_df = nhood_adata.obs[[annot_key, "logFC", "SpatialFDR"]].copy()
    anno_df["is_signif"] = anno_df["SpatialFDR"] < alpha
    anno_df = anno_df[anno_df[annot_key] != "nan"]

    cmap_df = pd.DataFrame(mdata[feature_key].var.groupby(annot_key)['logFC'].mean(), columns=['logFC'])
    cmap = np.full(np.shape(mdata[feature_key].var.groupby(annot_key)['logFC'].mean())[0], 'lightgrey', dtype='object')
    cmap[mdata[feature_key].var.groupby(annot_key)['logFC'].mean() <= -1] = colors[0]
    cmap[mdata[feature_key].var.groupby(annot_key)['logFC'].mean() > 1] = colors[1]
    cmap_df['cmap'] = cmap
    
    fig = plt.figure(figsize=figsize)
    gs = GridSpec(1, 1)
    
    ax0 = plt.subplot(gs[0])

    sns.boxplot(
        data=anno_df,
        y=annot_key,
        x="logFC",
        order=sorted_annos,
        orient="h",
        palette=cmap_df.loc[sorted_annos]['cmap'].values,
        ax=ax0,
        width=0.6,
        fliersize=0.8
    )

    ax0.set_xlabel('Log2(FC) Abundance', fontsize=fontsize)
    ax0.set_ylabel('Cell Type Neighborhoods', fontsize=fontsize)
    if xlim is not None:
        ax0.set_xlim(xlim[0], xlim[1])

    y_ticks = ax0.get_yticks()
    y_tick_labels = [tick.get_text() for tick in ax0.get_yticklabels()]
    ax0.set_yticklabels(y_tick_labels, fontsize=fontsize, ha='right', position=(-0.035, 0))

    ellipse_size = 0.4

    for y, cell_type in zip(y_ticks, y_tick_labels):
        if cell_type in colors_dict:
            color = colors_dict[cell_type]
            ellipse = Ellipse(
                (circle_x, y), width=circle_width, height=ellipse_size, color=color, transform=ax0.transData,
                clip_on=False, zorder=10
            )
            ax0.add_patch(ellipse)

    plt.tick_params(labelsize=fontsize)

    ax0.axvline(x=0, ymin=0, ymax=1, color="black", linestyle="--", linewidth=1)
    plt.subplots_adjust(wspace=0.12)

    if filename_save is not None:
        plt.savefig(os.path.join(save_directory, filename_save + '.pdf'), bbox_inches='tight')

    plt.show()

def plot_infiltrate(adata, cluster_key, patient_key, cell_type_list, cutoff, figsize, save_directory, filename_save):
    """
    Plots a heatmap of infiltration with chi-square test annotations for column-wise comparisons.

    Parameters:
        adata: AnnData object.
        cluster_key: Key in adata.obs for cell clusters.
        patient_key: Key in adata.obs for patient IDs.
        cell_type_list: List of cell types to include.
        cutoff: Cutoff for binary infiltration (0 or 1).
        figsize: Tuple for figure size.
        save_directory: Directory to save the figure.
        filename_save: Name of the saved file (without extension).
    """
    from scipy.stats import chisquare
    # Group data by patient and cluster, then filter for specified cell types
    df = adata.obs.groupby([patient_key, cluster_key]).size().unstack().loc[:, cell_type_list]
    df[df < cutoff] = 0
    df[df.isna()] = 0
    df[df >= cutoff] = 1

    chi_square_results = {}
    for i in range(len(cell_type_list) - 1):
        col1 = cell_type_list[i]
        col2 = cell_type_list[i + 1]

        # Count occurrences for each column
        counts_col1 = df[col1].value_counts()
        counts_col2 = df[col2].value_counts()

        if len(counts_col1) > 1 and len(counts_col2) > 1:
            aligned_counts = pd.DataFrame({'col1': counts_col1, 'col2': counts_col2})#.fillna(0)
            _, p_value = chisquare(f_obs=aligned_counts['col1'], f_exp=aligned_counts['col2'])
            chi_square_results[(col1, col2)] = p_value
        else:
            chi_square_results[(col1, col2)] = None

    fig, ax = plt.subplots(1, 1, sharex=True, sharey=True, figsize=figsize, gridspec_kw={'wspace': 0, 'bottom': 0.15})
    g = sns.heatmap(df.sort_values(by=cell_type_list, ascending=True), vmin=0, vmax=1, cmap='cividis', yticklabels=False, ax=ax, linewidths=0)
    g.set_ylabel(f'Patients (N = {len(df)})', fontsize=12)
    g.set_xlabel('')
    g.set_xticklabels(g.get_xticklabels(), rotation=45, horizontalalignment='right')

    # Annotate significant chi-square test results
    for i in range(len(cell_type_list) - 1):
        col1 = cell_type_list[i]
        col2 = cell_type_list[i + 1]
        p_value = chi_square_results.get((col1, col2))
        if p_value is not None and p_value < 0.05:  # Significant
            ax.text(i + 1, -5, '*', ha='center', va='center', fontsize=15, color='k', transform=ax.transData)

    plt.savefig(os.path.join(save_directory, filename_save + '.pdf'), bbox_inches='tight')

def plot_infiltrate_grid(adata, cluster_key, patient_key, cell_type_list, outcome_key, cutoff, figsize, save_directory, filename_save):
    """
    Plots a 2x2 grid showing infiltration heatmaps for outcome groups side by side.
    
    Parameters:
        adata: AnnData object.
        cluster_key: Key in adata.obs for cell clusters.
        patient_key: Key in adata.obs for patient IDs.
        cell_type_list: List of cell types to include.
        outcome_key: Key in adata.obs for the outcome variable.
        cutoff: Cutoff for binary infiltration (0 or 1).
        figsize: Tuple for figure size.
        save_directory: Directory to save the figure.
        filename_save: Name of the saved file (without extension).
    """
    df = adata.obs.groupby([patient_key, cluster_key]).size().unstack().loc[:, cell_type_list]
    df[df < cutoff] = 0
    df[df.isna()] = 0
    df[df >= cutoff] = 1

    outcomes = adata.obs[outcome_key].unique()
    n_outcomes = len(outcomes)
    n_rows, n_cols = 1, 2 
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=figsize, gridspec_kw={'wspace': 0.1, 'hspace': 0.2})
    axes = axes.flatten()

    for i, outcome in enumerate(outcomes):
        if i >= len(axes):
            break

        df_outcome = df.loc[adata.obs[adata.obs[outcome_key] == outcome][patient_key].unique()]

        sns.heatmap(
            df_outcome.sort_values(by=cell_type_list, ascending=True),
            vmin=0, vmax=1, cmap='cividis', yticklabels=False, ax=axes[i],linewidths=0
        )
        axes[i].set_title(f"{outcome}", fontsize=12)
        axes[i].set_ylabel(f"Patients (N = {len(df_outcome)})", fontsize=10)
        axes[i].set_xlabel('')

    for j in range(i + 1, len(axes)):
        fig.delaxes(axes[j])
    
    for ax in axes[:n_outcomes]:
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45, horizontalalignment='right')
    
    plt.savefig(os.path.join(save_directory, filename_save + '.pdf'), bbox_inches='tight')

def plot_correlation(mdata, save_directory, filename_save):
    nsig = mdata['quiche'].var[['SpatialFDR', 'Patient_ID']][mdata['quiche'].var['SpatialFDR'] < 0.05].groupby('Patient_ID').size()
    nsig.index = nsig.index.astype('str')
    counts = mdata['expression'].obs.groupby('Patient_ID').size()
    counts.index =counts.index.astype('str')

    test_df = pd.concat([nsig, counts], axis = 1)
    test_df.columns = ['number_hits', 'number_cells']

    correlation = test_df['number_cells'].corr(test_df['number_hits'])

    plt.figure(figsize=(4, 4))
    g = sns.scatterplot(x='number_cells', y='number_hits', data=test_df, s = 50, color = '#A68CA6')

    plt.text(0.95, 0.95, f'r = {correlation:.2f}', 
            horizontalalignment='right', verticalalignment='top', 
            transform=plt.gca().transAxes, fontsize=12, color='k')

    plt.xlabel('Number of cells', fontsize = 12)
    plt.ylabel('Number of significant \n QUICHE niche cells', fontsize = 12)
    g.tick_params(labelsize=10)
    plt.savefig(os.path.join(save_directory, filename_save+'.pdf'), bbox_inches = 'tight')

def plot_nhood_distances(mdata, niche_list, save_directory, filename_save):
    df = mdata['quiche'].var[['quiche_niche', 'cell_cluster']].copy()
    knn_mat = mdata['spatial_nhood'].obsp['connectivities']
    dist_mat = mdata['spatial_nhood'].obsp['distances']

    # Extract the data from the sparse matrix
    row_indices, col_indices = dist_mat.nonzero()
    distances = dist_mat.data
    nonzero_mask = distances > 0
    distances = distances[nonzero_mask]

    # Create a DataFrame to store the distances and corresponding labels
    df_distances = pd.DataFrame({
        'cell_index': row_indices,
        'neighbor_index': col_indices,
        'distance': distances,
        'cell_label': list(df['quiche_niche'].values[row_indices]),
        'neighbor_label': list(df['quiche_niche'].values[col_indices])
    })
    df_distances['same_label'] = df_distances['cell_label'] == df_distances['neighbor_label']
    intra_niche_means = []
    inter_niche_means = []
    niches = []

    for niche in niche_list:
        # Filter distances for cells belonging to the current quiche_niche
        df_niche = df_distances[df_distances['cell_label'] == niche]

        # Calculate mean intra-niche and inter-niche distances
        mean_intra_niche_distance = df_niche[df_niche['same_label']].groupby('cell_index')['distance'].mean().mean()
        mean_inter_niche_distance = df_niche[~df_niche['same_label']].groupby('cell_index')['distance'].mean().mean()
        
        # Append the results to the lists
        intra_niche_means.append(mean_intra_niche_distance)
        inter_niche_means.append(mean_inter_niche_distance)
        niches.append(niche)

    df_plot = pd.DataFrame({'Quiche Niche': niches,
        'Intra-Niche Distance': intra_niche_means,
        'Inter-Niche Distance': inter_niche_means})

    df_plot_nonans = df_plot[~df_plot['Intra-Niche Distance'].isnull()].copy()

    df_melted = df_plot_nonans.melt(id_vars='Quiche Niche', value_vars=['Intra-Niche Distance', 'Inter-Niche Distance'],
                            var_name='Distance Type', value_name='Mean Distance')

    plt.figure(figsize=(4, 6))
    sns.boxplot(x='Distance Type', y='Mean Distance', data=df_melted, width = 0.6, fliersize = 0, palette='Set2')
    sns.stripplot(x='Distance Type', y='Mean Distance', data=df_melted, jitter=False, color='black', size=4)
    for i in range(len(df_plot)):
        plt.plot([0, 1], [df_plot['Intra-Niche Distance'][i], df_plot['Inter-Niche Distance'][i]], 'k-', lw=0.5, alpha=0.2)
    plt.ylabel('Mean Distance')
    plt.xlabel('')
    plt.savefig(os.path.join(save_directory, filename_save+'.pdf'), bbox_inches = 'tight')

def plot_prop(mdata, niche_list, save_directory, filename_save):
    sns.set_style('ticks')
    df_prop = mdata['spatial_nhood'].to_df()
    df_prop['quiche_niche'] =  mdata['spatial_nhood'].obs['quiche_niche']
    average_proportions = df_prop.groupby('quiche_niche').mean()
    n_plots = len(niche_list)
    n_rows = 10 
    n_cols = 5 
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(20, 30))
    axes = axes.flatten()
    for i, niche in enumerate(niche_list):
        if i >= n_plots:
            break
        average_proportions.loc[niche].plot(kind='bar', color='skyblue', ax=axes[i])
        axes[i].tick_params(axis='x', rotation=90, labelsize=7)
        axes[i].set_title(f'{niche}')
        axes[i].set_xlabel('cell_cluster')
        axes[i].set_ylabel('Average proportion')
        axes[i].tick_params(axis='x', rotation=45)
        axes[i].set_xticklabels(axes[i].get_xticklabels(), ha='right')
    for j in range(i + 1, len(axes)):
        fig.delaxes(axes[j])

    plt.tight_layout()
    plt.savefig(os.path.join(save_directory, filename_save+'.pdf'), bbox_inches = 'tight')

def compute_abundance(adata, patient_key, cell_type_key, outcome_key, outcome_dict, id1, id2):
    counts = adata.obs.groupby([patient_key, cell_type_key]).size().unstack().copy()
    norm_counts = counts.div(counts.sum(axis=1), axis=0).reset_index()
    norm_counts[outcome_key] = pd.Series(norm_counts[patient_key]).map(outcome_dict)

    cell_type_columns = counts.columns
    results = {cell_type_key: [], 'p_value': [], 'Log2FC': []}
    
    for cell_type in cell_type_columns:
        group1 = norm_counts[norm_counts[outcome_key] == id1][cell_type]
        group2 = norm_counts[norm_counts[outcome_key] == id2][cell_type]

        mean_group1 = group1.mean()
        mean_group2 = group2.mean()
        if mean_group2 > 0:
            log2fc = np.log2((mean_group1 + 1e-10) / (mean_group2 + 1e-10))
        else:
            log2fc = np.nan    
        
        results['Log2FC'].append(log2fc)
        
        if len(group1.unique()) > 1 and len(group2.unique()) > 1:
            try:
                _, p_value = ranksums(group1, group2)
                results[cell_type_key].append(cell_type)
                results['p_value'].append(p_value)
            except ValueError:
                results[cell_type_key].append(cell_type)
                results['p_value'].append(np.nan)
        else:
            results[cell_type_key].append(cell_type)
            results['p_value'].append(np.nan)
    p_values_df = pd.DataFrame(results)
    fdr_corrected = multipletests(results['p_value'], method='fdr_bh')[1]
    p_values_df['FDR_p_value'] = fdr_corrected
    p_values_df['-log10(Adj. p-value)'] = -1*np.log10(fdr_corrected)
    return norm_counts, p_values_df

def plot_cell_type_abundance_grid(norm_counts, p_values_df, outcome_key, cell_type_key, fdr_column='FDR_p_value', threshold = 0.05, order = [0, 1], save_directory = 'figures', filename_save = 'boxplot'):
    cell_types = p_values_df[cell_type_key].tolist()
    fdr_p_values = p_values_df[fdr_column].tolist()
    n_cell_types = len(cell_types)
    n_cols = 7 
    n_rows = (n_cell_types + n_cols - 1) // n_cols

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(20, n_rows * 2.5))
    axes = axes.flatten()
    for i, (cell_type, fdr_p_value) in enumerate(zip(cell_types, fdr_p_values)):
        ax = axes[i]
        g = sns.boxplot(data=norm_counts, x=outcome_key, y=cell_type, ax=ax, width = 0.5, palette="Set2", fliersize=0, order=order)
        g = sns.stripplot(data=norm_counts, x=outcome_key, y=cell_type, ax=ax, palette="Set2", linewidth=0.2, edgecolor='gray', order=order)
        g.tick_params(labelsize=10)
        ax.set_title(f"{cell_type}", fontsize = 12)
        ax.set_xlabel("Outcome", fontsize = 12)
        ax.set_ylabel("Normalized abundance", fontsize = 12)
        fdr_text = f"FDR={fdr_p_value:.2g}"
        ax.text(0.75, 0.9, fdr_text,
            ha='center', va='center', transform=ax.transAxes,
            fontsize=10, color="red" if fdr_p_value >= -1*np.log10(threshold) else "black")
    for j in range(i + 1, len(axes)):
        fig.delaxes(axes[j])
    plt.tight_layout()
    plt.savefig(os.path.join(save_directory, filename_save+'.pdf'), bbox_inches = 'tight', dpi = 400)

def plot_auc(auc_score, save_directory = 'figures', filename_save = 'auc'):
    sns.set_style('ticks')
    _, axes = plt.subplots(1,1, figsize = (4,4), gridspec_kw={'hspace': 0.65, 'wspace': 0.3, 'bottom':0.15}, dpi = 600)
    g = sns.boxplot(x = 'Label', y = 'AUC', data = auc_score, linewidth = 1,fliersize=0, width = 0.6, ax = axes, palette = ['#55D6BE', 'lightgrey'])

    g = sns.stripplot(x = 'Label', y = 'AUC', data = auc_score, linewidth = 0.8, 
                        size=5, edgecolor="black", jitter = True, ax = axes, palette = ['#55D6BE', 'lightgrey'])
    g.tick_params(labelsize=14)
    g.set_xlabel('', fontsize = 16)
    g.set_ylabel('AUC', fontsize = 16)
    plt.ylim(0, 1.0)
    plt.savefig(os.path.join(save_directory, filename_save+'.pdf'), bbox_inches = 'tight')