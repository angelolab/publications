import json
import math
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter, MaxNLocator
import natsort as ns
import numpy as np
import os
import pandas as pd
import pathlib
from scipy.ndimage import gaussian_filter
import seaborn as sns
import skimage.io as io
from skimage.measure import label
from skimage import morphology
import xarray as xr
import python_files.utils as cancer_mask_utils

from os import PathLike
from PIL import Image, ImageDraw, ImageFont
from typing import Dict, List, Optional, Tuple, TypedDict, Union, Literal
from skimage.io import imread
import geopandas as gpd
from ark.utils import plot_utils, data_utils
from alpineer.io_utils import list_folders, list_files, remove_file_extensions, validate_paths
from alpineer.load_utils import load_imgs_from_tree, load_imgs_from_dir
from alpineer.misc_utils import verify_in_list
from .utils import  QuantileNormalization
from python_files.utils import compare_populations

# import multipletests library from statsmodels
from statsmodels.stats.multitest import multipletests

# from .utils import remove_ticks,

ACQUISITION_ORDER_INDICES = [
    11, 12, 13, 14, 15, 17, 18, 20, 22, 23, 24, 28, 29, 30, 31, 32, 33, 34, 35,
    36, 39, 40, 41, 42, 43, 44, 45, 46, 47
]
MEAN_PANEL_NORM = {
    "CD11c": 0.006193546577196089, "CD14": 0.01984555177226498, "CD163": 0.026490620436259955,
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
    "Vim": 0.06671803387741954
}

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

# generate stitching/annotation function, used by panel validation and acquisition order tiling
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
    annotation_labels: List[Any] = list(image_data.coords[image_data.dims[0]].values)

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


# panel validation helpers
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
    fig.savefig(
        pathlib.Path(save_dir) / f"functional_marker_thresholds.png",
        dpi=300
    )


# acquisition order helpers
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
    post_norm_subdir: str = "", padding: int = 25, font_size: int = 100, step: int = 1
):
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
    pre_norm_tiled: Image = stitch_and_annotate_padded_img(
        pre_norm_data, padding, font_size, step=step
    )
    post_norm_tiled: Image = stitch_and_annotate_padded_img(
        post_norm_data, padding, font_size, step=step
    )

    pre_norm_tiled.save(pre_norm_stitched_path)
    post_norm_tiled.save(post_norm_stitched_path)


def run_functional_marker_positivity_tuning_tests(
    cell_table: pd.DataFrame, save_dir: Union[str, pathlib.Path],
    marker_info: Dict[str, MarkerDict], threshold_mults: List[float]
):
    """At different multipliers of each functional marker threshold, visualize across all how they
    change using a line plot

    Args:
        cell_table (pd.DataFrame):
            Cell table with clustered cell populations
        save_dir (Union[str, pathlib.Path]):
            The directory to save the line plots showing changes in number of positive cells
            per functional marker
        marker_info (str):
            For each marker, define the populations, threshold, x-range, and x-tick locations
            NOTE: used as convenience with `functional_marker_thresholding`, only threshold value
            needed
        threshold_mults (List[str]):
            The multipliers to use and visualize across all markers
    """
    marker_threshold_data = {}

    for marker in marker_info:
        marker_threshold_data[marker] = {}

        for threshold in threshold_mults:
            multiplied_threshold = marker_info[marker]["threshold"] * threshold
            marker_threshold_data[marker][threshold] = {
                "multiplied_threshold": multiplied_threshold,
                "num_positive_cells": np.sum(cell_table[marker].values >= multiplied_threshold),
                "num_positive_cells_norm": np.sum(
                    cell_table[marker].values >= multiplied_threshold
                ) / np.sum(
                    cell_table[marker].values >= marker_info[marker]["threshold"]
                )
            }

    # # plot the num positive cells normalized by 1x per marker
    # fig = plt.figure()
    # ax = fig.add_subplot(1, 1, 1)
    # threshold_mult_strs = [str(np.round(np.log2(tm), 3)) for tm in threshold_mults]

    # for i, marker in enumerate(marker_threshold_data):
    #     mult_data = [
    #         mtd["num_positive_cells_norm"] for mtd in marker_threshold_data[marker].values()
    #     ]
    #     _ = ax.plot(threshold_mult_strs, mult_data, color="gray", label=marker)

    # _ = ax.set_title(
    #     f"Positive cells per threshold, normalized by 1x",
    #     fontsize=11
    # )
    # _ = ax.yaxis.set_major_formatter(ScalarFormatter(useOffset=False))
    # _ = ax.yaxis.get_major_formatter().set_scientific(False)
    # _ = ax.set_xlabel("log2(threshold multiplier)", fontsize=7)
    # _ = ax.set_ylabel("Positive cell counts, normalized by 1x", fontsize=7)
    # _ = ax.tick_params(axis="both", which="major", labelsize=7)

    # # save the figure to save_dir
    # _ = fig.savefig(
    #     pathlib.Path(extraction_pipeline_tuning_dir) /
    #     f"functional_marker_threshold_experiments_norm.pdf",
    #     dpi=300
    # )

    # plot the raw num positive cells per marker
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    threshold_mult_strs = [np.round(np.log2(tm), 3) for tm in threshold_mults]
    threshold_mult_labels = [str(t) for t in threshold_mult_strs]

    for i, marker in enumerate(marker_threshold_data):
        mult_data = [
            mtd["num_positive_cells"] for mtd in marker_threshold_data[marker].values()
        ]
        _ = ax.plot(threshold_mult_strs, mult_data, color="gray", label=marker)

    _ = ax.set_title(f"Positive cells per threshold", fontsize=11)
    _ = ax.yaxis.set_major_formatter(ScalarFormatter(useOffset=False))
    _ = ax.yaxis.get_major_formatter().set_scientific(False)
    _ = ax.set_xlabel("log2(threshold multiplier)", fontsize=7)
    _ = ax.set_ylabel("Positive cell counts", fontsize=7)
    _ = ax.tick_params(axis="both", which="major", labelsize=7)

    # Set explicit x-axis ticks based on your calculated strings
    _ = ax.set_xticks(threshold_mult_strs)
    _ = ax.set_xticklabels(threshold_mult_labels, rotation=90)

    _ = plt.tight_layout()

    # save the figure to save_dir
    _ = fig.savefig(
        pathlib.Path(save_dir) / f"functional_marker_threshold_experiments.pdf",
        dpi=300
    )


def run_min_cell_feature_gen_fovs_dropped_tests(
    cluster_broad_df: pd.DataFrame, min_cell_params: List[int],
    compartments: List[str], metrics: List[str], save_dir: Union[str, pathlib.Path]
):
    """For the feature generation, visualize a strip plot of how varying min_cell thresholds
    determine which FOVs get dropped

    Args:
        cluster_broad_df (pd.DataFrame):
            The data in `"cluster_core_per_df.csv"`. While the process for functional, morph,
            diversity, and distance varies, the way the FOVs are deselected is the same across all
            since this dataset is always used.
        min_cell_params (List[int]):
            The min_cell thresholds to use
        compartments (List[str]):
            The compartments to visualize
        metrics (List[str]):
            The specific metric to facet the strip plot on (ex. `"cluster_count_broad"`)
        save_dir (Union[str, pathlib.Path]):
            The directory to save the strip plot(s)
    """
    total_fovs_dropped = {}
    for metric in metrics:
        total_fovs_dropped[metric] = {}

    for compartment in compartments:
        for min_cells in min_cell_params:
            for metric in metrics:
                total_fovs_dropped[metric][min_cells] = {}
                count_df = cluster_broad_df[cluster_broad_df.metric == metric]
                count_df = count_df[count_df.subset == compartment]
                all_fovs = count_df.fov.unique()

                for cell_type in count_df.cell_type.unique():
                    keep_df = count_df[count_df.cell_type == cell_type]
                    keep_df = keep_df[keep_df.value >= min_cells]
                    keep_fovs = keep_df.fov.unique()
                    total_fovs_dropped[metric][min_cells][cell_type] = \
                        len(all_fovs) - len(keep_fovs)

        for metric in metrics:
            # visualize a strip plot of FOVs dropped per min_cell test per cluster broad assignment
            fovs_dropped_dict = {metric: total_fovs_dropped[metric]}
            df = pd.DataFrame(
                [
                    {
                        "min_cells": min_cells,
                        "num_fovs_dropped": value,
                        "feature": metric,
                        "cell_type": cell_type
                    }
                    for feature, min_cells_dict in fovs_dropped_dict.items()
                    for min_cells, cell_types in min_cells_dict.items()
                    for cell_type, value in cell_types.items()
                ]
            )
            plot = sns.catplot(
                x="min_cells",
                y="num_fovs_dropped",
                hue="cell_type",
                data=df,
                kind="strip",
                palette="Set2",
                dodge=False
            )
            plot.fig.subplots_adjust(top=0.9)
            plot.fig.suptitle("Distribution of FOVs dropped across min_cells trials")
            plt.savefig(
                pathlib.Path(save_dir) /
                f"{compartment}_min_cells_{metric}_fovs_dropped_stripplot.pdf",
                dpi=300
            )


def run_cancer_mask_inclusion_tests(
    cell_table_clusters: pd.DataFrame, channel_dir: pathlib.Path, seg_dir: pathlib.Path,
    threshold_mults: List[float], save_dir: Union[str, pathlib.Path], base_sigma: int = 10,
    base_channel_thresh: float = 0.0015, base_min_mask_size: int = 7000,
    base_max_hole_size: int = 1000, base_border_size: int = 50
):
    """Create box plots showing how much of the intermediate cancer mask in `create_cancer_boundary`
    after Gaussian blurring, channel thresholding, small object removal, and hole filling remain
    at individual tuning of these params.

    Finally, create box plot showing how much of the cancer mask gets defined as cancer or stroma 
    boundary at different border_size params. Note that the default sigma, channel threhsolding, 
    small object removal, and hole filling params are used

    Args:
        cell_table_clusters (pd.DataFrame):
            The data contained in "cell_table_clusters.csv"
        channel_dir (pathlib.Path):
            The path to the "samples" folder containing all the data
        seg_dir (pathlib.Path):
            The path to the segmentation data
        threshold_mults (List[float]):
            What value to multiply the base params to test in this function
        save_dir (Union[str, pathlib.Path]):
            The directory to save the box plots
        base_sigma (int):
            The sigma value currently used for Gaussian blurring
        base_channel_thresh (float):
            The base threshold value to use after Gaussian blurring
        base_min_mask_size (int):
            The base minimum value for what objects get kept in the image
        base_max_hole_size (int):
            The base maximum size of a hole that doesn't get filled in the image
        base_border_size (int):
            The base border size to use for erosion and dilation that defines the boundaries
    """
    folders = list_folders(channel_dir)
    threshold_mult_strs = [str(np.round(np.log2(tm), 3)) for tm in threshold_mults]

    cell_boundary_sigmas = [int(tm * base_sigma) for tm in threshold_mults]
    cell_boundary_channel_threshes = [tm * base_channel_thresh for tm in threshold_mults]
    cell_boundary_min_mask_sizes = [int(tm * base_min_mask_size) for tm in threshold_mults]
    # cell_boundary_max_hole_sizes = [int(tm * max_hole_size) for tm in threshold_mults]
    cell_boundary_border_sizes = [int(tm * base_border_size) for tm in threshold_mults]

    cell_boundary_sigma_data = {s: [] for s in cell_boundary_sigmas}
    cell_boundary_channel_thresh_data = {ct: [] for ct in cell_boundary_channel_threshes}
    cell_boundary_min_mask_size_data = {mms: [] for mms in cell_boundary_min_mask_sizes}
    # cell_boundary_max_hole_size_data = {mhs: [] for mhs in cell_boundary_max_hole_sizes}
    cell_boundary_border_size_data = {bs: [] for bs in cell_boundary_border_sizes}

    i = 0
    for folder in folders:
        ecad = io.imread(os.path.join(channel_dir, folder, "ECAD.tiff"))

        # generate cancer/stroma mask by combining segmentation mask with ECAD channel
        seg_label = io.imread(os.path.join(seg_dir, folder + "_whole_cell.tiff"))[0]
        seg_mask = cancer_mask_utils.create_cell_mask(
            seg_label, cell_table_clusters, folder, ["Cancer"]
        )

        for s in cell_boundary_sigmas:
            img_smoothed = gaussian_filter(ecad, sigma=s)
            img_mask = img_smoothed > base_channel_thresh

            # clean up mask prior to analysis
            img_mask = np.logical_or(img_mask, seg_mask)
            label_mask = label(img_mask)
            label_mask = morphology.remove_small_objects(
                label_mask, min_size=base_min_mask_size
            )
            label_mask = morphology.remove_small_holes(
                label_mask, area_threshold=base_max_hole_size
            )

            percent_hit = np.sum(label_mask) / label_mask.size
            cell_boundary_sigma_data[s].append(percent_hit)

        img_smoothed = gaussian_filter(ecad, sigma=base_sigma)
        for ct in cell_boundary_channel_threshes:
            img_mask = img_smoothed > ct

            # clean up mask prior to analysis
            img_mask = np.logical_or(img_mask, seg_mask)
            label_mask = label(img_mask)
            label_mask = morphology.remove_small_objects(
                label_mask, min_size=base_min_mask_size)
            label_mask = morphology.remove_small_holes(
                label_mask, area_threshold=base_max_hole_size
            )

            percent_hit = np.sum(label_mask) / label_mask.size
            cell_boundary_channel_thresh_data[ct].append(percent_hit)

        img_smoothed = gaussian_filter(ecad, sigma=base_sigma)
        img_mask = img_smoothed > base_channel_thresh
        img_mask = np.logical_or(img_mask, seg_mask)
        label_mask_base = label(img_mask)
        for mms in cell_boundary_min_mask_sizes:
            label_mask = morphology.remove_small_objects(
                label_mask_base, min_size=mms
            )
            label_mask = morphology.remove_small_holes(
                label_mask, area_threshold=base_max_hole_size
            )

            percent_hit = np.sum(label_mask) / label_mask.size
            cell_boundary_min_mask_size_data[mms].append(percent_hit)

        # for mhs in cell_boundary_max_hole_sizes:
        #     label_mask = morphology.remove_small_objects(
        #         label_mask_base, min_size=base_max_hole_size
        #     )
        #     label_mask = morphology.remove_small_holes(
        #         label_mask, area_threshold=mhs
        #     )

        #     percent_hit = np.sum(label_mask) / label_mask.size
        #     cell_boundary_max_hole_size_data[mhs].append(percent_hit)

        label_mask = morphology.remove_small_objects(
            label_mask_base, min_size=base_min_mask_size
        )
        label_mask = morphology.remove_small_holes(
            label_mask, area_threshold=base_max_hole_size
        )
        for bs in cell_boundary_border_sizes:
            # define external borders
            external_boundary = morphology.binary_dilation(label_mask)

            for _ in range(bs):
                external_boundary = morphology.binary_dilation(external_boundary)

            external_boundary = external_boundary.astype(int) - label_mask.astype(int)
            # create interior borders
            interior_boundary = morphology.binary_erosion(label_mask)

            for _ in range(bs):
                interior_boundary = morphology.binary_erosion(interior_boundary)

            interior_boundary = label_mask.astype(int) - interior_boundary.astype(int)

            combined_mask = np.ones_like(ecad)
            combined_mask[label_mask] = 4
            combined_mask[external_boundary > 0] = 2
            combined_mask[interior_boundary > 0] = 3

            percent_border = np.sum(
                (combined_mask == 2) | (combined_mask == 3)
            ) / combined_mask.size
            cell_boundary_border_size_data[bs].append(percent_border)

        i += 1
        if i % 10 == 0:
            print(f"Processed {i} folders")

    # plot the sigma experiments
    data_sigma = []
    labels_sigma = []
    for i, (_, values_s) in enumerate(cell_boundary_sigma_data.items()):
        data_sigma.extend(values_s)
        labels_sigma.extend([threshold_mult_strs[i]] * len(values_s))

    plt.figure(figsize=(10, 6))
    sns.boxplot(x=labels_sigma, y=data_sigma)
    plt.title("Distribution of % mask included in cancer across sigma (1x = 10)")
    plt.xlabel("log2(sigma multiplier)")
    plt.ylabel("% of mask included in cancer")
    plt.savefig(
        pathlib.Path(save_dir) /
        f"sigma_cancer_mask_inclusion_box.pdf",
        dpi=300
    )

    # plot the channel thresh experiments
    data_channel_thresh = []
    labels_channel_thresh = []
    for i, (_, values_ct) in enumerate(cell_boundary_channel_thresh_data.items()):
        data_channel_thresh.extend(values_ct)
        labels_channel_thresh.extend([threshold_mult_strs[i]] * len(values_ct))

    plt.figure(figsize=(10, 6))
    sns.boxplot(x=labels_channel_thresh, y=data_channel_thresh)
    plt.title("Distribution of % mask included in cancer across smoothing thresholds (1x = 0.0015)")
    plt.xlabel("log2(smooth thresh multiplier)")
    plt.ylabel("% of mask included in cancer")
    plt.savefig(
        pathlib.Path(save_dir) /
        f"smooth_thresh_cancer_mask_inclusion_box.pdf",
        dpi=300
    )

    # plot the min mask size experiments
    data_min_mask_size = []
    labels_min_mask_size = []
    for i, (_, values_mms) in enumerate(cell_boundary_min_mask_size_data.items()):
        data_min_mask_size.extend(values_mms)
        labels_min_mask_size.extend([threshold_mult_strs[i]] * len(values_mms))

    # Creating the boxplot
    plt.figure(figsize=(10, 6))
    sns.boxplot(x=labels_min_mask_size, y=data_min_mask_size)
    plt.title("Distribution of % mask included in cancer across min mask sizes (1x = 7000)")
    plt.xlabel("log2(min mask size multiplier)")
    plt.ylabel("% of mask included in cancer")
    plt.savefig(
        pathlib.Path(save_dir) /
        f"min_mask_size_cancer_mask_inclusion_box.pdf",
        dpi=300
    )

    # # plot the max hole size experiments
    # data_max_hole_size = []
    # labels_max_hole_size = []
    # for i, (_, values_mhs) in enumerate(cell_boundary_max_hole_size_data.items()):
    #     data_max_hole_size.extend(values_mhs)
    #     labels_max_hole_size.extend([threshold_mult_strs[i]] * len(values_mhs))

    # plt.figure(figsize=(10, 6))
    # sns.boxplot(x=labels_max_hole_size, y=data_max_hole_size)
    # plt.title("Distribution of % mask included in cancer across max hole sizes (1x = 1000)")
    # plt.xlabel("log2(max hole size multiplier)")
    # plt.ylabel("% of mask included in cancer")
    # plt.savefig(
    #     pathlib.Path(save_dir) /
    #     f"max_hole_size_cancer_mask_inclusion_box.pdf",
    #     dpi=300
    # )

    # plot the border size experiments
    data_border_size = []
    labels_border_size = []
    for i, (_, values_bs) in enumerate(cell_boundary_border_size_data.items()):
        data_border_size.extend(values_bs)
        labels_border_size.extend([threshold_mult_strs[i]] * len(values_bs))

    plt.figure(figsize=(10, 6))
    sns.boxplot(x=labels_border_size, y=data_border_size)
    plt.title("Distribution of % cancer boundary across border sizes (1x = 50)")
    plt.xlabel("log2(border size multiplier)")
    plt.ylabel("% of mask identified as cancer boundary")
    plt.savefig(
        pathlib.Path(save_dir) / f"border_size_cancer_region_percentages_box.pdf",
        dpi=300
    )


def random_feature_generation(combined_df, seed_num, tonic_features_df, feature_metadata):
    """
    Shuffle the Patient_ID (and thus Clinical_benefit) in the TONIC feature dataframe to
    re-generate the comparison of features.
    """
    patients = np.unique(combined_df.Patient_ID)
    np.random.seed(seed_num)
    p_rand = patients.copy()
    np.random.shuffle(p_rand)

    outcome = dict(zip(combined_df.Patient_ID, combined_df.Clinical_benefit))
    patient_map = pd.DataFrame({'Patient_ID': patients, 'Patients_rand': p_rand})

    combined_df_alt = combined_df.merge(patient_map, on=['Patient_ID'])
    outcome_map = pd.DataFrame({'Patients_rand': outcome.keys(), 'Clinical_benefit_rand': outcome.values()})

    combined_df_alt2 = combined_df_alt.merge(outcome_map, on=['Patients_rand'])
    combined_df_alt2 = combined_df_alt2.rename(columns={'Patient_ID': 'og_Patient_ID', 'Patients_rand': 'Patient_ID',
                                                        'Clinical_benefit': 'og_Clinical_benefit',
                                                        'Clinical_benefit_rand': 'Clinical_benefit'})

    # settings for generating hits
    method = 'ttest'
    total_dfs = []

    for comparison in combined_df_alt2.Timepoint.unique():
        population_df = compare_populations(feature_df=combined_df_alt2, pop_col='Clinical_benefit',
                                            timepoints=[comparison], pop_1='No', pop_2='Yes', method=method)

        if np.sum(~population_df.log_pval.isna()) == 0:
            continue
        long_df = population_df[['feature_name_unique', 'log_pval', 'mean_diff', 'med_diff']]
        long_df['comparison'] = comparison
        long_df = long_df.dropna()
        long_df['pval'] = 10 ** (-long_df.log_pval)
        long_df['fdr_pval'] = multipletests(long_df.pval, method='fdr_bh')[1]
        total_dfs.append(long_df)

    # summarize hits from all comparisons
    ranked_features_df = pd.concat(total_dfs)
    ranked_features_df['log10_qval'] = -np.log10(ranked_features_df.fdr_pval)

    # create importance score
    # get ranking of each row by log_pval
    ranked_features_df['pval_rank'] = ranked_features_df.log_pval.rank(ascending=False)
    ranked_features_df['cor_rank'] = ranked_features_df.med_diff.abs().rank(ascending=False)
    ranked_features_df['combined_rank'] = (ranked_features_df.pval_rank.values + ranked_features_df.cor_rank.values) / 2

    # generate importance score
    max_rank = len(~ranked_features_df.med_diff.isna())
    normalized_rank = ranked_features_df.combined_rank / max_rank
    ranked_features_df['importance_score'] = 1 - normalized_rank

    ranked_features_df = ranked_features_df.sort_values('importance_score', ascending=False)
    # ranked_features_df = ranked_features_df.sort_values('fdr_pval', ascending=True)

    # generate signed version of score
    ranked_features_df['signed_importance_score'] = ranked_features_df.importance_score * np.sign(
        ranked_features_df.med_diff)

    # add feature type
    ranked_features_df = ranked_features_df.merge(feature_metadata, on='feature_name_unique', how='left')

    feature_type_dict = {'functional_marker': 'phenotype', 'linear_distance': 'interactions',
                         'density': 'density', 'cell_diversity': 'diversity', 'density_ratio': 'density',
                         'mixing_score': 'interactions', 'region_diversity': 'diversity',
                         'compartment_area_ratio': 'compartment', 'density_proportion': 'density',
                         'morphology': 'phenotype', 'pixie_ecm': 'ecm', 'fiber': 'ecm', 'ecm_cluster': 'ecm',
                         'compartment_area': 'compartment', 'ecm_fraction': 'ecm'}
    ranked_features_df['feature_type_broad'] = ranked_features_df.feature_type.map(feature_type_dict)

    # identify top features
    ranked_features_df['top_feature'] = False
    ranked_features_df.iloc[:100, -1] = True

    ranked_features_df['full_feature'] = ranked_features_df.feature_name_unique + '-' + ranked_features_df.comparison + '_' + ranked_features_df.compartment
    tonic_features_df['full_feature'] = tonic_features_df.feature_name_unique + '-' + tonic_features_df.comparison + '_' + tonic_features_df.compartment

    random_features = ranked_features_df[:len(tonic_features_df)].full_feature
    tonic_features = tonic_features_df.full_feature

    # compare against TONIC features
    intersection_of_features = set(tonic_features).intersection(set(random_features))
    union_of_features = set(tonic_features).union(set(random_features))
    jaccard_score = len(intersection_of_features) / len(union_of_features)

    return intersection_of_features, jaccard_score, ranked_features_df[:len(tonic_features_df)]

class MembraneMarkersSegmentationPlot:
    def __init__(
        self,
        fov: str,
        image_data: PathLike,
        segmentation_dir: PathLike,
        membrane_channels: List[str],
        overlay_channels: str | List[str],
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
            segmentation_dir=self.seg_dir / "deepcell_output",
            data_dir=self.seg_dir / "deepcell_input",
            img_overlay_chans=self.overlay_channels,
            seg_overlay_comp="whole_cell",
        )

        self.fov_cell_segmentation = load_imgs_from_dir(
            data_dir=self.seg_dir / "deepcell_output",
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
        overlay_channels: str | List[str] = ["nuclear_channel", "membrane_channel"],
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
            segmentation_dir=self.seg_dir / "deepcell_output",
            data_dir=self.seg_dir / "deepcell_input",
            img_overlay_chans=self.overlay_channels,
            seg_overlay_comp="whole_cell",
        )
        self.fov_cell_segmentation = load_imgs_from_dir(
            data_dir=self.seg_dir / "deepcell_output",
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
            segmentation_dir=self.seg_dir / "deepcell_output",
            data_dir=self.seg_dir / "deepcell_input",
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


def calculate_feature_corr(timepoint_features,
                            top_features,
                            remaining_features,
                            top: bool = True,
                            n_iterations: int = 1000):
    """Compares the correlation between
            1. response-associated features to response-associated features
            2. response-associated features to remaining features
        by randomly sampling features with replacement.

    Parameters
    timepoint_features: pd.DataFrame
        dataframe containing the feature values for every patient (feature_name_unique, normalized_mean, Patient_ID, Timepoint)
    top_features: pd.DataFrame
        dataframe containing the top response-associated features (feature_name_unique, Timepoint)
    remaining features: pd.DataFrame
        dataframe containing non response-associated features (feature_name_unique, Timepoint)
    top: bool (default = True)
        boolean indicating if the comparison 1. (True) or 2. (False)
    n_iterations: int (default = 1000)
        number of features randomly selected for comparison
    ----------
    Returns
    corr_arr: np.array
        array containing the feature correlation values
    ----------
    """
    corr_arr = []
    for _ in range(n_iterations):
        #select feature 1 as a random feature from the top response-associated feature list
        rand_sample1 = top_features.sample(n = 1)
        f1 = timepoint_features.iloc[np.where((timepoint_features['feature_name_unique'] == rand_sample1['feature_name_unique'].values[0]) & (timepoint_features['Timepoint'] == rand_sample1['Timepoint'].values[0]))[0], :]
        if top == True:
            #select feature 2 as a random feature from the top response-associated list, ensuring f1 != f2
            rand_sample2 = rand_sample1
            while (rand_sample2.values == rand_sample1.values).all():
                rand_sample2 = top_features.sample(n = 1)
        else:
            #select feature 2 as a random feature from the remaining feature list
            rand_sample2 = remaining_features.sample(n = 1)

        f2 = timepoint_features.iloc[np.where((timepoint_features['feature_name_unique'] == rand_sample2['feature_name_unique'].values[0]) & (timepoint_features['Timepoint'] == rand_sample2['Timepoint'].values[0]))[0], :]
        merged_features = pd.merge(f1, f2, on = 'Patient_ID') #finds Patient IDs that are shared across timepoints to compute correlation
        corrval = np.abs(merged_features['normalized_mean_x'].corr(merged_features['normalized_mean_y'])) #regardless of direction
        corr_arr.append(corrval)

    return np.array(corr_arr)
