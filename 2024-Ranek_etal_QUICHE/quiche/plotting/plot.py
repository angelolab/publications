import os
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
import seaborn as sns
import scanpy as sc
import quiche as qu
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgba, Normalize
import matplotlib.colors as colors
from matplotlib.cm import ScalarMappable
from matplotlib import colormaps
from matplotlib import gridspec, patches, cm
from matplotlib.figure import Figure
from matplotlib.colorbar import ColorbarBase
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Wedge, FancyArrowPatch
import matplotlib
import matplotlib.patheffects as patheffects
from matplotlib.patches import Rectangle
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.cm as cm
import pathlib
import shutil
from typing import Dict, List, Literal, Optional, Tuple, Union
import natsort
from pandas.core.groupby.generic import DataFrameGroupBy
import skimage
from matplotlib.axes import Axes
import xarray as xr
from alpineer import image_utils, io_utils, load_utils, misc_utils
from skimage import io
from skimage.segmentation import find_boundaries
from skimage.util import img_as_ubyte
import networkx as nx
from tqdm.auto import tqdm
from collections import Counter
import anndata
from ark import settings
from ark.utils.data_utils import ClusterMaskData, erode_mask, generate_cluster_mask, save_fov_mask, map_segmentation_labels
matplotlib.use('Agg')
plt.rcParams['agg.path.chunksize'] = 10000
plt.rcParams['path.simplify'] = True
plt.rcParams['path.simplify_threshold'] = 0.0

def generate_colors(cmap="viridis", n_colors=3, alpha=.4):
    """Generate colors from matplotlib colormap; pass list to use exact colors"""
    if not isinstance(n_colors, int) or (n_colors < 2) or (n_colors > 6):
        raise ValueError("n_colors must be an integer between 2 and 6")
    if isinstance(cmap, list):
        colors = [to_rgba(color, alpha=alpha) for color in cmap]
    else:
        scalar_mappable = ScalarMappable(cmap=cmap)
        colors = scalar_mappable.to_rgba(range(n_colors), alpha=alpha).tolist()
    return colors[:n_colors]

def plot_grid_enrichment(df, hue_order, colors_dict, selected_grid_list, num_grids_x, num_grids_y, save_directory, filename_save):
    _, axes = plt.subplots(nrows=1, ncols=1, figsize=(6,6))
    sns.scatterplot(x = 'x', y = 'y', hue = 'group', data =df, alpha=0.5, palette=colors_dict, hue_order = hue_order, ax = axes)
    axes.tick_params(labelsize=10)
    plt.xlabel('Y', fontsize = 12)
    plt.ylabel('X', fontsize = 12)
    plt.xlim(0, 1)
    plt.ylim(0,1)
    plt.legend(loc='upper right', bbox_to_anchor=(1.15, 1.0))
    for selected_grid in selected_grid_list:
        selected_grid_rect = Rectangle(
            ((selected_grid[0] - 1) * (1 / num_grids_x), (selected_grid[1] - 1) * (1 / num_grids_y)),
            (1 / num_grids_x), (1 / num_grids_y),
            edgecolor='k', linewidth=1.5, fill=False
        )
        plt.gca().add_patch(selected_grid_rect)
    plt.savefig(os.path.join(save_directory, filename_save+'.pdf'), bbox_inches = 'tight')

def plot_purity(purity_score,
                figsize = (4,4),
                colors_dict = None,
                save_directory = os.path.join('figures', 'simulated'),
                filename_save = 'simulated',
                annot_list = 'all'):

    _, axes = plt.subplots(1, 1, figsize = figsize, gridspec_kw={'hspace': 0.45, 'wspace': 0.4, 'bottom':0.15})

    if annot_list is 'all':
        color_list = [colors_dict[key] for key in list(purity_score.index)] 
        g = purity_score.transpose().plot(kind='bar', stacked=True, ax = axes, color = color_list)
    
    else:
        color_list = [colors_dict[key] for key in purity_score.transpose().loc[annot_list, :].columns]
        g = purity_score.transpose().loc[annot_list, :].plot(kind='bar', stacked=True, ax = axes, color = color_list)

    g.tick_params(labelsize=12)
    g.set_ylabel('purity score', fontsize = 12)
    g.set_xlabel('')
    g.legend(bbox_to_anchor=(1.05,1.05), prop={'size':10})
    plt.savefig(os.path.join(save_directory, filename_save + '.pdf'), bbox_inches = "tight")

def plot_cell_types_per_group(proportion_df,
                            figsize = (6,4),
                            fov_key = 'Patient_ID',
                            annot_key = 'kmeans_cluster_10',
                            cluster_id = '8',
                            labels_key = 'mask_name',
                            colors_dict = None,
                            save_directory = os.path.join('figures', 'simulated'),
                            filename_save = 'simulated',
                            ylabel = 'cell_counts',
                            condition_list = ['cancer_core', 'cancer_border']):

    subset_df = proportion_df[proportion_df[annot_key] == cluster_id]
    
    ordered_ids = sorted(subset_df[fov_key].unique(), 
                        key=lambda x: (not x.endswith(condition_list[0]), not x.endswith(condition_list[1])))

    pivot_df = subset_df.pivot(index=fov_key, columns=labels_key, values='cell_count').loc[ordered_ids]

    color_list = [colors_dict[key] for key in list(pivot_df.columns)] 
    _, axes = plt.subplots(1, 1, figsize = figsize, gridspec_kw={'hspace': 0.45, 'wspace': 0.4, 'bottom':0.15})
    g = pivot_df.plot(kind='bar', stacked=True, ax = axes, color = color_list)
    g.tick_params(labelsize=12)
    g.set_ylabel(ylabel, fontsize = 12)
    g.set_xlabel(fov_key, fontsize = 12)
    g.set_title(f'cluster {cluster_id}', fontsize = 12)
    g.legend(bbox_to_anchor=(1.0,1.0), prop={'size':10})
    plt.savefig(os.path.join(save_directory, filename_save + '.pdf'), bbox_inches = "tight")

def plot_masks_group(adata_niche,
            condition_key = 'condition',
            save_directory = os.path.join('figures', 'simulated', 'kmeans', 'overlays'),
            condition_list = ['cancer_core', 'cancer_border'],
            annot_key = 'kmeans_cluster_10',
            fov_key = 'fov',
            seg_dir = r'/Volumes/Shared/Noah Greenwald/TNBC_Cohorts/SPAIN/segmentation/samples/deepcell_output',
            colors_dict = None):

    df = adata_niche.to_df()
    df[list(adata_niche.obs.columns)] = adata_niche.obs
    df[annot_key] = pd.Categorical(df[annot_key])

    if colors_dict is None:
        colors = sns.color_palette("tab20")
        colors_dict = dict(zip(adata_niche.obs[annot_key].cat.categories, [colors[i] for i in range(0, len(adata_niche.obs[annot_key].cat.categories))]))
        colors_dict[np.nan] = 'black'

    colormap = pd.DataFrame({annot_key: list(colors_dict.keys()),
                            'color': list(colors_dict.values())})

    for condition in condition_list:

        df_cond = df[df[condition_key] == condition]
        df_cond.loc[:, annot_key] = pd.Categorical(df_cond.loc[:, annot_key])

        # Directory to save overlays
        overlay_out_dir = os.path.join(save_directory, str(condition))
        if not os.path.exists(overlay_out_dir):
            os.makedirs(overlay_out_dir)


        qu.tl.cohort_cluster_plot(
            fovs=list(df_cond.fov.unique()),
            seg_dir=seg_dir,
            save_dir=overlay_out_dir,
            cell_data=df_cond,
            erode=True,
            fov_col= fov_key,
            label_col='label',
            cluster_col=annot_key,
            seg_suffix="_whole_cell.tiff",
            cmap=colormap,
            display_fig=True,
        )
    
def set_minimum_color_for_colormap(cmap, default=(0, 0, 0, 1)):
    """ Changes minimum value in provided colormap to black (#000000) or provided color

    This is useful for instances where zero-valued regions of an image should be
    distinct from positive regions (i.e transparent or non-colormap member color)

    Args:
        cmap (matplotlib.colors.Colormap):
            matplotlib color map
        default (Iterable):
            RGBA color values for minimum color. Default is black, (0, 0, 0, 1).

    Returns:
        matplotlib.colors.Colormap:
            corrected colormap
    """
    cmapN = cmap.N
    corrected = cmap(np.arange(cmapN))
    corrected[0, :] = list(default)
    return colors.ListedColormap(corrected)

def save_colored_mask(
    fov: str,
    save_dir: str,
    suffix: str,
    data: np.ndarray,
    cmap: colors.ListedColormap,
    norm: colors.BoundaryNorm,
) -> None:
    """Saves the colored mask to the provided save directory.

    Args:
        fov (str):
            The name of the FOV.
        save_dir (str):
            The directory where the colored mask will be saved.
        suffix (str):
            The suffix to append to the FOV name.
        data (np.ndarray):
            The mask to save.
        cmap (colors.ListedColormap):
            The colormap to use for the mask.
        norm (colors.BoundaryNorm):
            The normalization to use for the mask.
    """

    # Create the save directory if it does not exist
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)

    # Create the colored mask
    colored_mask = img_as_ubyte(cmap(norm(data)))

    # Save the image
    image_utils.save_image(
        fname=os.path.join(save_dir, f"{fov}{suffix}"),
        data=colored_mask,
    )

def plot_continuous_variable(
    image: np.ndarray,
    name: str,
    stat_name: str,
    cmap: Union[colors.Colormap, str],
    norm: colors.Normalize = None,
    cbar_visible: bool = True,
    dpi: int = 300,
    figsize: tuple[int, int] = (10, 10),
) -> Figure:
    """

    Plots an image measuring some type of continuous variable with a user provided colormap.

    Args:
        image (np.ndarray):
            An array representing an image to plot.
        name (str):
            The name of the image.
        stat_name (str):
            The name of the statistic to plot, this will be the colormap's label.
        cmap (colors.Colormap, str, optional): A colormap to plot the array with.
            Defaults to "viridis".
        cbar_visible (bool, optional): A flag for setting the colorbar on or not.
            Defaults to True.
        norm (colors.Normalize, optional): A normalization to apply to the colormap.
        dpi (int, optional):
            The resolution of the image. Defaults to 300.
        figsize (tuple[int, int], optional):
            The size of the image. Defaults to (10, 10).

    Returns:
        Figure : The Figure object of the image.
    """
    fig: Figure = plt.figure(figsize=figsize, dpi=dpi)
    fig.set_layout_engine(layout="tight")
    gs = gridspec.GridSpec(nrows=1, ncols=1, figure=fig)
    fig.suptitle(f"{name}")

    # Image axis
    ax: Axes = fig.add_subplot(gs[0, 0])
    ax.axis("off")
    ax.grid(visible=False)

    im = ax.imshow(
        X=image,
        cmap=cmap,
        norm=norm,
        origin="upper",
        aspect="equal",
        interpolation="none",
    )

    if cbar_visible:
        # Manually set the colorbar
        divider = make_axes_locatable(fig.gca())
        cax = divider.append_axes(position="right", size="5%", pad="3%")

        fig.colorbar(mappable=im, cax=cax, orientation="vertical",
                     use_gridspec=True, pad=0.1, shrink=0.9, drawedges=False, label=stat_name)

    return fig

def color_segmentation_by_stat(
    fovs: List[str],
    data_table: pd.DataFrame,
    seg_dir: Union[pathlib.Path, str],
    save_dir: Union[pathlib.Path, str],
    fov_col: str = settings.FOV_ID,
    label_col: str = settings.CELL_LABEL,
    stat_name: str = settings.CELL_TYPE,
    cmap: str = "viridis",
    reverse: bool = False,
    seg_suffix: str = "_whole_cell.tiff",
    cbar_visible: bool = True,
    style: str = "seaborn-v0_8-paper",
    erode: bool = False,
    display_fig: bool = False,
    fig_file_type: str = "png",
    vmin = None, 
    vmax = None,
    figsize: tuple = (10, 10),
    dpi: int = 300,
):
    """
    Colors segmentation masks by a given continuous statistic.

    Args:
        fovs: (List[str]):
            A list of FOVs to plot.
        data_table (pd.DataFrame):
            A DataFrame containing FOV and segmentation label identifiers
            as well as a collection of statistics for each label in a segmentation
            mask such as:

                - `fov_id` (identifier)
                - `label` (identifier)
                - `area` (statistic)
                - `fiber` (statistic)
                - etc...

        seg_dir (Union[pathlib.Path, str]):
            Path to the directory containing segmentation masks.
        save_dir (Union[pathlib.Path, str]):
            Path to the directory where the colored segmentation masks will be saved.
        fov_col: (str, optional):
            The name of the column in `data_table` containing the FOV identifiers.
            Defaults to "fov".
        label_col (str, optional):
            The name of the column in `data_table` containing the segmentation label identifiers.
            Defaults to "label".
        stat_name (str):
            The name of the statistic to color the segmentation masks by. This should be a column
            in `data_table`.
        seg_suffix (str, optional):
            The suffix of the segmentation file and it's file extension. Defaults to
            "_whole_cell.tiff".
        cmap (str, optional): The colormap for plotting. Defaults to "viridis".
        reverse (bool, optional):
            A flag to reverse the colormap provided. Defaults to False.
        cbar_visible (bool, optional):
            A flag to display the colorbar. Defaults to True.
        erode (bool, optional): Option to "thicken" the cell boundary via the segmentation label
            for visualization purposes. Defaults to False.
        style (str, optional): Set the matplotlib style image style. Defaults to 
            "seaborn-v0_8-paper".
            View the available styles here: 
            https://matplotlib.org/stable/gallery/style_sheets/style_sheets_reference.html
            Or run matplotlib.pyplot.style.available in a notebook to view all the styles.
        display_fig: (bool, optional):
            Option to display the cluster mask plots as they are generated. Defaults to False.
        fig_file_type (str, optional): The file type to save figures as. Defaults to 'png'.
        figsize (tuple, optional):
            The size of the figure to display. Defaults to (10, 10).
        dpi (int, optional):
            The resolution of the image to use for saving. Defaults to 300.
    """
    plt.style.use(style)

    if not isinstance(seg_dir, pathlib.Path):
        seg_dir = pathlib.Path(seg_dir)

    if not isinstance(save_dir, pathlib.Path):
        save_dir = pathlib.Path(save_dir)

    io_utils.validate_paths([seg_dir])

    try:
        io_utils.validate_paths([save_dir])
    except FileNotFoundError:
        save_dir.mkdir(parents=True, exist_ok=True)

    misc_utils.verify_in_list(
        statistic_name=[fov_col, label_col, stat_name],
        data_table_columns=data_table.columns,
    )

    if not (save_dir / "continuous_plots").exists():
        (save_dir / "continuous_plots").mkdir(parents=True, exist_ok=True)
    if not (save_dir / "colored").exists():
        (save_dir / "colored").mkdir(parents=True, exist_ok=True)

    # filter the data table to only include the FOVs we want to plot
    data_table = data_table[data_table[fov_col].isin(fovs)]

    data_table_subset_groups: DataFrameGroupBy = (
        data_table[[fov_col, label_col, stat_name]]
        .sort_values(by=[fov_col, label_col], key=natsort.natsort_keygen())
        .groupby(by=fov_col)
    )

    norm = colors.Normalize(vmin=vmin, vmax=vmax)

    if reverse:
        # Adding the suffix "_r" will reverse the colormap
        cmap = f"{cmap}_r"

    # Prepend black to the colormap
    color_map = set_minimum_color_for_colormap(
        cmap=colormaps[cmap], default=(0, 0, 0, 1)
    )

    with tqdm(
        total=len(data_table_subset_groups),
        desc=f"Generating {stat_name} Plots",
        unit="FOVs",
    ) as pbar:
        for fov, fov_group in data_table_subset_groups:
            pbar.set_postfix(FOV=fov)

            label_map: np.ndarray = io.imread(seg_dir / f"{fov}{seg_suffix}")

            if erode:
                label_map = erode_mask(
                    label_map, connectivity=2, mode="thick", background=0
                )

            mapped_seg_image: np.ndarray = map_segmentation_labels(
                labels=fov_group[label_col],
                values=fov_group[stat_name],
                label_map=label_map,
            )

            fig = plot_continuous_variable(
                image=mapped_seg_image,
                name=fov,
                stat_name=stat_name,
                norm=norm,
                cmap=color_map,
                cbar_visible=cbar_visible,
                figsize=figsize,
                dpi=dpi,
            )
            fig.savefig(fname=os.path.join(save_dir, "continuous_plots", f"{fov}.{fig_file_type}"))

            save_colored_mask(
                fov=fov,
                save_dir=save_dir / "colored",
                suffix=".tiff",
                data=mapped_seg_image,
                cmap=color_map,
                norm=norm,
            )
            if display_fig:
                fig.show(warn=False)
            else:
                plt.close(fig)

            pbar.update(1)

def plot_stacked_covariate(mdata,
                            niches = None,
                            feature_key = 'quiche',
                            annot_key = 'quiche_niche',
                            patient_key = 'Patient_ID',
                            design_key = 'Relapse',
                            patient_niche_threshold = 5,
                            figsize = (3, 12),
                            color = ['#377eb8', '#e41a1c'],
                            save_directory = 'figures',
                            filename_save = None):

    cov_count_df = qu.tl.compute_patient_proportion(mdata, niches = niches, feature_key = feature_key, annot_key = annot_key, patient_key = patient_key, design_key = design_key, patient_niche_threshold = patient_niche_threshold)
    #Sort the dataframe
    _, axes = plt.subplots(1, 1, figsize = figsize, gridspec_kw={'hspace': 0.45, 'wspace': 0.4, 'bottom':0.15})

    #Pivot the dataframe for stacked barplot
    df_pivot = cov_count_df.pivot_table(index=annot_key, columns=design_key, values='prop_cov')
    df_pivot = df_pivot.loc[cov_count_df[[annot_key, 'mean_logFC']].sort_values(by = 'mean_logFC', ascending = False).drop_duplicates()[annot_key].values]

    df_pivot.plot(kind='barh', stacked=True, ax = axes, color = color)
    plt.ylabel('prop_cov')
    if filename_save is not None:   
        plt.savefig(os.path.join(save_directory, filename_save + '.pdf'))

def plot_masks_stat(adata_niche,
            condition_key = 'condition',
            save_directory = os.path.join('figures', 'simulated', 'quiche', 'overlays'),
            condition_list = ['cancer_core', 'cancer_border'],
            fov_key = 'fov',
            cmap = 'greys_r',
            vmin = None,
            vmax = None,
            figsize = (4,4),
            stat_key = 'logSpatialFDR',
            seg_dir = r'/Volumes/Shared/Noah Greenwald/TNBC_Cohorts/SPAIN/segmentation/samples/deepcell_output'):

    df = adata_niche.to_df()
    df[list(adata_niche.obs.columns)] = adata_niche.obs

    for condition in condition_list:

        df_cond = df[df[condition_key] == condition].copy()
        df_cond['fov'] = list(df_cond['fov'].values)

        # Directory to save overlays
        overlay_out_dir = os.path.join(save_directory, str(condition), str(stat_key))
        if not os.path.exists(overlay_out_dir):
            os.makedirs(overlay_out_dir)

        qu.pl.color_segmentation_by_stat(fovs = list(np.unique(df_cond[fov_key])),
            data_table = df_cond,
            seg_dir = seg_dir,
            save_dir =overlay_out_dir,
            fov_col = fov_key,
            label_col = 'label',
            stat_name = stat_key,
            cmap = cmap,
            seg_suffix = "_whole_cell.tiff",
            cbar_visible = True,
            style = "seaborn-v0_8-paper",
            erode = True,
            vmin = vmin, 
            vmax = vmax,
            display_fig = True,
            fig_file_type = "pdf",
            figsize = figsize,
            dpi = 600)

def beeswarm_prev(mdata,
            feature_key = "quiche",
            alpha: float = 0.05,
            niches = None,
            figsize = (6, 12),
            annot_key = 'quiche_niche',
            design_key = 'condition',
            patient_key = 'sample',
            xlim = None,
            fontsize = 8,
            xlim_prev =[-1, 1],
            colors_dict = {'0':'#377eb8', '1':'#e41a1c'},
            save_directory = 'figures',
            filename_save = None):
        """Plot beeswarm plot of logFC against niche neighborhoods. Modified from pertpy: https://pertpy.readthedocs.io/en/latest/usage/tools/pertpy.tools.Milo.html
        """
        sns.set_style('ticks')
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
            raise RuntimeError(
                f"{annot_key} not defined"
            ) from None

        try:
            nhood_adata.obs["logFC"]
        except KeyError:
            raise RuntimeError(
                "Unable to find 'logFC' in mdata.uns['nhood_adata'].obs. Run 'core.da_nhoods(adata)' first."
            ) from None

        sorted_annos = (
            nhood_adata.obs[[annot_key, "logFC"]].groupby(annot_key).mean().sort_values("logFC", ascending=True).index
        )

        anno_df = nhood_adata.obs[[annot_key, "logFC", "SpatialFDR"]].copy()
        anno_df["is_signif"] = anno_df["SpatialFDR"] < alpha
        anno_df = anno_df[anno_df[annot_key] != "nan"]

        cmap_df = pd.DataFrame(mdata[feature_key].var.groupby(annot_key)['logFC'].mean(), columns = ['logFC'])
        cmap = np.full(np.shape(mdata[feature_key].var.groupby(annot_key)['logFC'].mean())[0], 'lightgrey', dtype = 'object')
        cmap[mdata[feature_key].var.groupby(annot_key)['logFC'].mean() <= 0] = list(colors_dict.values())[0]
        cmap[mdata[feature_key].var.groupby(annot_key)['logFC'].mean() > 0] = list(colors_dict.values())[1]
        cmap_df['cmap'] = cmap
        fig = plt.figure(figsize=figsize)
        gs = GridSpec(1, 2, width_ratios=[1, 0.4])  # 2 columns with equal width

        # Plot the strip plot in the first column

        ax0 = plt.subplot(gs[0])
        ax1 = plt.subplot(gs[1])
        
        g = sns.violinplot(
                data=anno_df,
                y=annot_key,
                x="logFC",
                order=sorted_annos,
                inner=None,
                orient="h",
                palette= cmap_df.loc[sorted_annos]['cmap'].values,
                linewidth=0,
                scale="width",
                alpha = 0.8,
                ax=ax0
            )

        g = sns.stripplot(
            data=anno_df,
            y=annot_key,
            x="logFC",
            order=sorted_annos,
            hue="is_signif",
            palette=["grey", "black"],
            size=2,
            orient="h",
            alpha=0.5,
            ax = ax0
        )

        ax0.set_xlabel('log2(fold change)', fontsize = fontsize)
        ax0.set_ylabel(annot_key +' neighborhoods', fontsize = fontsize)

        cov_count_df = qu.tl.compute_patient_proportion(mdata, niches=niches, feature_key=feature_key, annot_key=annot_key, patient_key=patient_key, design_key=design_key, patient_niche_threshold=3)

        # Pivot the dataframe for the stacked barplot
        df_pivot = cov_count_df.pivot_table(index=annot_key, columns=design_key, values='prop_cov')
        df_pivot = df_pivot.loc[cov_count_df[[annot_key, 'med_logFC']].sort_values(by='med_logFC', ascending=False).drop_duplicates()[annot_key].values]
        
        min_name = mdata[feature_key].var.groupby(design_key)['logFC'].mean().idxmin()
        max_name = mdata[feature_key].var.groupby(design_key)['logFC'].mean().idxmax()
        df_pivot.loc[:, min_name] = df_pivot.loc[:, min_name] * -1
        df_pivot = df_pivot.loc[sorted_annos, :]

        #works for binary alone
        ax1.barh(df_pivot.index, df_pivot.loc[:, min_name],  height = 0.5,color=colors_dict[min_name], edgecolor='none', label=df_pivot.iloc[:, 0].name)
        ax1.barh(df_pivot.index, df_pivot.loc[:, max_name], height = 0.5, color=colors_dict[max_name], edgecolor='none', left=0, label=df_pivot.iloc[:, 1].name)

        ax1.set_xlim(xlim_prev[0], xlim_prev[1])

        # Align y-ticks
        ax1.set_yticks(ax0.get_yticks())
        ax1.set_yticklabels(ax0.get_yticklabels())
        ax1.set_ylim(ax0.get_ylim())
        ax1.set_yticklabels([])
        ax1.set_yticks([])

        ax1.set_xlabel('proportion of patients', fontsize=fontsize)
        ax1.set_ylabel('')
        ax1.legend(loc="upper left", title=design_key, frameon=False, bbox_to_anchor=(1.05, 1), prop={'size': fontsize}, markerscale=1, fontsize=fontsize)
        ax0.tick_params(labelsize=fontsize)
        ax1.tick_params(labelsize=fontsize)

        # # Change x-tick labels to be positive
        xticks = ax1.get_xticks()
        ax1.set_xticklabels([str(abs(np.round(x, 2))) for x in xticks])

        if xlim is not None:
            ax0.set_xlim(xlim[0], xlim[1])

        ax0.axvline(x=0, ymin=0, ymax=1, color="black", linestyle="--", linewidth = 1)
        ax1.axvline(x=0, ymin=0, ymax=1, color="black", linestyle="--", linewidth = 1)
        ax0.legend(loc="upper right", title=f"< {int(alpha * 100)}% spatial FDR", bbox_to_anchor=(1, 1), frameon=False, prop={'size':fontsize}, markerscale=1, fontsize = fontsize)

        ax0.set_title(f'{min_name} vs {max_name}', fontsize = fontsize)
        if filename_save is not None:
            plt.savefig(os.path.join(save_directory, filename_save+'.pdf'), bbox_inches = 'tight')

def beeswarm(mdata,
            feature_key = "quiche",
            alpha: float = 0.05,
            niches = None,
            figsize = (6, 12),
            annot_key = 'quiche_niche',
            design_key = 'condition',
            patient_key = 'sample',
            xlim = None,
            save_directory = 'figures',
            ax = None,
            label = True,
            filename_save = None):
        """Plot beeswarm plot of logFC against niche neighborhoods. Modified from pertpy: https://pertpy.readthedocs.io/en/latest/usage/tools/pertpy.tools.Milo.html
        """
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
            raise RuntimeError(
                f"{annot_key} not defined"
            ) from None

        try:
            nhood_adata.obs["logFC"]
        except KeyError:
            raise RuntimeError(
                "Unable to find 'logFC' in mdata.uns['nhood_adata'].obs. Run 'core.da_nhoods(adata)' first."
            ) from None

        sorted_annos = (
            nhood_adata.obs[[annot_key, "logFC"]].groupby(annot_key).mean().sort_values("logFC", ascending=True).index
        )

        anno_df = nhood_adata.obs[[annot_key, "logFC", "SpatialFDR"]].copy()
        anno_df["is_signif"] = anno_df["SpatialFDR"] < alpha
        anno_df = anno_df[anno_df[annot_key] != "nan"]

        cmap_df = pd.DataFrame(mdata[feature_key].var.groupby(annot_key)['logFC'].mean(), columns = ['logFC'])
        cmap = np.full(np.shape(mdata[feature_key].var.groupby(annot_key)['logFC'].mean())[0], 'lightgrey', dtype = 'object')
        cmap[mdata[feature_key].var.groupby(annot_key)['logFC'].mean() < 0] = '#377eb8'
        cmap[mdata[feature_key].var.groupby(annot_key)['logFC'].mean() > 0] = '#e41a1c'
        cmap_df['cmap'] = cmap
        if ax is None:
            _, ax = plt.subplots(1, 1, figsize = figsize, gridspec_kw={'hspace': 0.45, 'wspace': 0.4, 'bottom':0.15})
        g = sns.violinplot(
                data=anno_df,
                y=annot_key,
                x="logFC",
                order=sorted_annos,
                inner=None,
                orient="h",
                palette= cmap_df.loc[sorted_annos]['cmap'].values,
                linewidth=0,
                scale="width",
                alpha = 0.8,
                ax=ax
            )

        g = sns.stripplot(
            data=anno_df,
            y=annot_key,
            x="logFC",
            order=sorted_annos,
            hue="is_signif",
            palette=["grey", "black"],
            size=2,
            orient="h",
            alpha=0.5,
            ax = ax
        )
        g.tick_params(labelsize=12)
        g.set_xlabel('log2(fold change)', fontsize = 12)
        if label== True:
            g.set_ylabel('annotated niche neighborhoods', fontsize = 12)
        else:
            g.set_yticks([])
            g.set_ylabel('', fontsize = 12)
        if xlim is not None:
            ax.set_xlim(xlim[0], xlim[1])
        ax.axvline(x=0, ymin=0, ymax=1, color="black", linestyle="--")
        ax.legend(loc="upper right", title=f"< {int(alpha * 100)}% spatial FDR", bbox_to_anchor=(1, 1), frameon=False, prop={'size':10}, markerscale=1)
        min_name = mdata[feature_key].var.groupby(design_key)['logFC'].mean().idxmin()
        max_name = mdata[feature_key].var.groupby(design_key)['logFC'].mean().idxmax()
        plt.title(f'{min_name} vs {max_name}', fontsize = 12)
        if filename_save is not None:
            plt.savefig(os.path.join(save_directory, filename_save+'.pdf'), bbox_inches = 'tight')

def plot_niche_network(G = None, figsize = (6, 6), save_directory = 'figures',filename_save = None):
    # Now prepare the attributes for drawing
    colors = [node[1]['color'] for node in G.nodes(data=True)]
    edge_weights = [G[u][v]['weight']*0.8 for u, v in G.edges()]

    _, axes = plt.subplots(1, 1, figsize = figsize, gridspec_kw={'hspace': 0.45, 'wspace': 0.4, 'bottom':0.15})

    pos = nx.spring_layout(G)

    # Draw the graph
    nx.draw(G, pos, with_labels=True, node_color=colors, font_size = 12, node_size = 500, font_weight = 'semibold', width=edge_weights, ax = axes, edge_color = 'grey')
    axes.set_xlim([1.2*x for x in axes.get_xlim()])
    axes.set_ylim([1.2*y for y in axes.get_ylim()])
    plt.tight_layout()
    if filename_save is not None:
        plt.savefig(os.path.join(save_directory, filename_save + '.pdf'), bbox_inches = 'tight')

def plot_niche_network_donut(
    G=None, 
    figsize=(12, 12), 
    save_directory='figures', 
    filename_save=None, 
    node_order=None, 
    font_size=10, 
    buffer=1.5, 
    weightscale=0.4,
    min_node_size=50,
    max_node_size=500,
    donut_radius_inner=1.2, 
    donut_radius_outer=1.4,
    centrality_measure='degree',
    edge_color='#666A8F',  
    colors_dict=None,
    lineage_dict=None,
    curvature=0.2,
    edge_cmap=cm.viridis,
    vmin=1,
    vmax=None,
    edge_label='Patients'
):
    """
    Plots a network graph with nodes arranged in a circle and curved edges, with an edge weight color bar,
    lineage labels, and a lineage donut plot. The lineage labels are placed and oriented based on the average
    position of all nodes (wedges) belonging to that lineage, ensuring labels are never upside-down.
    """

    # Validate inputs
    if node_order is None or not isinstance(node_order, list):
        raise ValueError("node_order must be a list of nodes (e.g., cell types)")

    if lineage_dict is None or not isinstance(lineage_dict, dict):
        raise ValueError("lineage_dict must be a dictionary mapping node names to lineages")

    # Ensure all nodes in node_order are in the graph
    for node in node_order:
        if node not in G.nodes:
            G.add_node(node, lineage='Unknown', color='lightgrey')

    # Assign lineage to each node
    for node in G.nodes():
        lineage = lineage_dict.get(node, 'Unknown')
        G.nodes[node]['lineage'] = lineage

    num_nodes = len(node_order)
    pos = {node: [np.cos(2 * np.pi * i / num_nodes), np.sin(2 * np.pi * i / num_nodes)] 
           for i, node in enumerate(node_order)}

    # Calculate edge weights and normalization
    edge_weights = [G[u][v].get('weight', 1.0) * weightscale for u, v in G.edges()]
    max_patient_count = max(edge_weights) / weightscale if vmax is None else vmax
    norm = Normalize(vmin=vmin, vmax=max_patient_count)

    # Compute centrality based on the selected measure
    if centrality_measure == 'degree':
        centrality = {node: sum(G[node][neighbor].get('weight', 1.0) for neighbor in G.neighbors(node)) 
                      for node in G.nodes()}
    elif centrality_measure == 'betweenness':
        centrality = nx.betweenness_centrality(G, weight='weight')
    elif centrality_measure == 'closeness':
        centrality = nx.closeness_centrality(G, distance='weight')
    elif centrality_measure == 'eigenvector':
        try:
            centrality = nx.eigenvector_centrality(G, max_iter=1000, weight='weight')
        except nx.NetworkXException as e:
            raise ValueError(f"Error computing eigenvector centrality: {e}")
    else:
        raise ValueError(f"Unsupported centrality measure: {centrality_measure}")

    centrality_values = list(centrality.values())
    min_centrality = min(centrality_values)
    max_centrality = max(centrality_values)

    # Calculate node sizes based on centrality
    node_sizes = [
        min_node_size + (centrality[node] - min_centrality) / (max_centrality - min_centrality) * 
        (max_node_size - min_node_size) if centrality[node] >= 0.001 else 0
        for node in G.nodes()
    ]

    fig, ax = plt.subplots(figsize=figsize)

    # Draw edges with curvature and color based on weight
    for (u, v, data) in G.edges(data=True):
        x1, y1 = pos[u]
        x2, y2 = pos[v]
        edge_weight = data.get('weight', 1.0) * weightscale
        color = edge_cmap(norm(edge_weight / weightscale))
        arrow = FancyArrowPatch(
            (x1, y1), (x2, y2),
            connectionstyle=f"arc3,rad={curvature}",
            color=color,
            linewidth=edge_weight,
            antialiased=True,
            arrowstyle='-',
        )
        ax.add_patch(arrow)

    # Draw nodes
    nx.draw_networkx_nodes(
        G, pos, 
        node_color=[G.nodes[node].get('color', 'lightgrey') for node in G.nodes()], 
        node_size=node_sizes, 
        alpha=0.9, 
        edgecolors='black',
        linewidths=0.5,
        ax=ax
    )

    # Node labels aligned tangentially
    angle_tolerance = 5  # Degrees
    for node, (x, y) in pos.items():
        label_color = '#B2B0BF' if centrality[node] < 0.001 else 'black'

        # Calculate angle in radians relative to the origin
        angle_rad = np.arctan2(y, x)
        # Convert angle to degrees for rotation and normalization
        angle_deg = (np.degrees(angle_rad)) % 360

                # Define angle ranges with a tolerance to handle floating-point precision
        is_near_90 = (90 - angle_tolerance) <= angle_deg <= (90 + angle_tolerance)
        is_near_270 = (270 - angle_tolerance) <= angle_deg <= (270 + angle_tolerance)

        # Determine offset distance
        if is_near_90 or is_near_270:
            # Increase offset distance for labels near 90° and 270°
            additional_offset = 0.05  # Adjust this value as needed
        else:
            additional_offset = 0.05  # Base offset for other labels

        # Calculate offset distance
        offset_distance = donut_radius_outer + additional_offset

        # Calculate offset positions
        offset_x = offset_distance * np.cos(angle_rad)
        offset_y = offset_distance * np.sin(angle_rad)

        # Adjust alignment and rotation for labels near 90° and 270°
        if is_near_90:
            ha = 'left'
            va = 'center'
            rotation_angle = 90  # Rotate text vertically upwards
        elif is_near_270:
            ha = 'left'
            va = 'bottom'
            rotation_angle = 270  # Rotate text vertically downwards
        else:
            # General alignment for other angles
            ha = 'left' if (angle_deg < 90 or angle_deg > 270) else 'right'
            va = 'bottom' if y < 0 else 'top'

            # Adjust rotation for readability
            if 90 < angle_deg < 270:
                rotation_angle = (angle_deg + 180) % 360
            else:
                rotation_angle = angle_deg

        # ha = 'left' if (angle_deg < 90 or angle_deg > 270) else 'right'
        # va = 'bottom' if y < 0 else 'top'

        ax.text(
            offset_x, offset_y, node,
            fontsize=font_size, ha=ha, va=va,
            color=label_color,
            rotation=rotation_angle, rotation_mode='anchor'
        )

    ax.set_xlim(-1 * buffer, 1 * buffer)
    ax.set_ylim(-1 * buffer, 1 * buffer)
    ax.set_aspect('equal')
    ax.axis('off')

    # Lineage donut plot
    lineage_counts = Counter([G.nodes[node]['lineage'] for node in node_order])
    unique_lineages = list(lineage_counts.keys())
    lineage_colors_donut = {
        lineage: colors_dict.get(lineage, 'grey') if colors_dict else edge_cmap(i / len(unique_lineages))
        for i, lineage in enumerate(unique_lineages)
    }

    angles_per_node = {node: 360 * i / num_nodes for i, node in enumerate(node_order)}
    for node in node_order:
        lineage = G.nodes[node].get('lineage', 'Unknown')
        angle = angles_per_node[node]
        theta1 = angle - (360 / num_nodes) / 2
        theta2 = angle + (360 / num_nodes) / 2
        wedge = Wedge(
            center=(0, 0),
            r=donut_radius_outer,
            theta1=theta1,
            theta2=theta2,
            width=donut_radius_outer - donut_radius_inner,
            facecolor=lineage_colors_donut[lineage],
            edgecolor='none'
        )
        ax.add_patch(wedge)

    # Determine lineage label angles based on all wedges for each lineage
    lineage_angles = {}
    for lineage in unique_lineages:
        lineage_node_angles = [angles_per_node[node] for node in node_order if G.nodes[node]['lineage'] == lineage]
        mean_angle = np.mean(lineage_node_angles) % 360
        lineage_angles[lineage] = mean_angle

    # Place lineage labels around the donut
    for lineage, angle in lineage_angles.items():
        angle_rad = np.deg2rad(angle)
        label_x = (donut_radius_inner + donut_radius_outer) / 2 * np.cos(angle_rad)
        label_y = (donut_radius_inner + donut_radius_outer) / 2 * np.sin(angle_rad)

        # Baseline orientation: angle + 90° (tangent to circle)
        rotation_angle = angle + 90
        # If this puts the label upside-down (between 90° and 270°), flip 180°
        if 90 < rotation_angle <= 270:
            rotation_angle -= 180

        ax.text(
            label_x, label_y, lineage,
            ha='center', va='center', fontsize=font_size,
            rotation=rotation_angle, rotation_mode='anchor', color='white'
        )

    # Add color bar for edge weights
    cbar_ax = fig.add_axes([0.82, 0.1, 0.12, 0.02])
    cbar = ColorbarBase(cbar_ax, cmap=edge_cmap, norm=norm, orientation='horizontal')
    cbar.set_label(edge_label, fontsize=font_size)
    cbar.set_ticks([vmin, vmax])
    cbar.set_ticklabels([f"{int(vmin)}", f"{int(vmax)}"])

    plt.tight_layout()
    if filename_save is not None:
        os.makedirs(save_directory, exist_ok=True)
        plt.savefig(os.path.join(save_directory, filename_save + '.pdf'), bbox_inches='tight', dpi=1000)

    plt.show()

def cohort_cluster_plot(
    fovs: List[str],
    seg_dir: Union[pathlib.Path, str],
    save_dir: Union[pathlib.Path, str],
    cell_data: pd.DataFrame,
    fov_col: str = settings.FOV_ID,
    label_col: str = settings.CELL_LABEL,
    unassigned_col: str = '',
    cluster_col: str = settings.CELL_TYPE,
    seg_suffix: str = "_whole_cell.tiff",
    cmap: Union[str, pd.DataFrame] = "viridis",
    unassigned_color: np.ndarray = np.array([0.3, 0.3, 0.3, 1]),
    style: str = "seaborn-v0_8-paper",
    erode: bool = False,
    display_fig: bool = False,
    fig_file_type: str = "png",
    figsize: tuple = (10, 10),
    dpi: int = 300,
) -> None:
    """
    Saves the cluster masks for each FOV in the cohort as:
    - Cluster mask numbered 1-N, where N is the number of clusters (tiff)
    - Cluster mask colored by cluster with or without a colorbar (png)
    - Cluster mask colored by cluster (tiff).
    """

    plt.style.use(style)

    if isinstance(seg_dir, str):
        seg_dir = pathlib.Path(seg_dir)

    try:
        io_utils.validate_paths(seg_dir)
    except ValueError:
        raise ValueError(f"Could not find the segmentation directory at {seg_dir.as_posix()}")

    if isinstance(save_dir, str):
        save_dir = pathlib.Path(save_dir)
        if not save_dir.exists():
            save_dir.mkdir(parents=True, exist_ok=True)
    if isinstance(fovs, str):
        fovs = [fovs]

    # Create directories
    for sub_dir in ["cluster_masks", "cluster_masks_colored", "cluster_plots"]:
        (save_dir / sub_dir).mkdir(parents=True, exist_ok=True)

    cmd = ClusterMaskData(
        data=cell_data,
        fov_col=fov_col,
        label_col=label_col,
        cluster_col=cluster_col,
    )

    if isinstance(cmap, pd.DataFrame):
        unique_clusters: pd.DataFrame = cmd.mapping[[cmd.cluster_column,
                                                     cmd.cluster_id_column]].drop_duplicates()
        cmap_colors: pd.DataFrame = cmap.merge(
            right=unique_clusters,
            on=cmd.cluster_column
        ).sort_values(by="cluster_id")["color"].values
        colors_like: list[bool] = [colors.is_color_like(c) for c in cmap_colors]

        if not all(colors_like):
            bad_color_values = cmap_colors[~np.array(colors_like)]
            raise ValueError(
                ("Not all colors in the provided cmap are valid colors. "
                 f"The following colors are invalid: {bad_color_values}")
            )

        np_colors = colors.to_rgba_array(cmap_colors)
        color_map, norm = create_cmap(np_colors, n_clusters=cmd.n_clusters, unassigned_color=unassigned_color)

    elif isinstance(cmap, str):
        color_map, norm = create_cmap(cmap, n_clusters=cmd.n_clusters, unassigned_color=unassigned_color)

    # Generate cluster masks for each FOV
    with tqdm(total=len(fovs), desc="Cluster Mask Generation", unit="FOVs") as pbar:
        for fov in fovs:
            pbar.set_postfix(FOV=fov)

            # Generate the cell mask for the FOV
            cluster_mask: np.ndarray = generate_cluster_mask(
                fov=fov,
                seg_dir=seg_dir,
                cmd=cmd,
                seg_suffix=seg_suffix,
                erode=erode,
            )

            # Save the cluster mask (numbered)
            save_fov_mask(
                fov,
                data_dir=save_dir / "cluster_masks",
                mask_data=cluster_mask,
                sub_dir=None,
            )

            # Save the colored cluster mask
            save_colored_mask(
                fov=fov,
                save_dir=save_dir / "cluster_masks_colored",
                suffix=".tiff",
                data=cluster_mask,
                cmap=color_map,
                norm=norm,
            )

            cluster_labels = ["Background"] + cmd.cluster_names + ["Unassigned"]

            fig = plot_cluster(
                image=cluster_mask,
                fov=fov,
                cmap=color_map,
                norm=norm,
                cbar_visible=True,
                cbar_labels=cluster_labels,
                figsize=figsize,
                dpi=dpi,
            )

            fig.savefig(
                fname=os.path.join(save_dir, "cluster_plots", f"{fov}.{fig_file_type}"),
            )

            if display_fig:
                fig.show(warn=False)
            else:
                plt.close(fig)

            pbar.update(1)


def create_cmap(
    cmap: Union[np.ndarray, list[str], str],
    n_clusters: int,
    unassigned_color: np.ndarray = np.array([0.5, 0.5, 0.5, 0.5])
) -> tuple[colors.ListedColormap, colors.BoundaryNorm]:
    """
    Creates a discrete colormap and a boundary norm from the provided colors.

    Args:
        cmap (Union[np.ndarray, list[str], str]): The colormap, or set of colors to use.
        n_clusters (int): The number of clusters for the colormap.
        unassigned_color (np.ndarray, optional): RGBA array for unassigned pixels.
            Defaults to black ([0.0, 0.0, 0.0, 1.0]).

    Returns:
        tuple[colors.ListedColormap, colors.BoundaryNorm]:
            The generated colormap and boundary norm.
    """

    if isinstance(cmap, np.ndarray):
        if cmap.ndim != 2:
            raise ValueError(
                f"cmap array must be a 2D array, got {cmap.ndim}D array"
            )
        if cmap.shape[0] != n_clusters:
            raise ValueError(
                f"cmap array must have {n_clusters} colors, got {cmap.shape[0]} colors"
            )

        color_map = colors.ListedColormap(
            colors=_cmap_add_background_unassigned(cmap, unassigned_color=unassigned_color)
        )

    elif isinstance(cmap, list):
        if len(cmap) != n_clusters:
            raise ValueError(
                f"cmap list must have {n_clusters} colors, got {len(cmap)} colors"
            )
        np_colors = colors.to_rgba_array(cmap)
        color_map = colors.ListedColormap(
            colors=_cmap_add_background_unassigned(np_colors, unassigned_color=unassigned_color)
        )

    elif isinstance(cmap, str):
        try:
            color_map_obj = colormaps[cmap]
        except KeyError:
            raise KeyError(f"Colormap {cmap} not found.")

        colors_rgba: np.ndarray = color_map_obj(np.linspace(0, 1, n_clusters))
        color_map = colors.ListedColormap(
            colors=_cmap_add_background_unassigned(colors_rgba, unassigned_color=unassigned_color)
        )

    else:
        raise TypeError("cmap must be a np.ndarray, list[str], or str")

    # Create boundaries and a norm for the ListedColormap
    bounds = [i - 0.5 for i in np.linspace(0, color_map.N, color_map.N + 1)]
    norm = colors.BoundaryNorm(bounds, color_map.N)
    return color_map, norm

def plot_cluster(
        image: np.ndarray,
        fov: str,
        cmap: colors.ListedColormap,
        norm: colors.BoundaryNorm,
        cbar_visible: bool = True,
        cbar_labels: list[str] = None,
        dpi: int = 300,
        figsize: tuple[int, int] = None) -> Figure:
    """
    Plots the cluster image with the provided colormap and norm.

    Args:
        image (np.ndarray):
            The cluster image to plot.
        fov (str):
            The name of the clustered FOV.
        cmap (colors.ListedColormap):
            A colormap to use for the cluster image.
        norm (colors.BoundaryNorm):
            A normalization to use for the cluster image.
        cbar_visible (bool, optional):
            Whether or not to display the colorbar. Defaults to True.
        cbar_labels (list[str], optional):
            Colorbar labels for the clusters. Devaults to None, where
            the labels will be automatically generated.
        dpi (int, optional):
            The resolution of the image to use for saving. Defaults to 300.
        figsize (tuple, optional):
            The size of the image to display. Defaults to (10, 10).

    Returns:
        Figure: Returns the cluster image as a matplotlib Figure.
    """
    # Default colorbar labels
    if cbar_labels is None:
        cbar_labels = [f"Cluster {x}" for x in range(1, len(cmap.colors))]

    fig: Figure = plt.figure(figsize=figsize, dpi=dpi)
    fig.set_layout_engine(layout="tight")
    gs = gridspec.GridSpec(nrows=1, ncols=1, figure=fig)
    fig.suptitle(f"{fov}")

    # Image axis
    ax: Axes = fig.add_subplot(gs[0, 0])
    ax.axis("off")
    ax.grid(visible=False)

    ax.imshow(
        X=image,
        cmap=cmap,
        norm=norm,
        origin="upper",
        aspect="equal",
        interpolation="none",
    )

    if cbar_visible:
        # # Manually set the colorbar
        divider = make_axes_locatable(fig.gca())
        cax = divider.append_axes(position="right", size="5%", pad="3%")

        cbar = fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap),
                            cax=cax, orientation="vertical", use_gridspec=True, pad=0.1,
                            shrink=0.9, drawedges=True)
        cbar.ax.set_yticks(
            ticks=np.arange(len(cbar_labels)),
            labels=cbar_labels
        )
        cbar.minorticks_off()

    return fig

def _cmap_add_background_unassigned(
    cluster_colors: np.ndarray,
    unassigned_color: np.ndarray = np.array([0.3, 0.3, 0.3, 1]),
    background_color: np.ndarray = np.array([0.0, 0.0, 0.0, 1.0])
) -> np.ndarray:
    """
    Inserts background and unassigned colors into the colormap array.

    Args:
        cluster_colors (np.ndarray): An array of shape (N, 4), where N is the number of cluster colors.
        unassigned_color (np.ndarray): An RGBA array for the unassigned color.
        background_color (np.ndarray): An RGBA array for the background color.

    Returns:
        np.ndarray: A new color array including background and unassigned colors.
    """
    return np.vstack([background_color, cluster_colors, unassigned_color])


def save_colored_mask(
    fov: str,
    save_dir: str,
    suffix: str,
    data: np.ndarray,
    cmap: colors.ListedColormap,
    norm: colors.BoundaryNorm,
) -> None:
    """Saves the colored mask to the provided save directory.

    Args:
        fov (str):
            The name of the FOV.
        save_dir (str):
            The directory where the colored mask will be saved.
        suffix (str):
            The suffix to append to the FOV name.
        data (np.ndarray):
            The mask to save.
        cmap (colors.ListedColormap):
            The colormap to use for the mask.
        norm (colors.BoundaryNorm):
            The normalization to use for the mask.
    """

    # Create the save directory if it does not exist
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)

    # Create the colored mask
    colored_mask = img_as_ubyte(cmap(norm(data)))

    # Save the image
    image_utils.save_image(
        fname=os.path.join(save_dir, f"{fov}{suffix}"),
        data=colored_mask,
    )

def plot_diff_func(df, total, labels_key, functional_markers, figsize =(6, 4.5), cmap = 'vlag', vmin = -1, vmax = 1, vcenter = 0, filename_save = 'func'):
    # Group by cell_cluster and calculate mean expression for each marker
    subset_means = df.groupby(labels_key).mean().loc[:, functional_markers]
    total_means = total.groupby(labels_key).mean().loc[:, functional_markers]

    # Calculate the difference in means
    mean_diff = subset_means - total_means

    # Remove rows with any NaN values
    mean_diff.dropna(inplace=True)

    adata_run = anndata.AnnData(mean_diff)
    adata_run.obs[labels_key] = pd.Categorical(mean_diff.index)
    adata_run.obs_names = [f'c_{i}' for i in range(0, len(adata_run.obs_names))]

    fig, axes = plt.subplots(1, 1, figsize=figsize, dpi = 400)

    sc.pl.matrixplot(adata_run, 
                        var_names=functional_markers, 
                        groupby=labels_key,
                        dendrogram = True,
                        vmin=vmin, vmax=vmax, vcenter = vcenter, cmap=cmap, 
                        colorbar_title='standard deviation', 
                        ax = axes,
                        return_fig=False,
                        save = filename_save)
    
def plot_overlay(seg_dir, data_dir, fov, channels, nuclei_channels, channel_to_rgb, save_directory, filename_save):
    img = None
    inst_seg = np.squeeze(io.imread(os.path.join(seg_dir, fov + '_whole_cell.tiff'))).astype(np.float32)
    boundaries = find_boundaries(inst_seg, mode='inner')
    inst_seg[boundaries] = 0
    inst_seg = (inst_seg > 0).astype(np.float32)

    img = np.zeros((inst_seg.shape[0], inst_seg.shape[1], 3), dtype=np.float32)
    for idx, channel in enumerate(channels):
        data_path = os.path.join(data_dir, fov, channel + '.tiff')
        channel_img = np.squeeze(io.imread(data_path))
        channel_img = channel_img / np.quantile(channel_img, 0.99)
        channel_img = np.clip(channel_img, 0, 1)

        # color mapping
        img[..., 0] += channel_img * channel_to_rgb[idx, 0]  # Red component
        img[..., 1] += channel_img * channel_to_rgb[idx, 1]  # Green component
        img[..., 2] += channel_img * channel_to_rgb[idx, 2]  # Blue component

    nuclei_img = []
    for nuclei_chan in nuclei_channels:
        data_path = os.path.join(data_dir, fov, nuclei_chan + '.tiff')
        nuclei_img += [np.squeeze(io.imread(data_path))]

    nuclei_img = np.stack(nuclei_img, axis=-1)
    nuclei_img = nuclei_img / np.quantile(nuclei_img, 0.99)
    nuclei_img = nuclei_img.mean(axis=-1)

    fig, ax = plt.subplots(1, 1, figsize=(5, 5))
    ax.imshow(img, vmin=0, vmax=1.0, interpolation="none")
    ax.imshow(nuclei_img, cmap="gray", vmin=0, vmax=1.0, interpolation="none", alpha = 0.6)
    ax.axis("off")

    # position and style the labels
    label_spacing = 0.03  # vertical space
    label_x_spacing = 0.05  # horizontal space
    label_position_y = 1 - label_spacing  # y position
    # label_position_x = 0.01  # x position start
    label_position_x = 0.25

    ax.text(label_position_x, label_position_y, "Nuc", color="gray", transform=ax.transAxes,
            fontsize=12, fontweight='bold', ha='left', va='center', path_effects=[patheffects.withStroke(linewidth=0.5, foreground='black')])
    label_position_x += len("Nuc") * 0.025 + label_x_spacing  # Adjust spacing
    
    for idx, channel in enumerate(channels):
        color = channel_to_rgb[idx]
        ax.text(label_position_x, label_position_y, channel, color=color, transform=ax.transAxes,
                fontsize=12, fontweight='bold', ha='left', va='center', path_effects=[patheffects.withStroke(linewidth=0.5, foreground='black')])
        label_position_x += len(channel) * 0.025 + label_x_spacing  # Adjust spacing

    # scale bar
    scale_bar_length = int((nuclei_img.shape[0]*100)/800)  # 2048 pixels is 800 uM so 256 pixels corresponds to 100 microns
    scale_bar_height = 20  # Height of the scale bar rectangle
    scale_bar_color = 'white'
    text_y_offset = 30  # Offset of the text from the scale bar

    ax.add_patch(Rectangle((img.shape[1] - scale_bar_length - 80, img.shape[0] - scale_bar_height - 20),
                            scale_bar_length, scale_bar_height, linewidth=0, edgecolor=None, facecolor=scale_bar_color))

    ax.text(img.shape[1] - scale_bar_length / 2 - 80, img.shape[0] - scale_bar_height - text_y_offset,
            '100 µm', color=scale_bar_color, fontweight='regular', fontsize=12, ha='center', va='bottom', path_effects=[patheffects.withStroke(linewidth=0.5, foreground='black')])
    
    plt.savefig(os.path.join(save_directory, f'{filename_save}.svg'), dpi = 600)
    plt.show()
    plt.close()