import os
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc 
import quiche
from matplotlib.colors import to_rgba
from matplotlib.cm import ScalarMappable
import quiche as qu
import os
import pathlib
import shutil
import logging
from dataclasses import dataclass, field
from operator import contains
from typing import Dict, List, Literal, Optional, Tuple, Union
from matplotlib import gridspec
from matplotlib.axes import Axes

import matplotlib.cm as cm
import matplotlib.colors as colors
from matplotlib import colormaps, patches
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import natsort
import numpy as np
import pandas as pd
from pandas.core.groupby.generic import DataFrameGroupBy
import skimage
import xarray as xr
from alpineer import image_utils, io_utils, load_utils, misc_utils
from alpineer.settings import EXTENSION_TYPES
from mpl_toolkits.axes_grid1 import make_axes_locatable
from skimage.exposure import rescale_intensity
from skimage import io
from matplotlib.gridspec import GridSpec
import matplotlib.ticker as ticker
import networkx as nx
from skimage.util import img_as_ubyte
from tqdm.auto import tqdm
from ark import settings
from skimage.segmentation import find_boundaries
from ark.utils.data_utils import (
    ClusterMaskData,
    erode_mask,
    generate_cluster_mask,
    save_fov_mask,
    map_segmentation_labels,
)


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

def plot_neighbors_histogram_per_fov(adata,
                                        radius=200,
                                        p=2,
                                        figcomp = [6,5],
                                        xlim = [0,450],
                                        ylim = [0, 2000],
                                        figsize = (15, 16),
                                        fov_subset = None,
                                        spatial_key='X_spatial',
                                        fov_key='fov',
                                        save_directory = 'figures',
                                        filename_save = None):
    
    fig, axes = plt.subplots(nrows=figcomp[0], ncols=figcomp[1], figsize=figsize, gridspec_kw={'hspace': 0.5, 'wspace': 0.3, 'bottom':0.15})
    sns.set_style('ticks')
    # Flatten the 2D array of subplots into a 1D array
    axes_flat = axes.flatten()
    nn_counts_per_fov = []
    # Loop through unique fields of view (FOV)
    for i, fov in enumerate(fov_subset):
        adata_fov = adata[np.where(adata.obs[fov_key] == fov)[0], :].copy()
        spatial_kdTree = cKDTree(adata_fov.obsm[spatial_key])

        # Query the KD-tree to find neighbors within the specified radius
        nn = spatial_kdTree.query_ball_point(adata_fov.obsm[spatial_key], r=radius, p=p)

        # Count the number of neighbors for each cell
        nn_counts = [len(neighbors) for neighbors in nn]
        nn_counts_per_fov.append(nn_counts)
        # Plot the histogram for the current FOV
        axes_flat[i].hist(nn_counts)
        axes_flat[i].set_title(f'{fov}')
        axes_flat[i].set_xlabel('Number of Neighbors', fontsize = 8)
        axes_flat[i].set_ylabel('Frequency', fontsize = 8)
        axes_flat[i].tick_params(labelsize=10)
        axes_flat[i].set_xlim(xlim[0], xlim[1])
        axes_flat[i].set_ylim(ylim[0], ylim[1])
        mean_counts = np.mean(nn_counts)
        axes_flat[i].axvline(mean_counts, color='red', linestyle='--', linewidth=1, label=f'Mean: {mean_counts:.2f}')
        axes_flat[i].legend()

    # Adjust layout
    plt.tight_layout()
    if filename_save is not None:
        plt.savefig(os.path.join(save_directory, filename_save + '.pdf'), bbox_inches = 'tight')
    return nn_counts_per_fov

def plot_nhood_graph(
        mdata,
        alpha = 0.1,
        min_logFC = 0,
        min_size = 10,
        plot_edges = False,
        save = None,
        ax = None,
        vmin = -10,
        vmax = 10,
        save_directory = None,
        filename_save = None,
        xlim = None,
        ylim = None,
        **kwargs,
    ) -> None:
        """Visualize DA results on abstracted graph (wrapper around sc.pl.embedding)

        Args:
            mdata: MuData object
            alpha: Significance threshold. (default: 0.1)
            min_logFC: Minimum absolute log-Fold Change to show results. If is 0, show all significant neighbourhoods. (default: 0)
            min_size: Minimum size of nodes in visualization. (default: 10)
            plot_edges: If edges for neighbourhood overlaps whould be plotted. Defaults to False.
            title: Plot title. Defaults to "DA log-Fold Change".
            show: Show the plot, do not return axis.
            save: If `True` or a `str`, save the figure. A string is appended to the default filename.
                  Infer the filetype if ending on {`'.pdf'`, `'.png'`, `'.svg'`}.
            **kwargs: Additional arguments to `scanpy.pl.embedding`.
        """
        nhood_adata = mdata["milo"].T.copy()

        if "Nhood_size" not in nhood_adata.obs.columns:
            raise KeyError(
                'Cannot find "Nhood_size" column in adata.uns["nhood_adata"].obs -- \
                    please run milopy.utils.build_nhood_graph(adata)'
            )

        nhood_adata.obs["graph_color"] = nhood_adata.obs["logFC"]
        nhood_adata.obs.loc[nhood_adata.obs["SpatialFDR"] > alpha, "graph_color"] = np.nan
        nhood_adata.obs["abs_logFC"] = abs(nhood_adata.obs["logFC"])
        nhood_adata.obs.loc[nhood_adata.obs["abs_logFC"] < min_logFC, "graph_color"] = np.nan

        # Plotting order - extreme logFC on top
        nhood_adata.obs.loc[nhood_adata.obs["graph_color"].isna(), "abs_logFC"] = np.nan
        ordered = nhood_adata.obs.sort_values("abs_logFC", na_position="first").index
        nhood_adata = nhood_adata[ordered]

        vmax = np.max([nhood_adata.obs["graph_color"].max(), abs(nhood_adata.obs["graph_color"].min())])
        vmin = -vmax

        size = np.array(nhood_adata.obs["Nhood_size"] * min_size)
        size[np.isnan(nhood_adata.obs["Nhood_size"] * min_size)] = 0.5

        g = sc.pl.embedding(
            nhood_adata,
            "X_milo_graph",
            color="graph_color",
            cmap="RdBu_r",
            edges=plot_edges,
            neighbors_key="nhood",
            sort_order=False,
            frameon=False,
            size = size,
            vmax=vmax,
            vmin=vmin,
            save=save,
            ax=ax,
            show = False,
            **kwargs)
        
        if xlim is not None: 
            g.set_xlim(-7, 22.5)
        if ylim is not None:
            g.set_ylim(-10, 22.5)
        plt.savefig(os.path.join(save_directory, filename_save+'.pdf'), bbox_inches = 'tight')

def da_beeswarm(mdata,
                feature_key = "rna",
                alpha: float = 0.1,
                subset_nhoods = None,
                figsize = (6, 12),
                niche_key = None,
                design_key = 'condition',
                patient_key = 'sample',
                xlim = None,
                percentile = 70,
                save_directory = 'figures',
                filename_save = None):
        """Plot beeswarm plot of logFC against nhood labels

        Args:
            mdata: MuData object
            anno_col: Column in adata.uns['nhood_adata'].obs to use as annotation. (default: 'nhood_annotation'.)
            alpha: Significance threshold. (default: 0.1)
            subset_nhoods: List of nhoods to plot. If None, plot all nhoods. (default: None)
            palette: Name of Seaborn color palette for violinplots.
                     Defaults to pre-defined category colors for violinplots.

        Examples:
            >>> import pertpy as pt
            >>> import scanpy as sc
            >>> adata = pt.dt.bhattacherjee()
            >>> milo = pt.tl.Milo()
            >>> mdata = milo.load(adata)
            >>> sc.pp.neighbors(mdata["rna"])
            >>> milo.make_nhoods(mdata["rna"])
            >>> mdata = milo.count_nhoods(mdata, sample_col="orig.ident")
            >>> milo.da_nhoods(mdata, design="~label")
            >>> milo.annotate_nhoods(mdata, anno_col='cell_type')
            >>> pt.pl.milo.da_beeswarm(mdata)
        """
        try:
            nhood_adata = mdata["milo"].T.copy()
            nhood_adata.obs[[patient_key, design_key]] = mdata[feature_key][mdata[feature_key].obs['nhood_ixs_refined'] == 1].obs[[patient_key, design_key]].values
        except KeyError:
            raise RuntimeError(
                "mdata should be a MuData object with two slots: feature_key and 'milo'. Run 'milopy.count_nhoods(adata)' first."
            ) from None

        if subset_nhoods is not None:
            nhood_adata = nhood_adata[subset_nhoods]

        try:
            nhood_adata.obs[niche_key]
        except KeyError:
            raise RuntimeError(
                f"Unable to find {niche_key} in mdata.uns['nhood_adata']. Run 'milopy.utils.annotate_nhoods(adata, anno_col)' first"
            ) from None

        try:
            nhood_adata.obs["logFC"]
        except KeyError:
            raise RuntimeError(
                "Unable to find 'logFC' in mdata.uns['nhood_adata'].obs. Run 'core.da_nhoods(adata)' first."
            ) from None

        sorted_annos = (
            nhood_adata.obs[[niche_key, "logFC"]].groupby(niche_key).mean().sort_values("logFC", ascending=True).index
        )

        anno_df = nhood_adata.obs[[niche_key, "logFC", "SpatialFDR", patient_key, design_key]].copy()
        anno_df["is_signif"] = anno_df["SpatialFDR"] < alpha
        anno_df = anno_df[anno_df[niche_key] != "nan"]

        mdata['milo'].var[design_key] = mdata['rna'].obs[design_key].values


        adata_runner = mdata['rna'][mdata['rna'].obs['nhood_ixs_refined'] == 1].copy()
        adata_runner.obs[niche_key] = mdata['milo'].var.loc[:, niche_key].values
        adata_runner = adata_runner[np.isin(adata_runner.obs[niche_key], sorted_annos)]
        adata_runner.obs[niche_key] = pd.Categorical(adata_runner.obs[niche_key], categories=sorted_annos)
        perc  = quiche.pp.compute_percentile(np.abs(mdata['milo'].var.groupby(niche_key)['logFC'].mean()[sorted_annos]), p = percentile)
        cmap_df = pd.DataFrame(mdata['milo'].var.groupby(niche_key)['logFC'].mean(), columns = ['logFC'])
        cmap = np.full(np.shape(mdata['milo'].var.groupby(niche_key)['logFC'].mean())[0], 'lightgrey', dtype = 'object')
        cmap[mdata['milo'].var.groupby(niche_key)['logFC'].mean() < -1*perc] = '#377eb8'
        cmap[mdata['milo'].var.groupby(niche_key)['logFC'].mean() >= perc] = '#e41a1c'
        cmap_df['cmap'] = cmap
        # if percentile != 0:
        #     sorted_annos = sorted_annos[np.where((cmap_df.loc[sorted_annos]['logFC'] < -1*perc) | (cmap_df.loc[sorted_annos]['logFC'] >= perc))[0]]
        

        _, ax = plt.subplots(1, 1, figsize = figsize, gridspec_kw={'hspace': 0.45, 'wspace': 0.4, 'bottom':0.15})
        g = sns.violinplot(
                data=anno_df,
                y=niche_key,
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
            y=niche_key,
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
        g.set_ylabel('annotated niche neighborhoods', fontsize = 12)
        if xlim is not None:
            ax.set_xlim(xlim[0], xlim[1])
        ax.axvline(x=0, ymin=0, ymax=1, color="black", linestyle="--")
        ax.legend(loc="upper right", title=f"< {int(alpha * 100)}% spatial FDR", bbox_to_anchor=(1, 1), frameon=False, prop={'size':10}, markerscale=1)
        min_name = mdata['milo'].var.groupby(design_key)['logFC'].mean().idxmin()
        max_name = mdata['milo'].var.groupby(design_key)['logFC'].mean().idxmax()
        plt.title(f'{min_name} vs {max_name}', fontsize = 12)
        plt.savefig(os.path.join(save_directory, filename_save+'.pdf'), bbox_inches = 'tight')
        return sorted_annos

def plot_neighborhood_dist(mdata, save_directory, filename_save):
    _, axes = plt.subplots(1, 1, figsize = (4, 3.5), gridspec_kw={'hspace': 0.45, 'wspace': 0.4, 'bottom':0.15})
    sns.set_style('ticks')
    nhood_size = np.array(mdata['rna'].obsm["nhoods"].sum(0)).ravel()
    axes.hist(nhood_size, bins=50)
    axes.tick_params(labelsize=12)
    plt.xlabel('number of cells in neighborhood', fontsize = 12)
    plt.ylabel('number of neighborhoods', fontsize = 12)
    plt.savefig(os.path.join(save_directory, filename_save+'.pdf'), bbox_inches = 'tight')

def plot_grid_enrichment(df, hue_order, colors_dict, selected_grid_list, num_grids_x, num_grids_y, save_directory, filename_save):
    _, axes = plt.subplots(nrows=1, ncols=1, figsize=(6,6))
    # Plot the points with colors representing the labels
    sns.scatterplot(x = 'x', y = 'y', hue = 'group', data =df, alpha=0.5, palette=colors_dict, hue_order = hue_order, ax = axes)
    axes.tick_params(labelsize=10)
    # Set axis labels and legend
    plt.xlabel('Y', fontsize = 12)
    plt.ylabel('X', fontsize = 12)
    plt.xlim(0, 1)
    plt.ylim(0,1)
    plt.legend(loc='upper right', bbox_to_anchor=(1.15, 1.0))
        # Outline the selected grid box in red
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

    # Pivot the data for plotting
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

logging.getLogger('fontTools.subset').setLevel(logging.WARNING)

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
        """Plot beeswarm plot of logFC against nhood labels

        Args:
            mdata: MuData object
            anno_col: Column in adata.uns['nhood_adata'].obs to use as annotation. (default: 'nhood_annotation'.)
            alpha: Significance threshold. (default: 0.1)
            subset_nhoods: List of nhoods to plot. If None, plot all nhoods. (default: None)
            palette: Name of Seaborn color palette for violinplots.
                     Defaults to pre-defined category colors for violinplots.

        Examples:
            >>> import pertpy as pt
            >>> import scanpy as sc
            >>> adata = pt.dt.bhattacherjee()
            >>> milo = pt.tl.Milo()
            >>> mdata = milo.load(adata)
            >>> sc.pp.neighbors(mdata["rna"])
            >>> milo.make_nhoods(mdata["rna"])
            >>> mdata = milo.count_nhoods(mdata, sample_col="orig.ident")
            >>> milo.da_nhoods(mdata, design="~label")
            >>> milo.annotate_nhoods(mdata, anno_col='cell_type')
            >>> pt.pl.milo.da_beeswarm(mdata)
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
        df_pivot.iloc[:, 1] = df_pivot.iloc[:, 1] * -1 #Either 0 or 1 here for flipping bargraph 
        df_pivot = df_pivot.loc[sorted_annos, :]

        # Plot the proportions on the twin axis
        ax1.barh(df_pivot.index, df_pivot.iloc[:, 0],  height = 0.5,color=colors_dict[df_pivot.iloc[:, 0].name], edgecolor='none', label=df_pivot.iloc[:, 0].name)
        ax1.barh(df_pivot.index, df_pivot.iloc[:, 1], height = 0.5, color=colors_dict[df_pivot.iloc[:, 1].name], edgecolor='none', left=0, label=df_pivot.iloc[:, 1].name)

        # Adjust the twin axis limits to make sure bars are centered
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
        # ticks = np.arange(xlim_prev[0], xlim_prev[1]+0.2, 0.2)
        # ax1.set_xticks(ticks)

        # # Determine labels based on xlim_prev
        # labels = [str(round(t, 2)) if i % 2 == 0 else '' for i, t in enumerate(ticks)]

        # ax1.set_xticklabels(labels)

        # # Change x-tick labels to be positive
        xticks = ax1.get_xticks()
        ax1.set_xticklabels([str(abs(np.round(x, 2))) for x in xticks])

        if xlim is not None:
            ax0.set_xlim(xlim[0], xlim[1])

        ax0.axvline(x=0, ymin=0, ymax=1, color="black", linestyle="--", linewidth = 1)
        ax1.axvline(x=0, ymin=0, ymax=1, color="black", linestyle="--", linewidth = 1)
        ax0.legend(loc="upper right", title=f"< {int(alpha * 100)}% spatial FDR", bbox_to_anchor=(1, 1), frameon=False, prop={'size':fontsize}, markerscale=1, fontsize = fontsize)

        # plt.subplots_adjust(wspace=0.12)
        min_name = mdata[feature_key].var.groupby(design_key)['logFC'].mean().idxmin()
        max_name = mdata[feature_key].var.groupby(design_key)['logFC'].mean().idxmax()
        ax0.set_title(f'{min_name} vs {max_name}', fontsize = fontsize)
        if filename_save is not None:
            plt.savefig(os.path.join(save_directory, filename_save+'.pdf'), bbox_inches = 'tight')

# #beeswarm plots 
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
            filename_save = None):
        """Plot beeswarm plot of logFC against nhood labels

        Args:
            mdata: MuData object
            anno_col: Column in adata.uns['nhood_adata'].obs to use as annotation. (default: 'nhood_annotation'.)
            alpha: Significance threshold. (default: 0.1)
            subset_nhoods: List of nhoods to plot. If None, plot all nhoods. (default: None)
            palette: Name of Seaborn color palette for violinplots.
                     Defaults to pre-defined category colors for violinplots.

        Examples:
            >>> import pertpy as pt
            >>> import scanpy as sc
            >>> adata = pt.dt.bhattacherjee()
            >>> milo = pt.tl.Milo()
            >>> mdata = milo.load(adata)
            >>> sc.pp.neighbors(mdata["rna"])
            >>> milo.make_nhoods(mdata["rna"])
            >>> mdata = milo.count_nhoods(mdata, sample_col="orig.ident")
            >>> milo.da_nhoods(mdata, design="~label")
            >>> milo.annotate_nhoods(mdata, anno_col='cell_type')
            >>> pt.pl.milo.da_beeswarm(mdata)
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
        g.set_ylabel('annotated niche neighborhoods', fontsize = 12)
        if xlim is not None:
            ax.set_xlim(xlim[0], xlim[1])
        ax.axvline(x=0, ymin=0, ymax=1, color="black", linestyle="--")
        ax.legend(loc="upper right", title=f"< {int(alpha * 100)}% spatial FDR", bbox_to_anchor=(1, 1), frameon=False, prop={'size':10}, markerscale=1)
        min_name = mdata[feature_key].var.groupby(design_key)['logFC'].mean().idxmin()
        max_name = mdata[feature_key].var.groupby(design_key)['logFC'].mean().idxmax()
        plt.title(f'{min_name} vs {max_name}', fontsize = 12)
        if filename_save is not None:
            plt.savefig(os.path.join(save_directory, filename_save+'.pdf'), bbox_inches = 'tight')

def plot_niche_network(G=None, font_size=12, figsize=(6, 6), save_directory='figures', filename_save=None, k=0.1):
    # Now prepare the attributes for drawing
    colors = [node[1]['color'] for node in G.nodes(data=True)]
    edge_weights = [G[u][v]['weight']*0.8 for u, v in G.edges()]

    _, axes = plt.subplots(1, 1, figsize=figsize, gridspec_kw={'hspace': 0.45, 'wspace': 0.4, 'bottom': 0.15})

    # Adjust the spring layout with the specified k value for more spacing
    pos = nx.spring_layout(G, k=k)

    # Draw the graph with text
    nx.draw(G, pos, with_labels=True, node_color=colors, font_size=font_size, node_size=500, font_weight='semibold', width=edge_weights, ax=axes, edge_color='grey')
    axes.set_xlim([1.2*x for x in axes.get_xlim()])
    axes.set_ylim([1.2*y for y in axes.get_ylim()])
    plt.tight_layout()
    
    # Save the version with text
    if filename_save is not None:
        plt.savefig(os.path.join(save_directory, filename_save + '.svg'), bbox_inches='tight')

    # Clear the plot without closing the figure
    axes.clear()

    # Draw the graph without text
    nx.draw(G, pos, with_labels=False, node_color=colors, node_size=500, width=edge_weights, ax=axes, edge_color='grey')
    axes.set_xlim([1.2*x for x in axes.get_xlim()])
    axes.set_ylim([1.2*y for y in axes.get_ylim()])
    plt.tight_layout()

    # Save the version without text
    if filename_save is not None:
        plt.savefig(os.path.join(save_directory, filename_save + '_no_text.svg'), bbox_inches='tight')


def boxplot_only(mdata,
                 feature_key="quiche",
                 alpha: float = 0.05,
                 niches=None,
                 figsize=(6, 12),
                 annot_key='quiche_niche',
                 design_key='condition',
                 patient_key='sample',
                 xlim=None,
                 fontsize=8,
                 xlim_prev = [-1,1],
                 colors=['#377eb8', '#e41a1c'],
                 save_directory='figures',
                 filename_save=None,
                 csv_filename=None):
    """Plot beeswarm plot of logFC against nhood labels

    Args:
        mdata: MuData object
        anno_col: Column in adata.uns['nhood_adata'].obs to use as annotation. (default: 'nhood_annotation'.)
        alpha: Significance threshold. (default: 0.1)
        subset_nhoods: List of nhoods to plot. If None, plot all nhoods. (default: None)
        palette: Name of Seaborn color palette for violinplots.
                 Defaults to pre-defined category colors for violinplots.

    Examples:
        >>> import pertpy as pt
        >>> import scanpy as sc
        >>> adata = pt.dt.bhattacherjee()
        >>> milo = pt.tl.Milo()
        >>> mdata = milo.load(adata)
        >>> sc.pp.neighbors(mdata["rna"])
        >>> milo.make_nhoods(mdata["rna"])
        >>> mdata = milo.count_nhoods(mdata, sample_col="orig.ident")
        >>> milo.da_nhoods(mdata, design="~label")
        >>> milo.annotate_nhoods(mdata, anno_col='cell_type')
        >>> pt.pl.milo.da_beeswarm(mdata)
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
    gs = GridSpec(1, 2, width_ratios=[1, 0.4])  # 2 columns with different width ratios

    # Plot the boxplot in the first column
    ax0 = plt.subplot(gs[0])
    # ax1 = plt.subplot(gs[1])

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


    ax0.axvline(x=0, ymin=0, ymax=1, color="black", linestyle="--", linewidth = 1)
    plt.subplots_adjust(wspace=0.12)
    min_name = mdata[feature_key].var.groupby(design_key)['logFC'].mean().idxmin()
    max_name = mdata[feature_key].var.groupby(design_key)['logFC'].mean().idxmax()
    ax0.set_title(f'{min_name} vs {max_name}', fontsize=fontsize)
    if filename_save is not None:
        plt.savefig(os.path.join(save_directory, filename_save + '.pdf'), bbox_inches='tight')
    if csv_filename is not None:
        anno_df.to_csv(os.path.join(save_directory, csv_filename + '.csv'), index=False)

def boxplot_prev(mdata,
                 feature_key="quiche",
                 alpha: float = 0.05,
                 niches=None,
                 figsize=(6, 12),
                 annot_key='quiche_niche',
                 design_key='condition',
                 patient_key='sample',
                 xlim=None,
                 fontsize=8,
                 xlim_prev = [-1,1],
                 colors=['#377eb8', '#e41a1c'],
                 save_directory='figures',
                 filename_save=None):
    """Plot beeswarm plot of logFC against nhood labels

    Args:
        mdata: MuData object
        anno_col: Column in adata.uns['nhood_adata'].obs to use as annotation. (default: 'nhood_annotation'.)
        alpha: Significance threshold. (default: 0.1)
        subset_nhoods: List of nhoods to plot. If None, plot all nhoods. (default: None)
        palette: Name of Seaborn color palette for violinplots.
                 Defaults to pre-defined category colors for violinplots.

    Examples:
        >>> import pertpy as pt
        >>> import scanpy as sc
        >>> adata = pt.dt.bhattacherjee()
        >>> milo = pt.tl.Milo()
        >>> mdata = milo.load(adata)
        >>> sc.pp.neighbors(mdata["rna"])
        >>> milo.make_nhoods(mdata["rna"])
        >>> mdata = milo.count_nhoods(mdata, sample_col="orig.ident")
        >>> milo.da_nhoods(mdata, design="~label")
        >>> milo.annotate_nhoods(mdata, anno_col='cell_type')
        >>> pt.pl.milo.da_beeswarm(mdata)
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
    gs = GridSpec(1, 2, width_ratios=[1, 0.4])  # 2 columns with different width ratios

    # Plot the boxplot in the first column
    ax0 = plt.subplot(gs[0])
    ax1 = plt.subplot(gs[1])

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

    cov_count_df = qu.tl.compute_patient_proportion(mdata, niches=niches, feature_key=feature_key, annot_key=annot_key, patient_key=patient_key, design_key=design_key, patient_niche_threshold=0)

    # Pivot the dataframe for the stacked barplot
    df_pivot = cov_count_df.pivot_table(index=annot_key, columns=design_key, values='prop_cov')
    df_pivot = df_pivot.loc[cov_count_df[[annot_key, 'med_logFC']].sort_values(by='med_logFC', ascending=False).drop_duplicates()[annot_key].values]
    df_pivot.iloc[:, 0] = df_pivot.iloc[:, 0] * -1
    df_pivot = df_pivot.loc[sorted_annos, :]

    # Plot the proportions on the twin axis
    ax1.barh(df_pivot.index, df_pivot.iloc[:, 0], height = 0.5, color=colors[0], edgecolor='none', label=df_pivot.iloc[:, 0].name)
    ax1.barh(df_pivot.index, df_pivot.iloc[:, 1], height  = 0.5, color=colors[1], edgecolor='none', left=0, label=df_pivot.iloc[:, 1].name)
    # Adjust the twin axis limits to make sure bars are centered
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

    # Change x-tick labels to be positive
    xticks = ax1.get_xticks()
    ax1.set_xticklabels([str(abs(x)) for x in xticks])

    ax0.set_xlabel('log2(fold change)', fontsize=fontsize)
    ax0.set_ylabel('cell type neighborhoods', fontsize=fontsize)
    if xlim is not None:
        ax0.set_xlim(xlim[0], xlim[1])

    ticks = np.arange(xlim_prev[0], xlim_prev[1]+0.5, 0.5)
    ax1.set_xticks(ticks)

        # Determine labels based on xlim_prev
    labels = [str(round(t, 2)) if i % 2.5 == 0 else '' for i, t in enumerate(ticks)]

    ax1.set_xticklabels(labels)

        # Change x-tick labels to be positive
    xticks = ax1.get_xticks()
    ax1.set_xticklabels([str(abs(np.round(x, 2))) for x in xticks])


    ax0.axvline(x=0, ymin=0, ymax=1, color="black", linestyle="--", linewidth = 1)
    ax1.axvline(x=0, ymin=0, ymax=1, color="black", linestyle="--", linewidth = 1)
    plt.subplots_adjust(wspace=0.12)
    min_name = mdata[feature_key].var.groupby(design_key)['logFC'].mean().idxmin()
    max_name = mdata[feature_key].var.groupby(design_key)['logFC'].mean().idxmax()
    ax0.set_title(f'{min_name} vs {max_name}', fontsize=fontsize)
    if filename_save is not None:
        plt.savefig(os.path.join(save_directory, filename_save + '.pdf'), bbox_inches='tight')


def plot_stacked_bar(data, design_key, network_label_neg, network_label_pos, save_directory, file_name, exclude_label=None):
    # Filter the data for each status
    if exclude_label is not None:
        neg_data = data[(data[design_key] == network_label_neg) & (data['cell_meta_cluster_final'] != exclude_label)]
        pos_data = data[(data[design_key] == network_label_pos) & (data['cell_meta_cluster_final'] != exclude_label)]
    else:
        neg_data = data[data[design_key] == network_label_neg]
        pos_data = data[data[design_key] == network_label_pos]

    # Calculate the percentage of each cell_meta_cluster_final for 'neg'
    neg_percentage = neg_data['cell_meta_cluster_final'].value_counts(normalize=True) 
    neg_percentage = neg_percentage.sort_index()

    # Calculate the percentage of each cell_meta_cluster_final for 'pos'
    pos_percentage = pos_data['cell_meta_cluster_final'].value_counts(normalize=True)
    pos_percentage = pos_percentage.sort_index()

    # Create a DataFrame for plotting
    percentage_df = pd.DataFrame({
        'neg': neg_percentage,
        'pos': pos_percentage
    }).fillna(0)

    # Define a fixed color mapping for each group
    unique_clusters = data['cell_meta_cluster_final'].unique()
    num_clusters = len(unique_clusters)
    colors = plt.cm.get_cmap('tab20', num_clusters).colors
    color_mapping = {cluster: colors[i] for i, cluster in enumerate(unique_clusters)}

    # Map colors to the DataFrame's index
    color_list = [color_mapping[cluster] for cluster in percentage_df.index]

    # Plot the stacked bar graphs
    fig, ax = plt.subplots(figsize=(10, 7))

    percentage_df.T.plot(kind='bar', stacked=True, ax=ax, color=color_list)
    ax.set_title('Percentage of cell_meta_cluster_final')
    ax.set_xlabel('Status')
    ax.set_ylabel('Percentage')
    ax.legend(title='cell_meta_cluster_final', bbox_to_anchor=(1.05, 1), loc='upper left')

    plt.tight_layout()
    
    # Save the figure as a PDF
    pdf_path = os.path.join(save_directory, file_name)
    fig.savefig(pdf_path, bbox_inches='tight', pad_inches=0.1)
    plt.show()

def plot_circle_network(G=None, font_size=12, figsize=(6, 6), save_directory='figures', filename_save=None):
    # Prepare the attributes for drawing
    colors = [node[1]['color'] for node in G.nodes(data=True)]
    edge_weights = [G[u][v]['weight'] * 0.8 for u, v in G.edges()]

    _, axes = plt.subplots(1, 1, figsize=figsize)

    # Use circular layout for the graph
    pos = nx.circular_layout(G)

    # Draw the graph with text
    nx.draw(G, pos, with_labels=True, node_color=colors, font_size=font_size, node_size=500, font_weight='semibold',
            width=edge_weights, ax=axes, edge_color='grey')
    plt.tight_layout()

    # Save the version with text
    if filename_save is not None:
        plt.savefig(os.path.join(save_directory, filename_save + '.svg'), bbox_inches='tight')

    # # Clear the plot without closing the figure
    # axes.clear()

    # # Draw the graph without text
    # nx.draw(G, pos, with_labels=False, node_color=colors, node_size=500, width=edge_weights, ax=axes, edge_color='grey')
    # plt.tight_layout()

    # Save the version without text
    if filename_save is not None:
        plt.savefig(os.path.join(save_directory, filename_save + '_no_text.svg'), bbox_inches='tight')
# def boxplot_prev(mdata,
#             feature_key = "quiche",
#             alpha: float = 0.05,
#             niches = None,
#             figsize = (6, 12),
#             annot_key = 'quiche_niche',
#             design_key = 'condition',
#             patient_key = 'sample',
#             xlim = None,
#             fontsize = 8,
#             xlim_prev =[0, 0.5],
#             colors = ['#377eb8', '#e41a1c'],
#             save_directory = 'figures',
#             filename_save = None):
#         """Plot beeswarm plot of logFC against nhood labels

#         Args:
#             mdata: MuData object
#             anno_col: Column in adata.uns['nhood_adata'].obs to use as annotation. (default: 'nhood_annotation'.)
#             alpha: Significance threshold. (default: 0.1)
#             subset_nhoods: List of nhoods to plot. If None, plot all nhoods. (default: None)
#             palette: Name of Seaborn color palette for violinplots.
#                      Defaults to pre-defined category colors for violinplots.

#         Examples:
#             >>> import pertpy as pt
#             >>> import scanpy as sc
#             >>> adata = pt.dt.bhattacherjee()
#             >>> milo = pt.tl.Milo()
#             >>> mdata = milo.load(adata)
#             >>> sc.pp.neighbors(mdata["rna"])
#             >>> milo.make_nhoods(mdata["rna"])
#             >>> mdata = milo.count_nhoods(mdata, sample_col="orig.ident")
#             >>> milo.da_nhoods(mdata, design="~label")
#             >>> milo.annotate_nhoods(mdata, anno_col='cell_type')
#             >>> pt.pl.milo.da_beeswarm(mdata)
#         """
#         try:
#             nhood_adata = mdata[feature_key].T.copy()
#         except KeyError:
#             raise RuntimeError(
#                 "mdata should be a MuData object with two slots: feature_key and 'milo'. Run 'milopy.count_nhoods(adata)' first."
#             ) from None

#         if niches is not None:
#             nhood_adata = nhood_adata[np.isin(nhood_adata.obs[annot_key], niches)]

#         try:
#             nhood_adata.obs[annot_key]
#         except KeyError:
#             raise RuntimeError(
#                 f"{annot_key} not defined"
#             ) from None

#         try:
#             nhood_adata.obs["logFC"]
#         except KeyError:
#             raise RuntimeError(
#                 "Unable to find 'logFC' in mdata.uns['nhood_adata'].obs. Run 'core.da_nhoods(adata)' first."
#             ) from None

#         sorted_annos = (
#             nhood_adata.obs[[annot_key, "logFC"]].groupby(annot_key).median().sort_values("logFC", ascending=True).index
#         )

#         anno_df = nhood_adata.obs[[annot_key, "logFC", "SpatialFDR"]].copy()
#         anno_df["is_signif"] = anno_df["SpatialFDR"] < alpha
#         anno_df = anno_df[anno_df[annot_key] != "nan"]

#         cmap_df = pd.DataFrame(mdata[feature_key].var.groupby(annot_key)['logFC'].mean(), columns = ['logFC'])
#         cmap = np.full(np.shape(mdata[feature_key].var.groupby(annot_key)['logFC'].mean())[0], 'lightgrey', dtype = 'object')
#         cmap[mdata[feature_key].var.groupby(annot_key)['logFC'].mean() <= -1] = colors[0]
#         cmap[mdata[feature_key].var.groupby(annot_key)['logFC'].mean() > 1] = colors[1]
#         cmap_df['cmap'] = cmap
#         fig = plt.figure(figsize=figsize)
#         gs = GridSpec(1, 2, width_ratios=[1, 0.4])  # 2 columns with equal width

#         # Plot the strip plot in the first column

#         ax0 = plt.subplot(gs[0])
#         ax1 = plt.subplot(gs[1])
        
#         g = sns.boxplot(
#                 data=anno_df,
#                 y=annot_key,
#                 x="logFC",
#                 order=sorted_annos,
#                 orient="h",
#                 palette= cmap_df.loc[sorted_annos]['cmap'].values,
#                 ax=ax0,
#                 width = 0.6,
#                 fliersize=0.8
#             )

#         cov_count_df = qu.tl.compute_patient_proportion(mdata, niches = niches, feature_key = feature_key, annot_key = annot_key, patient_key = patient_key, design_key = design_key, patient_niche_threshold = 3)
#         #Sort the dataframe

#         #Pivot the dataframe for stacked barplot
#         df_pivot = cov_count_df.pivot_table(index=annot_key, columns=design_key, values='prop_cov')
#         df_pivot = df_pivot.loc[cov_count_df[[annot_key, 'med_logFC']].sort_values(by = 'med_logFC', ascending = False).drop_duplicates()[annot_key].values]
#         df_pivot['0'] = df_pivot['0']*-1
#         df_pivot = df_pivot.loc[sorted_annos, :]
#         # df_pivot.iloc[np.where((df_pivot < 0.05) & (df_pivot > 0))] = 0.05

#         print(df_pivot)

    
#         # Plot the proportions on the twin axis
#         ax1.barh(df_pivot.index, df_pivot['0'], color=colors[0], edgecolor='none', label='Relapse 0')
#         ax1.barh(df_pivot.index, df_pivot['1'], color=colors[1], edgecolor='none', left=0, label='Relapse 1')

#         # Adjust the twin axis limits to make sure bars are centered
#         ax1.set_xlim(-1, 1)

#         # df_pivot.plot(kind='barh', stacked=True, ax = ax1, color = colors)

#         ax1.set_yticklabels([])
#         ax1.set_yticks([])
#         ax1.set_xlabel('proportion of patients', fontsize = fontsize)
#         ax1.set_ylabel('')
#         # ax1.set_xlim([0, 0.6])
#         ax1.legend(loc="upper right", title=design_key, frameon=False, bbox_to_anchor=(1, 1), prop={'size':fontsize}, markerscale=1,fontsize = fontsize)
#         ax0.tick_params(labelsize=8)
#         ax1.tick_params(labelsize=8)

#         ax0.set_xlabel('log2(fold change)', fontsize = fontsize)
#         ax0.set_ylabel('cell type neighborhoods', fontsize = fontsize)
#         if xlim is not None:
#             ax0.set_xlim(xlim[0], xlim[1])

#         # if xlim_prev is not None:
#         #     ax1.set_xlim(xlim_prev[0], xlim_prev[1])

#         ax0.axvline(x=0, ymin=0, ymax=1, color="black", linestyle="--", linewidth = 1)
#         # ax0.axvline(x=-1, ymin=0, ymax=1, color="black", linestyle="--", linewidth = 1)
#         # ax0.axvline(x=1, ymin=0, ymax=1, color="black", linestyle="--", linewidth = 1)
#         plt.subplots_adjust(wspace=0.12)
#         min_name = mdata[feature_key].var.groupby(design_key)['logFC'].mean().idxmin()
#         max_name = mdata[feature_key].var.groupby(design_key)['logFC'].mean().idxmax()
#         ax0.set_title(f'{min_name} vs {max_name}', fontsize = fontsize)
#         if filename_save is not None:
#             plt.savefig(os.path.join(save_directory, filename_save+'.pdf'), bbox_inches = 'tight')