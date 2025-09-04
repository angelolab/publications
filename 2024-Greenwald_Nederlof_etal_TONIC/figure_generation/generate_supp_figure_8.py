import os

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import colors
import skimage.io as io
import seaborn as sns

import python_files.supplementary_plot_helpers as supplementary_plot_helpers
from ark.utils.plot_utils import color_segmentation_by_stat

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42


BASE_DIR = "/Volumes/Shared/Noah Greenwald/TONIC_Cohort/"
SUPPLEMENTARY_FIG_DIR = os.path.join(BASE_DIR, "figures_revision/supp_figs_panels")
INTERMEDIATE_DIR = os.path.join(BASE_DIR, "intermediate_files")
seg_dir = os.path.join(BASE_DIR, "segmentation_data/deepcell_output")
image_dir = os.path.join(BASE_DIR, "image_data/samples/")
dist_mat_dir = os.path.join(INTERMEDIATE_DIR, "spatial_analysis", "dist_mats")


## 2.9 Feature correlation plot ##
# cluster features together to identify modules
fov_data_df = pd.read_csv(os.path.join(BASE_DIR, 'analysis_files/combined_feature_data_filtered.csv'))
fov_data_wide = fov_data_df.pivot(index='fov', columns='feature_name_unique', values='normalized_value')
corr_df = fov_data_wide.corr(method='spearman')
corr_df = corr_df.fillna(0)

clustergrid = sns.clustermap(corr_df, cmap='vlag', vmin=-1, vmax=1, figsize=(20, 20))
matrix_order = clustergrid.dendrogram_row.reordered_ind

new_tick_positions = [(i+1)-0.5 for i in range(len(matrix_order)) if (i+1)%10==0]
new_tick_labels = [str(i+1) for i in range(len(matrix_order)) if (i+1)%10==0]
clustergrid.ax_heatmap.set_xticks(new_tick_positions)
clustergrid.ax_heatmap.set_yticks(new_tick_positions)
clustergrid.ax_heatmap.set_xticklabels(new_tick_labels)
clustergrid.ax_heatmap.set_yticklabels(new_tick_labels)

clustergrid.savefig(os.path.join(os.path.join(SUPPLEMENTARY_FIG_DIR, 'supp_figure_8a.pdf')), dpi=300)
plt.close()

# Feature parameter tuning (panels b, c, d)
#extraction_pipeline_tuning_dir = os.path.join(SUPPLEMENTARY_FIG_DIR, "supp_figure_7_robustness")
extraction_pipeline_tuning_dir = os.path.join(SUPPLEMENTARY_FIG_DIR, "supp_figure_8_robustness")
if not os.path.exists(extraction_pipeline_tuning_dir):
    os.makedirs(extraction_pipeline_tuning_dir)

# vary the features for each marker threshold
cell_table_full = pd.read_csv(
    os.path.join(BASE_DIR, "analysis_files/combined_cell_table_normalized_cell_labels_updated.csv")
)
supplementary_plot_helpers.run_functional_marker_positivity_tuning_tests(
    cell_table_full, extraction_pipeline_tuning_dir, supplementary_plot_helpers.MARKER_INFO,
    threshold_mults=[1/4, 1/2, 3/4, 7/8, 1, 8/7, 4/3, 2, 4]
)

# vary min cell param to see how many FOVs get kept or not
cluster_broad_df = pd.read_csv(os.path.join(BASE_DIR, "output_files/cluster_df_per_core.csv"))
supplementary_plot_helpers.run_min_cell_feature_gen_fovs_dropped_tests(
    cluster_broad_df, min_cell_params=[1, 3, 5, 10, 20], compartments=["all"],
    metrics=["cluster_broad_count"], save_dir=extraction_pipeline_tuning_dir
)

# vary params for cancer mask and boundary definition inclusion
cell_table_clusters = pd.read_csv(os.path.join(BASE_DIR, "analysis_files/cell_table_clusters.csv"))
supplementary_plot_helpers.run_cancer_mask_inclusion_tests(
    cell_table_clusters, channel_dir=image_dir, seg_dir=seg_dir,
    threshold_mults=[1/4, 1/2, 3/4, 7/8, 1, 8/7, 4/3, 2, 4],
    save_dir=extraction_pipeline_tuning_dir, base_sigma=10, base_channel_thresh=0.0015,
    base_min_mask_size=7000, base_max_hole_size=1000, base_border_size=50
)


# Neighborhood/diversity tuning (panels e, f)
neighbors_mat_dir = os.path.join(BASE_DIR, "supplementary_figs", "supp_figure_13_data")
diversity_mixing_pipeline_tuning_dir = SUPPLEMENTARY_FIG_DIR

# vary the parameters for each marker threshold
cell_table_full = pd.read_csv(os.path.join(BASE_DIR, "analysis_files/cell_table_clusters.csv"))
supplementary_plot_helpers.run_diversity_mixing_tuning_tests(
    cell_table_full, dist_mat_dir=dist_mat_dir,
    neighbors_mat_dir=neighbors_mat_dir,
    save_dir=diversity_mixing_pipeline_tuning_dir,
    threshold_mults=[1/4, 1/2, 3/4, 7/8, 1, 8/7, 4/3, 2, 4],
    mixing_info=supplementary_plot_helpers.MIXING_INFO,
    base_pixel_radius=50, cell_type_col="cell_cluster_broad"
)


### R CODE ##
# Make heatmaps for panel g
#
# library(data.table)
# library(pheatmap)
# library(viridis)
# library(RColorBrewer)
#
# som_mean_path = "pixie/ecm_061423_pixel_output_dir/pixel_channel_avg_som_cluster.csv"
# meta_mean_path = "pixie/ecm_061423_pixel_output_dir/pixel_channel_avg_meta_cluster.csv"
# colors_path = "cluster_colors.csv"
# cap = 3
#
# # Phenotype to color mapping
# colors_tab = fread(colors_path)
# mat_colors = colors_tab$color
# names(mat_colors) = colors_tab$pixel_meta_cluster_rename
# mat_colors = list(clust = mat_colors)
#
# # Blue-white-red colors
# rwb_cols = colorRampPalette(c("royalblue4","white","red4"))(99)
#
# pdf("pixelClustering_heatmaps.pdf",height=8,width=7)
# ## Heatmap of pixel clusters x markers, average across pixel clusters
# mean_dat = fread(som_mean_path)
# mean_dat = mean_dat[order(pixel_meta_cluster_rename)]
# marker_names = colnames(mean_dat)
# marker_names = marker_names[!marker_names %in% c("pixel_som_cluster","count","pixel_meta_cluster","pixel_meta_cluster_rename")]
# mat_dat = data.frame(mean_dat[, ..marker_names])
# rownames(mat_dat) = paste0("clust_", mean_dat$pixel_som_cluster)
# # Z-score and cap
# mat_dat = scale(mat_dat)
# mat_dat = pmin(mat_dat, cap)
# # Annotations
# mat_col = data.frame(clust = mean_dat$pixel_meta_cluster_rename)
# rownames(mat_col) = paste0("clust_", mean_dat$pixel_som_cluster)
# # Determine breaks
# breaks = seq(-cap, cap, length.out=100)
# # Make heatmap
# pheatmap(mat_dat,
#          color = rwb_cols,
#          breaks = breaks,
#          cluster_rows = FALSE,
#          show_rownames = FALSE,
#          annotation_row = mat_col,
#          annotation_colors = mat_colors,
#          main = "Average across 100 SOM clusters")
#
# ## Heatmap of pixel hierarchical cluster x markers, average across hierarchical clusters
# mean_dat = fread(meta_mean_path)
# mat_dat = data.frame(mean_dat[,..marker_names])
# rownames(mat_dat) = mean_dat$pixel_meta_cluster_rename
# # Z-score the columns
# mat_dat = scale(mat_dat)
# # Make annotations
# mat_col = data.frame(clust = mean_dat$pixel_meta_cluster_rename)
# rownames(mat_col) = mean_dat$pixel_meta_cluster_rename
# # Determine breaks
# breaks = seq(-cap, cap, length.out=100)
# # Make heatmap
# pheatmap(mat_dat,
#          color = rwb_cols,
#          breaks = breaks,
#          annotation_row = mat_col,
#          annotation_colors = mat_colors,
#          main = "Average across meta clusters")
# dev.off()


# ECM pixel cluster masks panel h
pixie_dir = "pixie/ecm_061423_pixel_output_dir"
pixel_mask_dir = os.path.join(pixie_dir, "pixel_masks")
clust_to_pheno_path = os.path.join(pixie_dir, "pixel_channel_avg_meta_cluster.csv")

output_dir_colored = "colored_masks"
colors_path = "cluster_colors.csv"
clust_coln = "pixel_meta_cluster_rename"

# Make output directory
if not os.path.exists(output_dir_colored):
    os.makedirs(output_dir_colored)

# Get phenotype mapping
clust_to_pheno = pd.read_csv(clust_to_pheno_path)
clust_to_pheno = clust_to_pheno[['pixel_meta_cluster', 'pixel_meta_cluster_rename']]
colors_tab = pd.read_csv(colors_path)
clust_to_color = pd.merge(clust_to_pheno, colors_tab, on='pixel_meta_cluster_rename')
clust_to_color = clust_to_color.sort_values(['pixel_meta_cluster'])
maxk = max(clust_to_color['pixel_meta_cluster'])

# Fill in missing
full_df = pd.DataFrame({'pixel_meta_cluster': list(range(1, maxk + 1))})
clust_to_color = pd.merge(clust_to_color, full_df, on='pixel_meta_cluster', how='right')
clust_to_color.fillna('#000000', inplace=True)

## Create custom cmap
mycols = list(clust_to_color["color"])
mycols.insert(0, '#000000')  # first color is black
# Make bounds
bounds = [i - 0.5 for i in np.linspace(0, maxk + 1, maxk + 2)]
colmap = colors.ListedColormap(mycols)
norm = colors.BoundaryNorm(bounds, colmap.N)

all_fovs = os.listdir(pixel_mask_dir)
all_fovs = [x for x in all_fovs if "_pixel_mask.tiff" in x]
all_fovs = [x.replace("_pixel_mask.tiff", "") for x in all_fovs]
for fov in all_fovs:
    # Read in pixel mask
    clust_array = np.array(io.imread(os.path.join(pixel_mask_dir, fov + "_pixel_mask.tiff")))

    # Save colored overlay
    image = colmap(norm(clust_array))
    plt.imsave(os.path.join(output_dir_colored, fov + "_pixel_mask_colored.tiff"), image)
    print(fov)


## 3.7 Fiber visualizations (i. j)
# color fibers by alignment & length stats
fiber_table = pd.read_csv(os.path.join(INTERMEDIATE_DIR, 'fiber_segmentation_processed_data/fiber_object_table.csv'))
fibseg_dir = os.path.join(BASE_DIR, 'supplementary_figs/review_figures/fiber_features/fiber_masks')

feature_fovs = {
    'alignment_score': ['TONIC_TMA14_R8C3', 'TONIC_TMA4_R10C1'],
    'major_axis_length': ['TONIC_TMA7_R9C4', 'TONIC_TMA13_R8C1']
}
for metric in feature_fovs:
    fov_list = feature_fovs[metric]
    fiber_table_sub = fiber_table[fiber_table.fov.isin(fov_list)]
    fiber_table_sub[f'{metric}_norm'] = np.log(fiber_table_sub[metric])
    if metric == 'alignment_score':
        fiber_table_sub[f'{metric}_norm'] = fiber_table_sub[f'{metric}_norm'] + 0.5
    metric = f'{metric}_norm'

    save_dir = os.path.join(SUPPLEMENTARY_FIG_DIR, 'supp_figure_8ij', f'colored_{metric}')
    os.makedirs(save_dir, exist_ok=True)
    color_segmentation_by_stat(
        fovs=fiber_table_sub.fov.unique(), data_table=fiber_table_sub, seg_dir=fibseg_dir, save_dir=save_dir,
        stat_name=metric, cmap="Blues", seg_suffix="_fiber_labels.tiff", erode=True, fig_file_type='pdf')
