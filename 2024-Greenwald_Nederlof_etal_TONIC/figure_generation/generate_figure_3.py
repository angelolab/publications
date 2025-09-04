import os

import numpy as np
import pandas as pd
import matplotlib

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

import matplotlib.pyplot as plt
import seaborn as sns
from ark.utils.plot_utils import cohort_cluster_plot
import ark.settings as settings
import skimage.io as io


base_dir = '/Volumes/Shared/Noah Greenwald/TONIC_Cohort'
metadata_dir = os.path.join(base_dir, 'intermediate_files/metadata')
seg_dir = os.path.join(base_dir, 'segmentation_data/deepcell_output')
image_dir = os.path.join(base_dir, 'image_data/samples/')
plot_dir = os.path.join(base_dir, 'figures')

# load files
harmonized_metadata = pd.read_csv(os.path.join(metadata_dir, 'harmonized_metadata.csv'))
feature_metadata = pd.read_csv(os.path.join(base_dir, 'analysis_files/feature_metadata.csv'))
study_fovs = harmonized_metadata.loc[harmonized_metadata.Timepoint.isin(['primary', 'baseline', 'pre_nivo', 'on_nivo']), 'fov'].values
ranked_features_all = pd.read_csv(os.path.join(base_dir, 'analysis_files/feature_ranking.csv'))
ranked_features = ranked_features_all.loc[ranked_features_all.comparison.isin(['primary', 'baseline', 'pre_nivo', 'on_nivo'])]


# plot total volcano
fig, ax = plt.subplots(figsize=(3,3))
# color pallete options: Greys, magma, vlag, icefire
sns.scatterplot(data=ranked_features, x='med_diff', y='log_pval', alpha=1, hue='importance_score', palette=sns.color_palette("icefire", as_cmap=True),
                s=2.5, edgecolor='none', ax=ax)
ax.set_xlim(-3, 3)
ax.set_ylim(0, 8)
sns.despine()

# add gradient legend
norm = plt.Normalize(ranked_features.importance_score.min(), ranked_features.importance_score.max())
sm = plt.cm.ScalarMappable(cmap="icefire", norm=norm)
ax.get_legend().remove()
ax.figure.colorbar(sm, ax=ax)
plt.tight_layout()

plt.savefig(os.path.join(plot_dir, 'Figure3a_volcano.pdf'))
plt.close()

text_df = ranked_features.loc[ranked_features.feature_name_unique.isin(['PDL1+__APC', 'T__cell_cluster_broad_density', 'Cancer_diversity']), :]


# compare ratio features to best individual feature that is part of the ratio
top_ratios = ranked_features.iloc[:100, :]
top_ratios = top_ratios.loc[top_ratios.feature_type.isin(['density_ratio']), :]

cell_types, rankings, feature_num, score = [], [], [], []
for idx, row in top_ratios.iterrows():
    ratio_cell_types = row.feature_name.split('__')[:2]

    # find individual densities with same compartment, timepoint, and cell types as the ratio
    candidate_features = ranked_features.loc[ranked_features.compartment == row.compartment, :]
    candidate_features = candidate_features.loc[candidate_features.feature_type == 'density', :]
    candidate_features = candidate_features.loc[candidate_features.comparison == row.comparison, :]
    candidate_features = candidate_features.loc[np.logical_or(candidate_features.feature_name.str.contains(ratio_cell_types[0]),
                                                              candidate_features.feature_name.str.contains(ratio_cell_types[1])), :]
    candidate_features = candidate_features.loc[candidate_features.cell_pop_level == 'cell_cluster_broad', :]

    # save values for best density
    feature_num.append(candidate_features.shape[0])
    best_rank = candidate_features.combined_rank.min()
    candidate_features = candidate_features.loc[candidate_features.combined_rank == best_rank, :]

    if candidate_features.shape[0] != 0:
        # just take first in case there are ties
        candidate_features = candidate_features.iloc[0:1, :]

        cell_types.append(candidate_features.feature_name_unique.values)
        rankings.append(candidate_features.combined_rank.values)
        score.append(candidate_features.importance_score.values)
    else:
        cell_types.append(['None'])
        rankings.append([0])
        score.append([0])

# create comparison DF with both ratios and densities
comparison_df = pd.DataFrame({'cell_type': np.concatenate(cell_types), 'rank': np.concatenate(rankings), 'feature_num': feature_num, 'density_score': np.concatenate(score)})
comparison_df['original_ranking'] = top_ratios.combined_rank.values
comparison_df['original_feature'] = top_ratios.feature_name_unique.values
comparison_df['ratio_score'] = top_ratios.importance_score.values
comparison_df = comparison_df.loc[comparison_df.feature_num == 2, :]
comparison_df['feature_id'] = np.arange(comparison_df.shape[0])

plot_df = pd.melt(comparison_df, id_vars=['feature_id'], value_vars=['density_score', 'ratio_score'])
plot_df['variable'] = pd.Categorical(plot_df['variable'], categories=['ratio_score', 'density_score'], ordered=True)
# remove border from dots
fig, ax = plt.subplots(figsize=(4, 3))
sns.lineplot(data=plot_df, x='variable', y='value', units='feature_id', estimator=None,
             color='grey', alpha=0.5, marker='o', markeredgewidth=0, markersize=5, ax=ax) # markeredgecolor='none'
sns.despine()
plt.ylim(0, 1)
plt.tight_layout()
plt.savefig(os.path.join(plot_dir, 'Figure3b_ratio_vs_density.pdf'))
plt.close()

# annotate specific features
text_df = comparison_df.loc[comparison_df.original_feature == 'T__Cancer__ratio__cancer_core', :]
text_df = ranked_features.loc[ranked_features.feature_name_unique == 'Cancer__cell_cluster_broad_density__cancer_core', :]

# look at enrichment by compartment
top_counts = ranked_features.iloc[:100, :].groupby('compartment').count().iloc[:, 0]
total_counts = feature_metadata.groupby('compartment').count().iloc[:, 0]

# calculate abundance of each compartment in the top 100 and across all features
top_prop = top_counts / np.sum(top_counts)
total_prop = total_counts / np.sum(total_counts)
top_ratio = top_prop / total_prop
top_ratio = np.log2(top_ratio)

# create df
ratio_df = pd.DataFrame({'compartment': top_ratio.index, 'ratio': top_ratio.values})
ratio_df = ratio_df.sort_values(by='ratio', ascending=False)

fig, ax = plt.subplots(figsize=(4, 3))
sns.barplot(data=ratio_df, x='compartment', y='ratio', color='grey', ax=ax)
sns.despine()
ax.set_ylim(-0.7, 1.5)
plt.savefig(os.path.join(plot_dir, 'Figure3c_enrichment_by_compartment.pdf'))
plt.close()

print(top_counts)

# look at enrichment of spatial features
spatial_features = ['mixing_score', 'cell_diversity', 'compartment_area_ratio', 'pixie_ecm',
                    'compartment_area', 'fiber', 'linear_distance', 'ecm_fraction', 'ecm_cluster', 'kmeans_cluster']
spatial_mask = np.logical_or(ranked_features.feature_type.isin(spatial_features), ranked_features.compartment != 'all')
ranked_features['spatial_feature'] = spatial_mask

spatial_mask_metadata = np.logical_or(feature_metadata.feature_type.isin(spatial_features), feature_metadata.compartment != 'all')
feature_metadata['spatial_feature'] = spatial_mask_metadata

# calculate proportion of spatial features in top 100 vs all features
top_count_spatial = ranked_features.iloc[:100, :].groupby('spatial_feature').count().iloc[:, 0]
total_counts_spatial = feature_metadata.groupby('spatial_feature').count().iloc[:, 0]
top_prop = top_count_spatial / np.sum(top_count_spatial)
total_prop = total_counts_spatial / np.sum(total_counts_spatial)

top_ratio = top_prop / total_prop
top_ratio = np.log2(top_ratio)
ratio_df = pd.DataFrame({'spatial_feature': top_ratio.index, 'ratio': top_ratio.values})
ratio_df = ratio_df.sort_values(by='ratio', ascending=False)

fig, ax = plt.subplots(figsize=(4, 3))
ax.set_ylim(-0.7, 1.5)
sns.barplot(data=ratio_df, x='spatial_feature', y='ratio', color='grey', ax=ax)
sns.despine()

plt.savefig(os.path.join(plot_dir, 'Figure3d_enrichment_by_spatial.pdf'))
plt.close()

print(top_count_spatial)

# plot top features
plot_features = ranked_features.copy()

# create csv file with necessary metadata to make it easy to reorganize based on categories
plot_features['ratio'] = plot_features.feature_type.isin(['density_ratio', 'density_proportion'])
plot_features['density'] = plot_features.feature_type == 'density'
plot_features['diversity'] = plot_features.feature_type.isin(['region_diversity', 'cell_diversity'])
plot_features['phenotype'] = plot_features.feature_type == 'functional_marker'
plot_features['sign'] = plot_features.med_diff > 0
plot_features = plot_features.iloc[:52, :]
plot_features = plot_features[['feature_name', 'feature_name_unique', 'compartment', 'ratio', 'density', 'diversity', 'phenotype', 'sign']]
plot_features = plot_features.drop_duplicates()
plot_features_sort = plot_features.sort_values(by='feature_name')
plot_features_sort.to_csv(os.path.join(plot_dir, 'Figure3e_top_hits.csv'))


# PDL1+__CD68 macs on nivo example
combined_df = pd.read_csv(os.path.join(base_dir, 'analysis_files/timepoint_combined_features_outcome_labels.csv'))

feature_name = 'PDL1+__CD68_Mac'
timepoint = 'on_nivo'

plot_df = combined_df.loc[(combined_df.feature_name_unique == feature_name) &
                                    (combined_df.Timepoint == timepoint), :]

fig, ax = plt.subplots(1, 1, figsize=(2, 4))
sns.stripplot(data=plot_df, x='Clinical_benefit', y='raw_mean', order=['Yes', 'No'],
                color='black', ax=ax)
sns.boxplot(data=plot_df, x='Clinical_benefit', y='raw_mean', order=['Yes', 'No'],
                color='grey', ax=ax, showfliers=False, width=0.3)
ax.set_title(feature_name + ' ' + timepoint)
ax.set_ylim([0, 1])
sns.despine()
plt.tight_layout()
plt.savefig(os.path.join(plot_dir, 'Figure3f_feature_{}_{}.pdf'.format(feature_name, timepoint)))
plt.close()

# example overlays
cell_table_func = pd.read_csv(os.path.join(base_dir, 'analysis_files/cell_table_func_single_positive.csv'))
cell_table_func['M1_plot'] = cell_table_func.cell_cluster
cell_table_func.loc[cell_table_func.cell_cluster != 'CD68_Mac', 'M1_plot'] = 'Other'
cell_table_func.loc[(cell_table_func.cell_cluster == 'CD68_Mac') & (cell_table_func.PDL1.values), 'M1_plot'] = 'CD68_Mac__PDL1+'

m1_colormap = pd.DataFrame({'M1_plot': ['CD68_Mac', 'Other', 'CD68_Mac__PDL1+'],
                         'color': ['blue','grey', 'lightsteelblue']})

# generate overlays for representative patients to identify CD68 Mac PDL1+ cells
# subset = plot_df.loc[plot_df.raw_mean < 0.08, :]
# pats = [37, 33, 59, 62, 64, 65] # responders
# pats = [24, 60, 87, 88, 107, 114] # nonresponders
#
# fovs = harmonized_metadata.loc[(harmonized_metadata.Patient_ID.isin(pats) & harmonized_metadata.MIBI_data_generated.values), 'fov'].unique()
#
# cell_table_subset = cell_table_subset.loc[(cell_table_subset.fov.isin(fovs)), :]
#

# m1_plot_dir = os.path.join(plot_dir, 'Figure3f_M1_overlays')
# if not os.path.exists(m1_plot_dir):
#     os.mkdir(m1_plot_dir)
#
#
# for pat in pats:
#     pat_fovs = harmonized_metadata.loc[(harmonized_metadata.Patient_ID == pat) & (harmonized_metadata.MIBI_data_generated.values) & (harmonized_metadata.Timepoint == 'on_nivo'), 'fov'].unique()
#     pat_df = cell_table_subset.loc[cell_table_subset.fov.isin(pat_fovs), :]
#
#     pat_dir = os.path.join(m1_plot_dir, 'patient_{}'.format(pat))
#     if not os.path.exists(pat_dir):
#         os.mkdir(pat_dir)
#
#     cohort_cluster_plot(
#         fovs=pat_fovs,
#         seg_dir=seg_dir,
#         save_dir=pat_dir,
#         cell_data=pat_df,
#         erode=True,
#         fov_col=settings.FOV_ID,
#         label_col=settings.CELL_LABEL,
#         cluster_col='M1_plot',
#         seg_suffix="_whole_cell.tiff",
#         cmap=m1_colormap,
#         display_fig=False,
#     )
#

# create crops for selected FOVs
fovs = ['TONIC_TMA6_R7C6', 'TONIC_TMA11_R7C4', 'TONIC_TMA11_R4C2', 'TONIC_TMA20_R2C3'] # patient 33, 62, 60, 114

subset_dir = os.path.join(plot_dir, 'Figure3f_CD68_overlays')
if not os.path.exists(subset_dir):
    os.mkdir(subset_dir)

cohort_cluster_plot(
    fovs=fovs,
    seg_dir=seg_dir,
    save_dir=subset_dir,
    cell_data=cell_table_func,
    erode=True,
    fov_col=settings.FOV_ID,
    label_col=settings.CELL_LABEL,
    cluster_col='M1_plot',
    seg_suffix="_whole_cell.tiff",
    cmap=m1_colormap,
    display_fig=False,
)


# select crops for visualization
fov1 = fovs[0]
row_start, col_start = 1448, 1300
row_len, col_len = 600, 600

fov1_image = io.imread(os.path.join(subset_dir, 'cluster_masks_colored', fov1 + '.tiff'))
fov1_image = fov1_image[row_start:row_start + row_len, col_start:col_start + col_len, :]
io.imsave(os.path.join(subset_dir, 'cluster_masks_colored', fov1 + '_crop.tiff'), fov1_image)


fov2 = fovs[1]
row_start, col_start = 900, 600
row_len, col_len = 600, 600

fov2_image = io.imread(os.path.join(subset_dir, 'cluster_masks_colored', fov2 + '.tiff'))
fov2_image = fov2_image[row_start:row_start + row_len, col_start:col_start + col_len, :]
io.imsave(os.path.join(subset_dir, 'cluster_masks_colored', fov2 + '_crop.tiff'), fov2_image)


fov3 = fovs[2]
row_start, col_start = 600, 1100
row_len, col_len = 600, 600

fov3_image = io.imread(os.path.join(subset_dir, 'cluster_masks_colored', fov3 + '.tiff'))
fov3_image = fov3_image[row_start:row_start + row_len, col_start:col_start + col_len, :]
io.imsave(os.path.join(subset_dir, 'cluster_masks_colored', fov3 + '_crop.tiff'), fov3_image)


fov4 = fovs[3]
row_start, col_start = 800, 0
row_len, col_len = 600, 600

fov4_image = io.imread(os.path.join(subset_dir, 'cluster_masks_colored', fov4 + '.tiff'))
fov4_image = fov4_image[row_start:row_start + row_len, col_start:col_start + col_len, :]
io.imsave(os.path.join(subset_dir, 'cluster_masks_colored', fov4 + '_crop.tiff'), fov4_image)


# diversity of cancer border on nivo
feature_name = 'cell_cluster_broad_diversity__cancer_border'
timepoint = 'on_nivo'

plot_df = combined_df.loc[(combined_df.feature_name_unique == feature_name) &
                                    (combined_df.Timepoint == timepoint), :]

fig, ax = plt.subplots(1, 1, figsize=(2, 4))
sns.stripplot(data=plot_df, x='Clinical_benefit', y='raw_mean', order=['Yes', 'No'],
                color='black', ax=ax)
sns.boxplot(data=plot_df, x='Clinical_benefit', y='raw_mean', order=['Yes', 'No'],
                color='grey', ax=ax, showfliers=False, width=0.3)
ax.set_title(feature_name + ' ' + timepoint)
ax.set_ylim([0, 2])
sns.despine()
plt.tight_layout()
plt.savefig(os.path.join(plot_dir, 'Figure3g_feature_{}_{}.pdf'.format(feature_name, timepoint)))
plt.close()



# corresponding overlays
cell_table_clusters = pd.read_csv(os.path.join(base_dir, 'analysis_files/cell_table_clusters.csv'))
annotations_by_mask = pd.read_csv(os.path.join(base_dir, 'intermediate_files/mask_dir', 'cell_annotation_mask.csv'))
annotations_by_mask = annotations_by_mask.rename(columns={'mask_name': 'tumor_region'})
cell_table_clusters = cell_table_clusters.merge(annotations_by_mask, on=['fov', 'label'], how='left')

# add column for cells in cancer border
cell_table_clusters['border_plot'] = cell_table_clusters.cell_cluster_broad
cell_table_clusters.loc[cell_table_clusters.tumor_region != 'cancer_border', 'border_plot'] = 'Other_region'

diversity_colormap = pd.DataFrame({'border_plot': ['Cancer', 'Structural', 'Mono_Mac', 'T', 'Other', 'Granulocyte', 'NK', 'B', 'Other_region'],
                             'color': ['white', 'darksalmon', 'red', 'yellow',  'yellowgreen', 'aqua', 'dodgerblue', 'darkviolet', 'dimgrey']})


# subset = plot_df.loc[plot_df.raw_mean > 0.8, :]
# #
# pats = [33, 59, 62, 65, 100, 115]
# pats = [7, 20, 50, 82, 106, 107, 112, 127]
# fovs = harmonized_metadata.loc[(harmonized_metadata.Patient_ID.isin(pats) & harmonized_metadata.MIBI_data_generated.values), 'fov'].unique()
#
# # 33, 62 previously included
#
# figure_dir = os.path.join(plot_dir, 'Figure3_border_diversity_test')
# if not os.path.exists(figure_dir):
#     os.mkdir(figure_dir)
#
# for pat in pats:
#     pat_dir = os.path.join(figure_dir, 'patient_{}'.format(pat))
#     if not os.path.exists(pat_dir):
#         os.mkdir(pat_dir)
#     pat_fovs = harmonized_metadata.loc[(harmonized_metadata.Patient_ID == pat) & (harmonized_metadata.MIBI_data_generated.values) & (harmonized_metadata.Timepoint == 'on_nivo'), 'fov'].unique()
#     pat_df = cell_table_clusters.loc[cell_table_clusters.fov.isin(pat_fovs), :]
#
#     cohort_cluster_plot(
#         fovs=pat_fovs,
#         seg_dir=seg_dir,
#         save_dir=pat_dir,
#         cell_data=pat_df,
#         erode=True,
#         fov_col=settings.FOV_ID,
#         label_col=settings.CELL_LABEL,
#         cluster_col='border_plot',
#         seg_suffix="_whole_cell.tiff",
#         cmap=diversity_colormap,
#         display_fig=False,
#     )

fovs = ['TONIC_TMA6_R7C6', 'TONIC_TMA14_R11C4'] # 33, 82

subset_dir = os.path.join(plot_dir, 'Figure3g_border_diversity')
if not os.path.exists(subset_dir):
    os.mkdir(subset_dir)

cohort_cluster_plot(
    fovs=fovs,
    seg_dir=seg_dir,
    save_dir=subset_dir,
    cell_data=cell_table_clusters,
    erode=True,
    fov_col=settings.FOV_ID,
    label_col=settings.CELL_LABEL,
    cluster_col='border_plot',
    seg_suffix="_whole_cell.tiff",
    cmap=diversity_colormap,
    display_fig=False,
)


# same thing for compartment masks
compartment_colormap = pd.DataFrame({'tumor_region': ['cancer_core', 'cancer_border', 'stroma_border', 'stroma_core', 'immune_agg'],
                         'color': ['blue', 'deepskyblue', 'lightcoral', 'firebrick', 'firebrick']})
subset_mask_dir = os.path.join(plot_dir, 'Figure3g_border_diversity_masks')
if not os.path.exists(subset_mask_dir):
    os.mkdir(subset_mask_dir)

cohort_cluster_plot(
    fovs=fovs,
    seg_dir=seg_dir,
    save_dir=subset_mask_dir,
    cell_data=cell_table_clusters,
    erode=True,
    fov_col=settings.FOV_ID,
    label_col=settings.CELL_LABEL,
    cluster_col='tumor_region',
    seg_suffix="_whole_cell.tiff",
    cmap=compartment_colormap,
    display_fig=False,
)

# crop overlays
fov1 = fovs[0]
row_start, col_start = 1200, 900
row_len, col_len = 600, 600

# scale bars: 100 um = 2048 pixels / 8 = 256 pixels, 100um = 256 / 600 = 0.426666666666666 of image

for dir in [subset_dir, subset_mask_dir]:
    fov1_image = io.imread(os.path.join(dir, 'cluster_masks_colored', fov1 + '.tiff'))
    fov1_image = fov1_image[row_start:row_start + row_len, col_start:col_start + col_len, :]
    io.imsave(os.path.join(dir, 'cluster_masks_colored', fov1 + '_crop.tiff'), fov1_image)

fov2 = fovs[1]
row_start, col_start = 100, 400
row_len, col_len = 600, 600

for dir in [subset_dir, subset_mask_dir]:
    fov2_image = io.imread(os.path.join(dir, 'cluster_masks_colored', fov2 + '.tiff'))
    fov2_image = fov2_image[row_start:row_start + row_len, col_start:col_start + col_len, :]
    io.imsave(os.path.join(dir, 'cluster_masks_colored', fov2 + '_crop.tiff'), fov2_image)