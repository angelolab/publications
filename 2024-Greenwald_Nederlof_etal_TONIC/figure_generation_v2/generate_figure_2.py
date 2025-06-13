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


base_dir = '/Volumes/Shared/Noah Greenwald/TONIC_Cohort'
metadata_dir = os.path.join(base_dir, 'intermediate_files/metadata')
image_dir = os.path.join(base_dir, 'image_data/samples/')
seg_dir = os.path.join(base_dir, 'segmentation_data/deepcell_output')
plot_dir = os.path.join(base_dir, 'figures')

harmonized_metadata = pd.read_csv(os.path.join(metadata_dir, 'harmonized_metadata.csv'))
study_fovs = harmonized_metadata.loc[harmonized_metadata.Timepoint.isin(['primary', 'baseline', 'pre_nivo', 'on_nivo']), 'fov'].values


# tumor compartment overlays
annotations_by_mask = pd.read_csv(os.path.join(base_dir, 'intermediate_files/mask_dir', 'cell_annotation_mask.csv'))

compartment_colormap = pd.DataFrame({'mask_name': ['cancer_core', 'cancer_border', 'stroma_border', 'stroma_core'],
                         'color': ['blue', 'deepskyblue', 'lightcoral', 'firebrick']})

compartment_plot_dir = os.path.join(plot_dir, 'Figure2a_compartment_overlays')
if not os.path.exists(compartment_plot_dir):
    os.mkdir(compartment_plot_dir)

# generate sampling of overlays across range of compartment frequencies to pick images
# # set thresholds
# fov_features = pd.read_csv(os.path.join(base_dir, 'analysis_files/fov_features_filtered.csv'))
# fov_features = fov_features.loc[fov_features.fov.isin(study_fovs), :]
# compartment_features = fov_features.loc[(fov_features.feature_name == 'cancer_core__proportion'), :]
#
# sns.histplot(data=compartment_features, x='raw_value',  bins=20, multiple='stack')
#
# # set thresholds for each group
# thresholds = {'low': [0.01, 0.1], 'mid': [0.3, .5], 'high': [0.6, 1]}
#
#
# # generate overlays

#
# for group in thresholds.keys():
#     group_dir = os.path.join(compartment_plot_dir, group)
#     if not os.path.exists(group_dir):
#         os.mkdir(group_dir)
#
#     min_val, max_val = thresholds[group]
#
#     # get fovs
#     fovs = compartment_features.loc[(compartment_features.raw_value > min_val) &
#                                         (compartment_features.raw_value < max_val), 'fov'].values
#
#     fovs = fovs[:10]
#     # generate overlays
#     cell_table_subset = annotations_by_mask.loc[(annotations_by_mask.fov.isin(fovs)), :]
#
#     cohort_cluster_plot(
#         fovs=fovs,
#         seg_dir=seg_dir,
#         save_dir=group_dir,
#         cell_data=annotations_by_mask,
#         erode=True,
#         fov_col=settings.FOV_ID,
#         label_col=settings.CELL_LABEL,
#         cluster_col='mask_name',
#         seg_suffix="_whole_cell.tiff",
#         cmap=compartment_colormap,
#         display_fig=False,
#     )

# scale bars: 800um image, 1/8 = 100um
# generate only selected overlays for figure
selected_compartment_fovs = ['TONIC_TMA3_R11C6', 'TONIC_TMA2_R7C3']
cell_table_subset = annotations_by_mask.loc[(annotations_by_mask.fov.isin(selected_compartment_fovs)), :]

cohort_cluster_plot(
    fovs=selected_compartment_fovs,
    seg_dir=seg_dir,
    save_dir=compartment_plot_dir,
    cell_data=cell_table_subset,
    erode=True,
    fov_col=settings.FOV_ID,
    label_col=settings.CELL_LABEL,
    cluster_col='mask_name',
    seg_suffix="_whole_cell.tiff",
    cmap=compartment_colormap,
    display_fig=False,
)

# cell proliferation
cell_table_func = pd.read_csv(os.path.join(base_dir, 'analysis_files/cell_table_func_single_positive.csv'))

cell_table_subset = cell_table_func.loc[(cell_table_func.fov.isin(study_fovs)), :]

ki67_colormap = pd.DataFrame({'proliferating': ['True', 'False'],
                            'color': ['red', 'dimgrey']})

cell_table_subset['proliferating'] = cell_table_subset.Ki67.astype(str)

ki67_plot_dir = os.path.join(plot_dir, 'Figure2a_ki67_overlays')
if not os.path.exists(ki67_plot_dir):
    os.mkdir(ki67_plot_dir)

# generate sampling of Ki67 proliferation levels to pick images

# # proliferating tumor cells
# ki67_features = fov_features.loc[(fov_features.feature_name == 'Ki67+__all') & (fov_features.compartment == 'all'), :]
#
# sns.histplot(data=ki67_features, x='raw_value',  bins=20, multiple='stack')
#
# # set thresholds for each group
# thresholds = {'low': [0.01, 0.2], 'mid': [0.2, .4], 'high': [0.5, 1]}

# for group in thresholds.keys():
#     group_dir = os.path.join(ki67_plot_dir, group)
#     if not os.path.exists(group_dir):
#         os.mkdir(group_dir)
#
#     min_val, max_val = thresholds[group]
#
#     # get fovs
#     fovs = ki67_features.loc[(ki67_features.raw_value > min_val) &
#                                         (ki67_features.raw_value < max_val), 'fov'].values
#
#     fovs = fovs[:10]
#
#     # generate overlays
#     cohort_cluster_plot(
#         fovs=fovs,
#         seg_dir=seg_dir,
#         save_dir=group_dir,
#         cell_data=cell_table_subset,
#         erode=True,
#         fov_col=settings.FOV_ID,
#         label_col=settings.CELL_LABEL,
#         cluster_col='proliferating',
#         seg_suffix="_whole_cell.tiff",
#         cmap=ki67_colormap,
#         display_fig=False,
#     )

# generate only selected overlays for figure
selected_ki67_fovs = ['TONIC_TMA14_R10C1', 'TONIC_TMA10_R1C5']
cohort_cluster_plot(
    fovs=selected_ki67_fovs,
    seg_dir=seg_dir,
    save_dir=ki67_plot_dir,
    cell_data=cell_table_subset,
    erode=True,
    fov_col=settings.FOV_ID,
    label_col=settings.CELL_LABEL,
    cluster_col='proliferating',
    seg_suffix="_whole_cell.tiff",
    cmap=ki67_colormap,
    display_fig=False,
)

# CD8T density
cell_table_clusters = pd.read_csv(os.path.join(base_dir, 'analysis_files/cell_table_clusters.csv'))
cell_table_clusters['CD8T_plot'] = cell_table_clusters.cell_cluster
cell_table_clusters.loc[cell_table_clusters.cell_cluster != 'CD8T', 'CD8T_plot'] = 'Other'

# set up plotting
CD8_colormap = pd.DataFrame({'CD8T_plot': ['CD8T', 'Other'],
                            'color': ['yellow', 'dimgrey']})

CD8_plot_dir = os.path.join(plot_dir, 'Figure2a_CD8_overlays')
if not os.path.exists(CD8_plot_dir):
    os.mkdir(CD8_plot_dir)

# generate sampling of CD8T levels to pick images
# CD8T_features = fov_features.loc[(fov_features.feature_name == 'CD8T__cluster_density') & (fov_features.compartment == 'all'), :]
#
# sns.histplot(data=CD8T_features, x='raw_value',  bins=20, multiple='stack')
#
#
# # set thresholds for each group
# thresholds = {'low': [0.001, 0.05], 'mid': [0.1, 0.2], 'high': [0.3, 0.6]}
#
# for group in thresholds.keys():
#     group_dir = os.path.join(CD8_plot_dir, group)
#     if not os.path.exists(group_dir):
#         os.mkdir(group_dir)
#
#     min_val, max_val = thresholds[group]
#
#     # get fovs
#     fovs = CD8T_features.loc[(CD8T_features.raw_value > min_val) &
#                              (CD8T_features.raw_value < max_val), 'fov'].values
#
#     fovs = fovs[:10]
#     # generate overlays
#     cell_table_subset = cell_table_clusters.loc[(cell_table_clusters.fov.isin(fovs)), :]
#
#     cohort_cluster_plot(
#         fovs=fovs,
#         seg_dir=seg_dir,
#         #save_dir=group_dir,
#         save_dir=CD8_plot_dir,
#         cell_data=cell_table_subset,
#         erode=True,
#         fov_col=settings.FOV_ID,
#         label_col=settings.CELL_LABEL,
#         cluster_col='CD8T_plot',
#         seg_suffix="_whole_cell.tiff",
#         cmap=CD8_colormap,
#         display_fig=False,
#     )

# generate only selected overlays for figure
selected_cd8_fovs = ['TONIC_TMA18_R4C6', 'TONIC_TMA10_R5C2']
cell_table_subset = cell_table_clusters.loc[(cell_table_clusters.fov.isin(selected_cd8_fovs)), :]

cohort_cluster_plot(
    fovs=selected_cd8_fovs,
    seg_dir=seg_dir,
    save_dir=CD8_plot_dir,
    cell_data=cell_table_subset,
    erode=True,
    fov_col=settings.FOV_ID,
    label_col=settings.CELL_LABEL,
    cluster_col='CD8T_plot',
    seg_suffix="_whole_cell.tiff",
    cmap=CD8_colormap,
    display_fig=False,
)

# Image diversity
diversity_plot_dir = os.path.join(plot_dir, 'Figure2a_diversity_overlays')
if not os.path.exists(diversity_plot_dir):
    os.mkdir(diversity_plot_dir)

diversity_colormap = pd.DataFrame({'cell_cluster_broad': ['Cancer', 'Structural', 'Mono_Mac', 'T',
                                                          'Other', 'Granulocyte', 'NK', 'B'],
                                   'color': ['dimgrey', 'darksalmon', 'red', 'yellow',
                                             'yellowgreen',
                                             'aqua', 'dodgerblue', 'darkviolet']})

# generate sampling of diversity levels to pick images
# diversity_features = fov_features.loc[(fov_features.feature_name == 'cluster_broad_diversity') &
#                                       (fov_features.compartment == 'all'), :]
#
# sns.histplot(data=diversity_features, x='raw_value',  bins=20, multiple='stack')
#
# # set thresholds for each group
# thresholds = {'low': [0.1, 0.5], 'mid': [1, 1.5], 'high': [2, 3]}
#
# for group in thresholds.keys():
#     group_dir = os.path.join(diversity_plot_dir, group)
#     if not os.path.exists(group_dir):
#         os.mkdir(group_dir)
#
#     min_val, max_val = thresholds[group]
#
#     # get fovs
#     fovs = diversity_features.loc[(diversity_features.raw_value > min_val) &
#                                         (diversity_features.raw_value < max_val), 'fov'].values
#
#     fovs = fovs[:10]
#     # generate overlays
#     cell_table_subset = cell_table_clusters.loc[(cell_table_clusters.fov.isin(fovs)), :]
#
#     cohort_cluster_plot(
#         fovs=fovs,
#         seg_dir=seg_dir,
#         save_dir=group_dir,
#         cell_data=cell_table_subset,
#         erode=True,
#         fov_col=settings.FOV_ID,
#         label_col=settings.CELL_LABEL,
#         cluster_col='cell_cluster_broad',
#         seg_suffix="_whole_cell.tiff",
#         cmap=diversity_colormap,
#         display_fig=False,
#     )

# generate only selected overlays for figure
selected_diversity_fovs = ['TONIC_TMA10_R5C6', 'TONIC_TMA15_R8C2']
cell_table_subset = cell_table_clusters.loc[(cell_table_clusters.fov.isin(selected_diversity_fovs)), :]

cohort_cluster_plot(
    fovs=selected_diversity_fovs,
    seg_dir=seg_dir,
    save_dir=diversity_plot_dir,
    cell_data=cell_table_subset,
    erode=True,
    fov_col=settings.FOV_ID,
    label_col=settings.CELL_LABEL,
    cluster_col='cell_cluster_broad',
    seg_suffix="_whole_cell.tiff",
    cmap=diversity_colormap,
    display_fig=False,
)

# summary plots for features
feature_metadata = pd.read_csv(os.path.join(base_dir, 'analysis_files/feature_metadata.csv'))
feature_classes = {'cell_abundance': ['density', 'density_ratio', 'density_proportion'],
                     'diversity': ['cell_diversity', 'region_diversity'],
                     'cell_phenotype': ['functional_marker', 'morphology', ],
                     'cell_interactions': ['mixing_score', 'linear_distance', 'kmeans_cluster'],
                   'structure': ['compartment_area_ratio', 'compartment_area', 'ecm_cluster', 'ecm_fraction', 'pixie_ecm', 'fiber']}

# label with appropriate high-level summary category
for feature_class in feature_classes.keys():
    feature_metadata.loc[feature_metadata.feature_type.isin(feature_classes[feature_class]), 'feature_class'] = feature_class

# add extra column to make stacked bar plotting work easily
feature_metadata_stacked = feature_metadata.copy()
feature_metadata_stacked['count'] = 1
feature_metadata_stacked = feature_metadata_stacked[['feature_class', 'count']].groupby(['feature_class']).sum().reset_index()
feature_metadata_stacked['feature'] = 'feature1'
feature_metadata_stacked['feature2'] = 'feature2'

feature_metadata_wide = pd.pivot(feature_metadata_stacked, index='feature', columns='feature_class', values='count')
feature_metadata_wide = feature_metadata_wide[['cell_phenotype', 'cell_abundance',  'structure', 'diversity', 'cell_interactions']]
feature_metadata_wide.plot.bar(stacked=True)
sns.despine()
plt.ylim(0, 1000)
plt.yticks(np.arange(0, 1001, 200))
plt.savefig(os.path.join(plot_dir, 'Figure2b_feature_class_counts_stacked.pdf'))
plt.close()

# cell types
feature_metadata = pd.read_csv(os.path.join(base_dir, 'analysis_files/feature_metadata.csv'))
feature_metadata['cell_class'] = np.nan
immune_cell_list = ['Mono_Mac'] + ['CD68_Mac', 'CD163_Mac', 'Mac_Other', 'Monocyte', 'APC'] + ['B'] + ['T'] + \
                   ['CD4T', 'CD8T', 'Treg', 'T_Other'] + ['Granulocyte'] + ['Neutrophil', 'Mast'] + ['NK']
cancer_cell_list = ['Cancer'] + ['Cancer_1', 'Cancer_2', 'Cancer_3']
structural_cell_list = ['Structural'] + ['Endothelium', 'CAF', 'Fibroblast', 'Smooth_Muscle']

for immune_cell in immune_cell_list:
    feature_metadata.loc[np.logical_and(feature_metadata.feature_type.isin(
        ['density', 'density_proportion', 'cell_diversity', 'region_diversity', 'morphology', 'functional_marker']),
        feature_metadata.feature_name.str.contains(immune_cell)), 'cell_class'] = 'Immune'

for cancer_cell in cancer_cell_list:
    feature_metadata.loc[np.logical_and(feature_metadata.feature_type.isin(
        ['density', 'density_proportion', 'cell_diversity', 'region_diversity', 'morphology', 'functional_marker']),
        feature_metadata.feature_name.str.contains(cancer_cell)), 'cell_class'] = 'Cancer'

for structural_cell in structural_cell_list:
    feature_metadata.loc[np.logical_and(feature_metadata.feature_type.isin(
        ['density', 'density_proportion', 'cell_diversity', 'region_diversity', 'morphology', 'functional_marker']),
        feature_metadata.feature_name.str.contains(structural_cell)), 'cell_class'] = 'Structural'

feature_metadata.loc[feature_metadata.feature_type.isin(
    ['density_ratio', 'linear_distance', 'kmeans_cluster', 'mixing_score', 'compartment_area', 'compartment_area_ratio']),
                     'cell_class'] = 'Multiple'

feature_metadata.loc[feature_metadata.feature_type.isn(['fiber', 'pixie_ecm', 'ecm_fraction', 'ecm_cluster']), 'cell_class'] = 'ECM'

feature_metadata.loc[feature_metadata.feature_name.str.contains('Immune_Other'), 'cell_class'] = 'Multiple'
feature_metadata.loc[feature_metadata.feature_name.str.contains('Other'), 'cell_class'] = 'Multiple'
feature_metadata.loc[np.logical_and(feature_metadata.feature_name.str.contains('all'),
                                    feature_metadata.cell_class.isna()), 'cell_class'] = 'Multiple'

# stacked bar version
feature_metadata_stacked = feature_metadata.copy()
feature_metadata_stacked['count'] = 1
feature_metadata_stacked = feature_metadata_stacked[['cell_class', 'count']].groupby(['cell_class']).sum().reset_index()
feature_metadata_stacked['feature'] = 'feature1'
feature_metadata_stacked['feature2'] = 'feature2'

feature_metadata_wide = pd.pivot(feature_metadata_stacked, index='feature', columns='cell_class', values='count')

# set order for stacked bar
feature_metadata_wide = feature_metadata_wide[['Immune', 'Multiple', 'Cancer', 'Structural', 'ECM']]
feature_metadata_wide.plot.bar(stacked=True)
sns.despine()
plt.ylim(0, 1000)
plt.yticks(np.arange(0, 1001, 200))

plt.savefig(os.path.join(plot_dir, 'Figure2b_cell_class_counts_stacked.pdf'))
plt.close()


# cluster features together to identify modules
fov_data_df = pd.read_csv(os.path.join(base_dir, 'analysis_files/combined_feature_data_filtered.csv'))
fov_data_df = fov_data_df.loc[fov_data_df.fov.isin(study_fovs), :]

# create wide df
fov_data_wide = fov_data_df.pivot(index='fov', columns='feature_name_unique', values='normalized_value')
corr_df = fov_data_wide.corr(method='spearman')
corr_df = corr_df.fillna(0)


clustergrid = sns.clustermap(corr_df, cmap='vlag', vmin=-1, vmax=1, figsize=(20, 20))
clustergrid.savefig(os.path.join(plot_dir, 'Figure2c_feature_clustermap_filtered.pdf'), dpi=300)
plt.close()

# # get names of features from clustergrid for annotating features within clusters
# feature_names = clustergrid.data2d.columns
# 
# start_idx = 700
# end_idx = 880
# clustergrid_small = sns.clustermap(corr_df.loc[feature_names[start_idx:end_idx], feature_names[start_idx:end_idx]], cmap='vlag', vmin=-1, vmax=1, figsize=(20, 20),
#                                    col_cluster=False, row_cluster=False)
# clustergrid_small.savefig(os.path.join(plot_dir, 'spearman_correlation_dp_functional_markers_clustermap_small_700.png'), dpi=300)
# plt.close()