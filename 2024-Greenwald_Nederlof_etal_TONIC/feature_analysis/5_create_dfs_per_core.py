import os

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.stats import spearmanr

from python_files.utils import create_long_df_by_functional, create_long_df_by_cluster

#
# This file creates plotting-ready data structures to enumerate the frequency, count, and density
# of cell populations. It also creates data structures with the frequency and count of functional
# marker positivity. Each of these dfs can be created across multiple levels of clustering
# granularity. For example, a broad classification might include T, B, Myeloid, Structural, and Cancer,
# whereas a more granular clustering scheme would separate out CD4T, CD8T, Tregs, etc.
#

base_dir = '/Volumes/Shared/Noah Greenwald/TONIC_Cohort'
intermediate_dir = os.path.join(base_dir, 'intermediate_files')
output_dir = os.path.join(base_dir, 'output_files')
analysis_dir = os.path.join(base_dir, 'analysis_files')

TIMEPOINT_NAMES = ['primary', 'baseline', 'pre_nivo', 'on_nivo']

if not os.path.exists(output_dir):
    os.makedirs(output_dir)

#
# Preprocess metadata to ensure all samples are present
#

study_name = 'TONIC'

# load relevant tables
core_metadata = pd.read_csv(os.path.join(intermediate_dir, 'metadata', f'{study_name}_data_per_core.csv'))
timepoint_metadata = pd.read_csv(os.path.join(intermediate_dir, 'metadata', f'{study_name}_data_per_timepoint.csv'))
harmonized_metadata = pd.read_csv(os.path.join(analysis_dir, 'harmonized_metadata.csv'))
cell_table_clusters = pd.read_csv(os.path.join(analysis_dir, 'cell_table_clusters.csv'))
cell_table_func = pd.read_csv(os.path.join(analysis_dir, 'cell_table_func_all.csv'))
cell_table_morph = pd.read_csv(os.path.join(analysis_dir, 'cell_table_morph.csv'))
cell_table_diversity = pd.read_csv(os.path.join(intermediate_dir, 'spatial_analysis/cell_neighbor_analysis/neighborhood_diversity_radius50.csv'))
cell_table_distances_broad = pd.read_csv(os.path.join(intermediate_dir, 'spatial_analysis/cell_neighbor_analysis/cell_cluster_broad_avg_dists-nearest_1.csv'))
area_df = pd.read_csv(os.path.join(intermediate_dir, 'mask_dir', 'fov_annotation_mask_area.csv'))
annotations_by_mask = pd.read_csv(os.path.join(intermediate_dir, 'mask_dir', 'cell_annotation_mask.csv'))
fiber_stats = pd.read_csv(os.path.join(intermediate_dir, 'fiber_segmentation_processed_data', 'fiber_stats_table.csv'))
fiber_tile_df = pd.read_csv(os.path.join(intermediate_dir, 'fiber_segmentation_processed_data/tile_stats_512', 'fiber_stats_table-tile_512.csv'))
kmeans_cell_table = pd.read_csv(os.path.join(intermediate_dir, 'spatial_analysis/neighborhood_analysis/cell_cluster_radius100_frequency_12', 'cell_table_clusters.csv'))
mixing_scores = pd.read_csv(os.path.join(intermediate_dir, 'spatial_analysis/mixing_score/cell_cluster_broad/homogeneous_mixing_scores.csv'))

# merge cell-level annotations
harmonized_annotations = annotations_by_mask
harmonized_annotations = harmonized_annotations.rename(columns={'mask_name': 'tumor_region'})
assert len(harmonized_annotations) == len(cell_table_clusters)

cell_table_clusters = cell_table_clusters.merge(harmonized_annotations, on=['fov', 'label'], how='left')
cell_table_func = cell_table_func.merge(harmonized_annotations, on=['fov', 'label'], how='left')
cell_table_morph = cell_table_morph.merge(harmonized_annotations, on=['fov', 'label'], how='left')
cell_table_diversity = cell_table_diversity.merge(harmonized_annotations, on=['fov', 'label'], how='left')
cell_table_distances_broad = cell_table_distances_broad.merge(harmonized_annotations, on=['fov', 'label'], how='left')

# check for FOVs present in imaged data that aren't in core metadata
missing_fovs = cell_table_clusters.loc[~cell_table_clusters.fov.isin(core_metadata.fov), 'fov'].unique()

# TONIC_TMA24_R11C2 and TONIC_TMA24_R11C3 are wrong tissue
cell_table_clusters = cell_table_clusters.loc[~cell_table_clusters.fov.isin(missing_fovs), :]
cell_table_func = cell_table_func.loc[~cell_table_func.fov.isin(missing_fovs), :]
cell_table_morph = cell_table_morph.loc[~cell_table_morph.fov.isin(missing_fovs), :]
cell_table_diversity = cell_table_diversity.loc[~cell_table_diversity.fov.isin(missing_fovs), :]
cell_table_distances_broad = cell_table_distances_broad.loc[~cell_table_distances_broad.fov.isin(missing_fovs), :]
fiber_stats = fiber_stats.loc[~fiber_stats.fov.isin(missing_fovs), :]
fiber_tile_df = fiber_tile_df.loc[~fiber_tile_df.fov.isin(missing_fovs), :]

#
# Generate counts and proportions of cell clusters per FOV
#

# Specify the types of cluster dfs to produce. Each row corresponds to a different way of
# summarizing the data. The first item in each list is the name for the df, the second is the
# name of the column to use to for the cluster labels, and the third is a boolean flag controlling
# whether normalized frequencies or total counts will be returned

cluster_df_params = [['cluster_broad_freq', 'cell_cluster_broad', True],
                     ['cluster_broad_count', 'cell_cluster_broad', False],
                     ['cluster_freq', 'cell_cluster', True],
                     ['cluster_count', 'cell_cluster', False],
                     ['meta_cluster_freq', 'cell_meta_cluster', True],
                     ['meta_cluster_count', 'cell_meta_cluster', False]]

cluster_dfs = []
for result_name, cluster_col_name, normalize in cluster_df_params:
    cluster_dfs.append(create_long_df_by_cluster(cell_table=cell_table_clusters,
                                                 result_name=result_name,
                                                 cluster_col_name=cluster_col_name,
                                                 subset_col='tumor_region',
                                                 normalize=normalize))


# create masks for dfs looking at only a subset of cells

# proportion of T cell subsets
tcell_mask = cell_table_clusters['cell_cluster'].isin(['Treg', 'CD8T', 'CD4T', 'T_Other'])

# proportion of immune cell subsets
immune_mask = cell_table_clusters['cell_cluster_broad'].isin(['Mono_Mac', 'T',
                                                              'Granulocyte', 'NK', 'B'])
immune_mask_2 = cell_table_clusters.cell_cluster == 'Immune_Other'
immune_mask = np.logical_or(immune_mask, immune_mask_2)

# proportion of structural subsets
structural_mask = cell_table_clusters['cell_cluster_broad'].isin(['Structural'])

# proportion of cancer subsets
cancer_mask = cell_table_clusters['cell_cluster_broad'].isin(['Cancer'])

cluster_mask_params = [['tcell_freq', 'cell_cluster', True, tcell_mask],
                       ['immune_freq', 'cell_cluster', True, immune_mask],
                       ['structural_freq', 'cell_meta_cluster', True, structural_mask],
                       ['cancer_freq', 'cell_meta_cluster', True, cancer_mask]]

for result_name, cluster_col_name, normalize, mask in cluster_mask_params:
    cluster_dfs.append(create_long_df_by_cluster(cell_table=cell_table_clusters.loc[mask, :],
                                                 result_name=result_name,
                                                 cluster_col_name=cluster_col_name,
                                                 subset_col='tumor_region',
                                                 normalize=normalize))

# calculate total number of cells per image
grouped_cell_counts = cell_table_clusters[['fov']].groupby('fov').value_counts()
grouped_cell_counts = pd.DataFrame(grouped_cell_counts)
grouped_cell_counts.columns = ['value']
grouped_cell_counts.reset_index(inplace=True)
grouped_cell_counts['metric'] = 'total_cell_count'
grouped_cell_counts['cell_type'] = 'all'
grouped_cell_counts['subset'] = 'all'

#

# calculate total number of cells per region per image
grouped_cell_counts_region = cell_table_clusters[['fov', 'tumor_region']].groupby(['fov', 'tumor_region']).value_counts()
grouped_cell_counts_region = pd.DataFrame(grouped_cell_counts_region)
grouped_cell_counts_region.columns = ['value']
grouped_cell_counts_region.reset_index(inplace=True)
grouped_cell_counts_region['metric'] = 'total_cell_count'
grouped_cell_counts_region.rename(columns={'tumor_region': 'subset'}, inplace=True)
grouped_cell_counts_region['cell_type'] = 'all'

# calculate proportions of cells per region per image
grouped_cell_freq_region = cell_table_clusters[['fov', 'tumor_region']].groupby(['fov'])
grouped_cell_freq_region = grouped_cell_freq_region['tumor_region'].value_counts(normalize=True)
grouped_cell_freq_region = pd.DataFrame(grouped_cell_freq_region)
grouped_cell_freq_region.columns = ['value']
grouped_cell_freq_region.reset_index(inplace=True)
grouped_cell_freq_region['metric'] = 'total_cell_freq'
grouped_cell_freq_region.rename(columns={'tumor_region': 'subset'}, inplace=True)
grouped_cell_freq_region['cell_type'] = 'all'

# add manually defined dfs to overall list
cluster_dfs.extend([grouped_cell_counts,
                    grouped_cell_counts_region,
                    grouped_cell_freq_region])

# create single df with appropriate metadata
total_df = pd.concat(cluster_dfs, axis=0)

# compute density of cells for counts-based metrics
count_metrics = total_df.metric.unique()
count_metrics = [x for x in count_metrics if 'count' in x]

count_df = total_df.loc[total_df.metric.isin(count_metrics), :]
area_df = area_df.rename(columns={'compartment': 'subset'})
count_df = count_df.merge(area_df, on=['fov', 'subset'], how='left')
count_df['value'] = count_df['value'] / count_df['area']
count_df['value'] = count_df['value'] * 1000

# rename metric from count to density
count_df['metric'] = count_df['metric'].str.replace('count', 'density')
count_df = count_df.drop(columns=['area'])
total_df = pd.concat([total_df, count_df], axis=0)

# check that all metadata from core_metadata succesfully transferred over
len_total_df = len(total_df)
total_df = total_df.merge(harmonized_metadata, on='fov', how='left')
assert len_total_df == len(total_df)


# save annotated cluster counts
total_df.to_csv(os.path.join(output_dir, 'cluster_df_per_core.csv'), index=False)

# create version aggregated by timepoint
total_df_grouped = total_df.groupby(['Tissue_ID', 'cell_type', 'metric', 'subset'])
total_df_timepoint = total_df_grouped['value'].agg([np.mean, np.std])
total_df_timepoint.reset_index(inplace=True)
total_df_timepoint = total_df_timepoint.merge(harmonized_metadata.drop(['fov', 'MIBI_data_generated'], axis=1).drop_duplicates(), on='Tissue_ID')

# save timepoint df
total_df_timepoint.to_csv(os.path.join(output_dir, 'cluster_df_per_timepoint.csv'), index=False)


#
# Create summary dataframe with proportions and counts of different functional marker populations
#


# Columns which are not thresholded (such as ratios between markers) can only be calculated for
# dfs looking at normalized expression, and need to be dropped when calculating counts
count_drop_cols = ['H3K9ac_H3K27me3_ratio', 'CD45RO_CD45RB_ratio']

# Create list to hold parameters for each df that will be produced
func_df_params = [['cluster_broad_count', 'cell_cluster_broad', False],
                  ['cluster_broad_freq', 'cell_cluster_broad', True],
                  ['cluster_count', 'cell_cluster', False],
                  ['cluster_freq', 'cell_cluster', True],
                  ['meta_cluster_count', 'cell_meta_cluster', False],
                  ['meta_cluster_freq', 'cell_meta_cluster', True]]


func_dfs = []
for result_name, cluster_col_name, normalize in func_df_params:
    # columns which are not functional markers need to be dropped from the df
    drop_cols = ['label']
    if not normalize:
        drop_cols.extend(count_drop_cols)

    # remove cluster_names except for the one specified for the df
    cluster_names = ['cell_meta_cluster', 'cell_cluster', 'cell_cluster_broad'] # , 'kmeans_labels']
    cluster_names.remove(cluster_col_name)
    drop_cols.extend(cluster_names)

    # create df
    func_dfs.append(create_long_df_by_functional(func_table=cell_table_func,
                                                 result_name=result_name,
                                                 cluster_col_name=cluster_col_name,
                                                 drop_cols=drop_cols,
                                                 normalize=normalize,
                                                 subset_col='tumor_region'))

# create combined df
total_df_func = pd.concat(func_dfs, axis=0)
total_df_func = total_df_func.merge(harmonized_metadata, on='fov', how='inner')

# save combined df
total_df_func.to_csv(os.path.join(output_dir, 'functional_df_per_core.csv'), index=False)

# create version aggregated by timepoint
total_df_grouped_func = total_df_func.groupby(['Tissue_ID', 'cell_type', 'functional_marker', 'metric', 'subset'])
total_df_timepoint_func = total_df_grouped_func['value'].agg([np.mean, np.std])
total_df_timepoint_func.reset_index(inplace=True)
total_df_timepoint_func = total_df_timepoint_func.merge(harmonized_metadata.drop(['fov', 'MIBI_data_generated'], axis=1).drop_duplicates(), on='Tissue_ID')

# save timepoint df
total_df_timepoint_func.to_csv(os.path.join(output_dir, 'functional_df_per_timepoint.csv'), index=False)

#
# Filter functional markers
#

# filter functional markers to only include FOVs with at least the specified number of cells
total_df = pd.read_csv(os.path.join(output_dir, 'cluster_df_per_core.csv'))
min_cells = 5

filtered_dfs = []
metrics = [['cluster_broad_count', 'cluster_broad_freq'],
           ['cluster_count', 'cluster_freq'],
           ['meta_cluster_count', 'meta_cluster_freq']]

for metric in metrics:
    # subset count df to include cells at the relevant clustering resolution
    for compartment in ['cancer_core', 'cancer_border', 'stroma_core', 'stroma_border',
                        'immune_agg', 'all']:
        count_df = total_df[total_df.metric == metric[0]]
        count_df = count_df[count_df.subset == compartment]

        # subset functional df to only include functional markers at this resolution
        func_df = total_df_func[total_df_func.metric.isin(metric)]
        func_df = func_df[func_df.subset == compartment]

        # for each cell type, determine which FOVs have high enough counts to be included
        for cell_type in count_df.cell_type.unique():
            keep_df = count_df[count_df.cell_type == cell_type]
            keep_df = keep_df[keep_df.value >= min_cells]
            keep_fovs = keep_df.fov.unique()

            # subset functional df to only include FOVs with high enough counts
            keep_markers = func_df[func_df.cell_type == cell_type]
            keep_markers = keep_markers[keep_markers.fov.isin(keep_fovs)]

            # append to list of filtered dfs
            filtered_dfs.append(keep_markers)

filtered_func_df = pd.concat(filtered_dfs)

# take subset for plotting average functional marker expression
sp_markers = [x for x in filtered_func_df.functional_marker.unique() if '__' not in x]
filtered_func_df_plot = filtered_func_df.loc[filtered_func_df.subset == 'all', :]
filtered_func_df_plot = filtered_func_df_plot.loc[filtered_func_df_plot.metric.isin(['cluster_broad_freq', 'cluster_freq', 'meta_cluster_freq']), :]
filtered_func_df_plot = filtered_func_df_plot.loc[filtered_func_df_plot.functional_marker.isin(sp_markers), :]

# save filtered df
filtered_func_df_plot.to_csv(os.path.join(output_dir, 'functional_df_per_core_filtered_plot.csv'), index=False)


#
# The code below is used to identify which functional markers should be evaluated in which cells
# based on frequency of positivity. Once these thresholds are set, the code doesn't need to be
# rerun again, and instead you can skip to the end and load the output file
#

# # identify combinations of markers and cell types to include in analysis based on threshold
# mean_percent_positive = 0.05
# broad_df = filtered_func_df[filtered_func_df.metric == 'cluster_broad_freq']
# broad_df = broad_df[broad_df.functional_marker.isin(sp_markers)]
# broad_df = broad_df[broad_df.subset == 'all']
# broad_df = broad_df[broad_df.Timepoint.isin(TIMEPOINT_NAMES)]
# broad_df_agg = broad_df[['functional_marker', 'cell_type', 'value']].groupby(['cell_type', 'functional_marker']).agg(np.mean)
#
# broad_df_agg = broad_df_agg.reset_index()
# broad_df = broad_df_agg.pivot(index='cell_type', columns='functional_marker', values='value')
# broad_df_include = broad_df > mean_percent_positive
#
# # include for all cells
# general_markers = ['Ki67', 'HLA1', 'Vim', 'H3K9ac_H3K27me3_ratio']
# broad_df_include[[general_markers]] = True
#
# # CD45 isoform ratios
# double_pos = np.logical_and(broad_df_include['CD45RO'], broad_df_include['CD45RB'])
# broad_df_include['CD45RO_CD45RB_ratio'] = double_pos
#
# # Cancer expression
# broad_df_include.loc['Cancer', ['HLADR', 'CD57']] = True
#
# broad_df_include.to_csv(os.path.join(intermediate_dir, 'post_processing', 'inclusion_matrix_broad.csv'))
#
# # apply thresholds to medium level clustering
# assignment_dict_med = {'Cancer': ['Cancer_1', 'Cancer_2', 'Cancer_3'],
#                      'Mono_Mac': ['CD68_Mac', 'CD163_Mac', 'Mac_Other', 'Monocyte', 'APC'],
#                      'B': ['B'],
#                      'T': ['CD4T', 'CD8T', 'Treg', 'T_Other', 'Immune_Other'],
#                      'Granulocyte': ['Neutrophil', 'Mast'],
#                      'Structural': ['Endothelium', 'CAF', 'Fibroblast', 'Smooth_Muscle'],
#                      'NK': ['NK'],
#                      'Other': ['Other']}
#
# # get a list of all cell types
# cell_types = np.concatenate([assignment_dict_med[key] for key in assignment_dict_med.keys()])
# cell_types.sort()
#
# med_df_include = pd.DataFrame(index=cell_types, columns=broad_df.columns)
#
# for key in assignment_dict_med.keys():
#     values = assignment_dict_med[key]
#     med_df_include.loc[values] = broad_df_include.loc[key].values
#
# # check if assignment makes sense
# med_df = filtered_func_df[filtered_func_df.metric == 'cluster_freq']
# med_df = med_df[med_df.functional_marker.isin(sp_markers)]
# med_df = med_df[med_df.subset == 'all']
# med_df = med_df[med_df.Timepoint.isin(TIMEPOINT_NAMES)]
# med_df_agg = med_df[['functional_marker', 'cell_type', 'value']].groupby(['cell_type', 'functional_marker']).agg(np.mean)
#
# med_df_agg.reset_index(inplace=True)
# med_df = med_df_agg.pivot(index='cell_type', columns='functional_marker', values='value')
# med_df_bin = med_df > mean_percent_positive
#
# # add CD38 signal
# med_df_include.loc[['Endothelium', 'Immune_Other'], 'CD38'] = True
#
# # add IDO signal
# med_df_include.loc[['APC'], 'IDO'] = True
#
# # compare to see where assignments disagree, to see if any others need to be added
# new_includes = (med_df_bin == True) & (med_df_include == False)
#
# med_df_include.to_csv(os.path.join(intermediate_dir, 'post_processing', 'inclusion_matrix_med.csv'))
#
# # do the same for the fine-grained clustering
# assignment_dict_meta = {'Cancer_1': ['Cancer_CD56', 'Cancer_CK17', 'Cancer_Ecad'],
#                    'Cancer_2': ['Cancer_SMA', 'Cancer_Vim'],
#                    'Cancer_3': ['Cancer_Other', 'Cancer_Mono'],
#                    'CD68_Mac': ['CD68'],
#                    'CD163_Mac': ['CD163'],
#                    'Mac_Other': ['CD68_CD163_DP'],
#                    'Monocyte': ['CD4_Mono', 'CD14'],
#                    'APC': ['CD11c_HLADR'],
#                    'B':  ['CD20'],
#                    'Endothelium': ['CD31', 'CD31_VIM'],
#                    'CAF': ['FAP', 'FAP_SMA'],
#                    'Fibroblast': ['Stroma_Collagen', 'Stroma_Fibronectin', 'VIM'],
#                    'Smooth_Muscle': ['SMA'],
#                    'NK': ['CD56'],
#                    'Neutrophil': ['Neutrophil'],
#                    'Mast': ['Mast'],
#                    'CD4T': ['CD4T', 'CD4T_HLADR'],
#                    'CD8T': ['CD8T'],
#                    'Treg': ['Treg'],
#                    'T_Other': ['CD3_DN','CD4T_CD8T_DP'],
#                    'Immune_Other': ['Immune_Other'],
#                    'Other': ['Other']}
#
# # get a list of all cell types
# cell_types = np.concatenate([assignment_dict_meta[key] for key in assignment_dict_meta.keys()])
# cell_types.sort()
#
# meta_df_include = pd.DataFrame(index=cell_types, columns=broad_df.columns)
#
# for key in assignment_dict_meta.keys():
#     values = assignment_dict_meta[key]
#     meta_df_include.loc[values] = med_df_include.loc[key].values
#
# # check if assignment makes sense
# meta_df = filtered_func_df[filtered_func_df.metric == 'meta_cluster_freq']
# meta_df = meta_df[meta_df.functional_marker.isin(sp_markers)]
# meta_df = meta_df[meta_df.subset == 'all']
# meta_df = meta_df[meta_df.Timepoint.isin(TIMEPOINT_NAMES)]
# meta_df_agg = meta_df[['functional_marker', 'cell_type', 'value']].groupby(['cell_type', 'functional_marker']).agg(np.mean)
#
# meta_df_agg.reset_index(inplace=True)
# meta_df = meta_df_agg.pivot(index='cell_type', columns='functional_marker', values='value')
# meta_df_bin = meta_df > mean_percent_positive
#
# # compare to see where assignments disagree
# new_includes = (meta_df_bin == True) & (meta_df_include == False)
#
# meta_df_include.to_csv(os.path.join(intermediate_dir, 'post_processing', 'inclusion_matrix_meta.csv'))

# process functional data so that only the specified cell type/marker combos are included

# load matrices
broad_df_include = pd.read_csv(os.path.join(intermediate_dir, 'post_processing', 'inclusion_matrix_broad.csv'), index_col=0)
med_df_include = pd.read_csv(os.path.join(intermediate_dir, 'post_processing', 'inclusion_matrix_med.csv'), index_col=0)
meta_df_include = pd.read_csv(os.path.join(intermediate_dir, 'post_processing', 'inclusion_matrix_meta.csv'), index_col=0)

# identify metrics and dfs that will be filtered
filtering = [['cluster_broad_count', 'cluster_broad_freq', broad_df_include],
           ['cluster_count', 'cluster_freq', med_df_include],
           ['meta_cluster_count', 'meta_cluster_freq', meta_df_include]]

combo_dfs = []

for filters in filtering:
    # get variables
    metric_names = filters[:2]
    metric_df = filters[2]

    # subset functional df to only include functional markers at this resolution
    func_df = filtered_func_df[filtered_func_df.metric.isin(metric_names)]

    # loop over each cell type, and get the corresponding markers
    for cell_type in metric_df.index:
        markers = metric_df.columns[metric_df.loc[cell_type] == True]

        # subset functional df to only include this cell type
        func_df_cell = func_df[func_df.cell_type == cell_type]

        # subset functional df to only include markers for this cell type
        func_df_cell = func_df_cell[func_df_cell.functional_marker.isin(markers)]

        # append to list of dfs
        combo_dfs.append(func_df_cell)

#
# The commented code below is the same as the commented code above, but generates thresholds
# for double positive functional markers. Once run, you can skip to the end and load in the
# generated output files
#

# # identify combinations of markers and cell types to include in analysis based on threshold
# dp_markers = [x for x in filtered_func_df.functional_marker.unique() if '__' in x]
# mean_percent_positive_dp = 0.05
# broad_df_dp = filtered_func_df[filtered_func_df.metric == 'cluster_broad_freq']
# broad_df_dp = broad_df_dp[broad_df_dp.subset == 'all']
# broad_df_dp = broad_df_dp[broad_df_dp.Timepoint.isin(TIMEPOINT_NAMES)]
# broad_df_dp = broad_df_dp[broad_df_dp.functional_marker.isin(dp_markers)]
# broad_df_dp_agg = broad_df_dp[['functional_marker', 'cell_type', 'value']].groupby(['cell_type', 'functional_marker']).agg(np.mean)
#
# broad_df_dp_agg = broad_df_dp_agg.reset_index()
# broad_df_dp = broad_df_dp_agg.pivot(index='cell_type', columns='functional_marker', values='value')
# broad_df_dp_include = broad_df_dp > mean_percent_positive_dp
#
# broad_df_dp_include.to_csv(os.path.join(intermediate_dir, 'post_processing', 'inclusion_matrix_broad_dp.csv'))
#
# # do the same for medium-level clustering
# med_df_dp = filtered_func_df[filtered_func_df.metric == 'cluster_freq']
# med_df_dp = med_df_dp[med_df_dp.subset == 'all']
# med_df_dp = med_df_dp[med_df_dp.Timepoint.isin(TIMEPOINT_NAMES)]
# med_df_dp = med_df_dp[med_df_dp.functional_marker.isin(dp_markers)]
# med_df_dp_agg = med_df_dp[['functional_marker', 'cell_type', 'value']].groupby(['cell_type', 'functional_marker']).agg(np.mean)
#
# # create matrix of cell types and markers
# med_df_dp_agg = med_df_dp_agg.reset_index()
# med_df_dp = med_df_dp_agg.pivot(index='cell_type', columns='functional_marker', values='value')
# med_df_dp_include = med_df_dp > mean_percent_positive_dp
#
# med_df_dp_include.to_csv(os.path.join(intermediate_dir, 'post_processing', 'inclusion_matrix_med_dp.csv'))
#
#
# # do the same for finest-level clustering
# meta_df_dp = filtered_func_df[filtered_func_df.metric == 'meta_cluster_freq']
# meta_df_dp = meta_df_dp[meta_df_dp.subset == 'all']
# meta_df_dp = meta_df_dp[meta_df_dp.Timepoint.isin(TIMEPOINT_NAMES)]
# meta_df_dp = meta_df_dp[meta_df_dp.functional_marker.isin(dp_markers)]
# meta_df_dp_agg = meta_df_dp[['functional_marker', 'cell_type', 'value']].groupby(['cell_type', 'functional_marker']).agg(np.mean)
#
# # create matrix of cell types and markers
# meta_df_dp_agg = meta_df_dp_agg.reset_index()
# meta_df_dp = meta_df_dp_agg.pivot(index='cell_type', columns='functional_marker', values='value')
# meta_df_dp_include = meta_df_dp > mean_percent_positive_dp
#
# meta_df_dp_include.to_csv(os.path.join(intermediate_dir, 'post_processing', 'inclusion_matrix_meta_dp.csv'))

# load inclusion matrices
broad_df_include_dp = pd.read_csv(os.path.join(intermediate_dir, 'post_processing', 'inclusion_matrix_broad_dp.csv'), index_col=0)
med_df_include_dp = pd.read_csv(os.path.join(intermediate_dir, 'post_processing', 'inclusion_matrix_med_dp.csv'), index_col=0)
meta_df_include_dp = pd.read_csv(os.path.join(intermediate_dir, 'post_processing', 'inclusion_matrix_meta_dp.csv'), index_col=0)

# identify metrics and dfs that will be filtered
filtering = [['cluster_broad_count', 'cluster_broad_freq', broad_df_include_dp],
           ['cluster_count', 'cluster_freq', med_df_include_dp],
           ['meta_cluster_count', 'meta_cluster_freq', meta_df_include_dp]]

for filters in filtering:
    # get variables
    metric_names = filters[:2]
    metric_df = filters[2]

    # subset functional df to only include functional markers at this resolution
    func_df = filtered_func_df[filtered_func_df.metric.isin(metric_names)]

    # loop over each cell type, and get the corresponding markers
    for cell_type in metric_df.index:
        markers = metric_df.columns[metric_df.loc[cell_type] == True]

        # subset functional df to only include this cell type
        func_df_cell = func_df[func_df.cell_type == cell_type]

        # subset functional df to only include markers for this cell type
        func_df_cell = func_df_cell[func_df_cell.functional_marker.isin(markers)]

        # append to list of dfs
        combo_dfs.append(func_df_cell)


# create manual df with total functional marker positivity across all cells in an image
# dp_markers = [x for x in filtered_func_df.functional_marker.unique() if '__' in x]
# func_table_small = cell_table_func.loc[:, ~cell_table_func.columns.isin(['cell_cluster', 'cell_cluster_broad', 'cell_meta_cluster', 'label', 'tumor_region'])]
# func_table_small = func_table_small.loc[:, ~func_table_small.columns.isin(dp_markers)]
#
# # group by specified columns
# grouped_table = func_table_small.groupby('fov')
# transformed = grouped_table.agg(np.mean)
# transformed.reset_index(inplace=True)
#
# # reshape to long df
# long_df = pd.melt(transformed, id_vars=['fov'], var_name='functional_marker')
# long_df['metric'] = 'total_freq'
# long_df['cell_type'] = 'all'
# long_df['subset'] = 'all'
#
# long_df = long_df.merge(harmonized_metadata, on='fov', how='inner')
#
# long_df.to_csv(os.path.join(output_dir, 'total_func_per_core.csv'), index=False)

# append to list of dfs
long_df = pd.read_csv(os.path.join(output_dir, 'total_func_per_core.csv'))
combo_dfs.append(long_df)

# combine
combo_df = pd.concat(combo_dfs)
combo_df.to_csv(os.path.join(output_dir, 'functional_df_per_core_filtered.csv'), index=False)

# create version of filtered df aggregated by timepoint
combo_df_grouped_func = combo_df.groupby(['Tissue_ID', 'cell_type', 'functional_marker', 'metric', 'subset'])
combo_df_timepoint_func = combo_df_grouped_func['value'].agg([np.mean, np.std])
combo_df_timepoint_func.reset_index(inplace=True)
combo_df_timepoint_func = combo_df_timepoint_func.merge(harmonized_metadata.drop(['fov', 'MIBI_data_generated'], axis=1).drop_duplicates(), on='Tissue_ID')

# save timepoint df
combo_df_timepoint_func.to_csv(os.path.join(output_dir, 'functional_df_per_timepoint_filtered.csv'), index=False)


#
# Remove double positive functional markers that are highly correlated with single positive scores
# This section can be skipped once it's been run once and just use the resulting output file
#

# cluster_resolution = [['cluster_broad_freq', 'cluster_broad_count'],
#            ['cluster_freq', 'cluster_count'],
#            ['meta_cluster_freq', 'meta_cluster_count']]
#
# all_markers = combo_df.functional_marker.unique()
# dp_markers = [x for x in all_markers if '__' in x]
#
# exclude_lists = []
# for cluster in cluster_resolution:
#
#     # subset functional df to only include functional markers at this resolution
#     func_df = combo_df[combo_df.metric.isin(cluster)]
#
#     # add unique identifier for cell + marker combo
#     func_df['feature_name'] = func_df['cell_type'] + '__' + func_df['functional_marker']
#
#     # subset the df further to look at just frequency and just one compartment
#     func_df_subset = func_df[func_df.metric == cluster[0]]
#     func_df_subset = func_df_subset[func_df_subset.subset == 'all']
#
#     # loop over each cell type, and each double positive functional marker
#     exclude_markers = []
#     cell_types = func_df_subset.cell_type.unique()
#     for cell_type in cell_types:
#         for marker in dp_markers:
#             # get the two markers that make up the double positive marker
#             marker_1, marker_2 = marker.split('__')
#
#             # subset to only include this cell type and these markers
#             current_df = func_df_subset.loc[func_df_subset.cell_type == cell_type, :]
#             current_df = current_df.loc[current_df.functional_marker.isin([marker, marker_1, marker_2]), :]
#
#             # these markers are not present in this cell type
#             if len(current_df) == 0:
#                 continue
#
#             # this double positive marker is not present in this cell type
#             if marker not in current_df.functional_marker.unique():
#                 continue
#
#             current_df_wide = current_df.pivot(index='fov', columns='functional_marker', values='value')
#
#             # the double positive marker is present, but both single positives are not; exclude it
#             if len(current_df_wide.columns) != 3:
#                 exclude_markers.append(marker)
#                 continue
#
#             corr_1, _ = spearmanr(current_df_wide[marker_1].values, current_df_wide[marker].values)
#             corr_2, _ = spearmanr(current_df_wide[marker_2].values, current_df_wide[marker].values)
#
#             if (corr_1 > 0.7) | (corr_2 > 0.7):
#                 exclude_markers.append(cell_type + '__' + marker)
#
#     # add to list
#     exclude_lists.append(exclude_markers)
#
# # construct df to hold list of exlcuded cells
# exclude_df_cluster = ['cluster_broad_freq'] * len(exclude_lists[0]) + ['cluster_freq'] * len(exclude_lists[1]) + ['meta_cluster_freq'] * len(exclude_lists[2])
# exclude_df_name = exclude_lists[0] + exclude_lists[1] + exclude_lists[2]
#
# exclude_df = pd.DataFrame({'metric': exclude_df_cluster, 'feature_name': exclude_df_name})
# exclude_df.to_csv(os.path.join(intermediate_dir, 'post_processing', 'exclude_double_positive_markers.csv'), index=False)

# use previously generated exclude list
exclude_df = pd.read_csv(os.path.join(intermediate_dir, 'post_processing', 'exclude_double_positive_markers.csv'))

dedup_dfs = []

cluster_resolution = [['cluster_broad_freq', 'cluster_broad_count'],
           ['cluster_freq', 'cluster_count'],
           ['meta_cluster_freq', 'meta_cluster_count']]

for cluster in cluster_resolution:
    # subset functional df to only include functional markers at this resolution
    func_df = combo_df[combo_df.metric.isin(cluster)]

    # add unique identifier for cell + marker combo
    func_df['feature_name'] = func_df['cell_type'] + '__' + func_df['functional_marker']

    exclude_names = exclude_df.loc[exclude_df.metric == cluster[0], 'feature_name'].values

    # remove double positive markers that are highly correlated with single positive scores
    func_df = func_df[~func_df.feature_name.isin(exclude_names)]
    dedup_dfs.append(func_df)


dedup_dfs.append(long_df)
deduped_df = pd.concat(dedup_dfs)
deduped_df = deduped_df.drop('feature_name', axis=1)

# save deduped df
deduped_df.to_csv(os.path.join(output_dir, 'functional_df_per_core_filtered_deduped.csv'), index=False)

# create version aggregated by timepoint
deduped_df_grouped = deduped_df.groupby(['Tissue_ID', 'cell_type', 'functional_marker', 'metric', 'subset'])
deduped_df_timepoint = deduped_df_grouped['value'].agg([np.mean, np.std])
deduped_df_timepoint.reset_index(inplace=True)
deduped_df_timepoint = deduped_df_timepoint.merge(harmonized_metadata.drop(['fov', 'MIBI_data_generated'], axis=1).drop_duplicates(), on='Tissue_ID')

# save timepoint df
deduped_df_timepoint.to_csv(os.path.join(output_dir, 'functional_df_per_timepoint_filtered_deduped.csv'), index=False)


# morphology metric summary

# Create list to hold parameters for each df that will be produced
morph_df_params = [['cluster_broad_freq', 'cell_cluster_broad'],
                  ['cluster_freq', 'cell_cluster'],
                  ['meta_cluster_freq', 'cell_meta_cluster']]

morph_dfs = []
for result_name, cluster_col_name in morph_df_params:

    # remove cluster_names except for the one specified for the df
    drop_cols = ['cell_meta_cluster', 'cell_cluster', 'cell_cluster_broad', 'label']
    drop_cols.remove(cluster_col_name)

    # create df
    morph_dfs.append(create_long_df_by_functional(func_table=cell_table_morph,
                                                 result_name=result_name,
                                                 cluster_col_name=cluster_col_name,
                                                 drop_cols=drop_cols,
                                                 normalize=True,
                                                 subset_col='tumor_region'))

# create combined df
total_df_morph = pd.concat(morph_dfs, axis=0)
total_df_morph = total_df_morph.merge(harmonized_metadata, on='fov', how='inner')
total_df_morph = total_df_morph.rename(columns={'functional_marker': 'morphology_feature'})

# save df
total_df_morph.to_csv(os.path.join(output_dir, 'morph_df_per_core.csv'), index=False)


# create manual df with total morphology marker average across all cells in an image
morph_table_small = cell_table_morph.loc[:, ~cell_table_morph.columns.isin(['cell_cluster', 'cell_cluster_broad', 'cell_meta_cluster', 'label', 'tumor_region'])]

# group by specified columns
grouped_table = morph_table_small.groupby('fov')
transformed = grouped_table.agg(np.mean)
transformed.reset_index(inplace=True)

# reshape to long df
long_df = pd.melt(transformed, id_vars=['fov'], var_name='morphology_feature')
long_df['metric'] = 'total_freq'
long_df['cell_type'] = 'all'
long_df['subset'] = 'all'

long_df = long_df.merge(harmonized_metadata, on='fov', how='inner')

long_df.to_csv(os.path.join(output_dir, 'total_morph_per_core.csv'), index=False)

# filter morphology markers to only include FOVs with at least the specified number of cells
total_df = pd.read_csv(os.path.join(output_dir, 'cluster_df_per_core.csv'))
min_cells = 5

filtered_dfs = []
metrics = [['cluster_broad_count', 'cluster_broad_freq'],
           ['cluster_count', 'cluster_freq'],
           ['meta_cluster_count', 'meta_cluster_freq']]
for metric in metrics:
    # subset count df to include cells at the relevant clustering resolution
    for compartment in ['cancer_core', 'cancer_border', 'stroma_core', 'stroma_border',
                        'immune_agg', 'all']:
        count_df = total_df[total_df.metric == metric[0]]
        count_df = count_df[count_df.subset == compartment]

        # subset morphology df to only include morphology metrics at this resolution
        morph_df = total_df_morph[total_df_morph.metric == metric[1]]
        morph_df = morph_df[morph_df.subset == compartment]

        # for each cell type, determine which FOVs have high enough counts to be included
        for cell_type in count_df.cell_type.unique():
            keep_df = count_df[count_df.cell_type == cell_type]
            keep_df = keep_df[keep_df.value >= min_cells]
            keep_fovs = keep_df.fov.unique()

            # subset morphology df to only include FOVs with high enough counts
            keep_features = morph_df[morph_df.cell_type == cell_type]
            keep_features = keep_features[keep_features.fov.isin(keep_fovs)]

            # append to list of filtered dfs
            filtered_dfs.append(keep_features)

filtered_dfs.append(long_df)
filtered_morph_df = pd.concat(filtered_dfs)

# save filtered df
filtered_morph_df.to_csv(os.path.join(output_dir, 'morph_df_per_core_filtered.csv'), index=False)

# create version aggregated by timepoint
filtered_morph_df_grouped = filtered_morph_df.groupby(['Tissue_ID', 'cell_type', 'morphology_feature', 'metric', 'subset'])
filtered_morph_df_timepoint = filtered_morph_df_grouped['value'].agg([np.mean, np.std])
filtered_morph_df_timepoint.reset_index(inplace=True)
filtered_morph_df_timepoint = filtered_morph_df_timepoint.merge(harmonized_metadata.drop(['fov', 'MIBI_data_generated'], axis=1).drop_duplicates(), on='Tissue_ID')

# save timepoint df
filtered_morph_df_timepoint.to_csv(os.path.join(output_dir, 'morph_df_per_timepoint_filtered.csv'), index=False)

# remove redundant morphology features
block1 = ['area', 'major_axis_length', 'minor_axis_length', 'perimeter', 'convex_area', 'equivalent_diameter']

block2 = ['area_nuclear', 'major_axis_length_nuclear', 'minor_axis_length_nuclear', 'perimeter_nuclear', 'convex_area_nuclear', 'equivalent_diameter_nuclear']

block3 = ['eccentricity', 'major_axis_equiv_diam_ratio']

block4 = ['eccentricity_nuclear', 'major_axis_equiv_diam_ratio_nuclear', 'perim_square_over_area_nuclear']

deduped_morph_df = filtered_morph_df.loc[~filtered_morph_df.morphology_feature.isin(block1[1:] + block2[1:] + block3[1:] + block4[1:]), :]

# only keep complex morphology features for cancer cells, remove everything except area and nc for others
cancer_clusters = ['Cancer_1', 'Cancer_2', 'Cancer_3', 'Cancer_CD56', 'Cancer_CK17',
                   'Cancer_Ecad', 'Cancer_Mono', 'Cancer_SMA', 'Cancer_Vim']
basic_morph_features = ['area', 'area_nuclear', 'nc_ratio']

deduped_morph_df = deduped_morph_df.loc[~(~(deduped_morph_df.cell_type.isin(cancer_clusters)) & ~(deduped_morph_df.morphology_feature.isin(basic_morph_features))), :]

# saved deduped
deduped_morph_df.to_csv(os.path.join(output_dir, 'morph_df_per_core_filtered_deduped.csv'), index=False)

# same for timepoints
deduped_morph_df_timepoint = filtered_morph_df_timepoint.loc[~filtered_morph_df_timepoint.morphology_feature.isin(block1[1:] + block2[1:] + block3[1:] + block4[1:]), :]
deduped_morph_df_timepoint = deduped_morph_df_timepoint.loc[~(~(deduped_morph_df_timepoint.cell_type.isin(cancer_clusters)) & ~(deduped_morph_df_timepoint.morphology_feature.isin(basic_morph_features))), :]
deduped_morph_df_timepoint.to_csv(os.path.join(output_dir, 'morph_df_per_timepoint_filtered_deduped.csv'), index=False)

#
# spatial features
#

# format mixing scores
cols = mixing_scores.columns.tolist()
keep_cols = [col for col in cols if 'mixing_score' in col]
mixing_scores = mixing_scores[['fov'] + keep_cols]

mixing_scores = pd.melt(mixing_scores, id_vars=['fov'], var_name='mixing_score', value_name='value')
mixing_scores.to_csv(os.path.join(output_dir, 'formatted_mixing_scores.csv'), index=False)

# compute local diversity scores per image

# Create list to hold parameters for each df that will be produced
diversity_df_params = [['cluster_broad_freq', 'cell_cluster_broad'],
                  ['cluster_freq', 'cell_cluster'],
                  ['meta_cluster_freq', 'cell_meta_cluster']]

diversity_dfs = []
for result_name, cluster_col_name in diversity_df_params:

    # remove cluster_names except for the one specified for the df
    drop_cols = ['cell_meta_cluster', 'cell_cluster', 'cell_cluster_broad', 'label']
    drop_cols.remove(cluster_col_name)

    # create df
    diversity_dfs.append(create_long_df_by_functional(func_table=cell_table_diversity,
                                                 result_name=result_name,
                                                 cluster_col_name=cluster_col_name,
                                                 drop_cols=drop_cols,
                                                 normalize=True,
                                                 subset_col='tumor_region'))

# create combined df
total_df_diversity = pd.concat(diversity_dfs, axis=0)
total_df_diversity = total_df_diversity.merge(harmonized_metadata, on='fov', how='inner')
total_df_diversity = total_df_diversity.rename(columns={'functional_marker': 'diversity_feature'})

# save df
total_df_diversity.to_csv(os.path.join(output_dir, 'diversity_df_per_core.csv'), index=False)

# filter diversity scores to only include FOVs with at least the specified number of cells
total_df = pd.read_csv(os.path.join(output_dir, 'cluster_df_per_core.csv'))
min_cells = 5

filtered_dfs = []
metrics = [['cluster_broad_count', 'cluster_broad_freq'],
           ['cluster_count', 'cluster_freq'],
           ['meta_cluster_count', 'meta_cluster_freq']]
for metric in metrics:
    # subset count df to include cells at the relevant clustering resolution
    for compartment in ['cancer_core', 'cancer_border', 'stroma_core',
                        'stroma_border', 'immune_agg', 'all']:
        count_df = total_df[total_df.metric == metric[0]]
        count_df = count_df[count_df.subset == compartment]

        # subset diversity df to only include diversity metrics at this resolution
        diversity_df = total_df_diversity[total_df_diversity.metric == metric[1]]
        diversity_df = diversity_df[diversity_df.subset == compartment]

        # for each cell type, determine which FOVs have high enough counts to be included
        for cell_type in count_df.cell_type.unique():
            keep_df = count_df[count_df.cell_type == cell_type]
            keep_df = keep_df[keep_df.value >= min_cells]
            keep_fovs = keep_df.fov.unique()

            # subset morphology df to only include FOVs with high enough counts
            keep_features = diversity_df[diversity_df.cell_type == cell_type]
            keep_features = keep_features[keep_features.fov.isin(keep_fovs)]

            # append to list of filtered dfs
            filtered_dfs.append(keep_features)

filtered_diversity_df = pd.concat(filtered_dfs)

# save filtered df
filtered_diversity_df.to_csv(os.path.join(output_dir, 'diversity_df_per_core_filtered.csv'), index=False)

# create version aggregated by timepoint
filtered_diversity_df_grouped = filtered_diversity_df.groupby(['Tissue_ID', 'cell_type', 'diversity_feature', 'metric', 'subset'])
filtered_diversity_df_timepoint = filtered_diversity_df_grouped['value'].agg([np.mean, np.std])
filtered_diversity_df_timepoint.reset_index(inplace=True)
filtered_diversity_df_timepoint = filtered_diversity_df_timepoint.merge(harmonized_metadata.drop(['fov', 'MIBI_data_generated'], axis=1).drop_duplicates(), on='Tissue_ID')

# save timepoint df
filtered_diversity_df_timepoint.to_csv(os.path.join(output_dir, 'diversity_df_per_timepoint_filtered.csv'), index=False)


# investigate correlation between diversity scores
fov_data = filtered_diversity_df.copy()
fov_data['feature_name_unique'] = fov_data['cell_type'] + '_' + fov_data['diversity_feature']
fov_data = fov_data.loc[(fov_data.subset == 'all') & (fov_data.metric == 'cluster_freq')]
fov_data = fov_data.loc[fov_data.diversity_feature != 'diversity_cell_meta_cluster']
fov_data_wide = fov_data.pivot(index='fov', columns='feature_name_unique', values='value')

corr_df = fov_data_wide.corr(method='spearman')

# replace Nans
corr_df = corr_df.fillna(0)
clustergrid = sns.clustermap(corr_df, cmap='vlag', vmin=-1, vmax=1, figsize=(20, 20))

# save deduped df that excludes cell meta cluster
deduped_diversity_df = filtered_diversity_df.loc[filtered_diversity_df.diversity_feature != 'diversity_cell_meta_cluster']
deduped_diversity_df.to_csv(os.path.join(output_dir, 'diversity_df_per_core_filtered_deduped.csv'), index=False)

# create version aggregated by timepoint
deduped_diversity_df_grouped = deduped_diversity_df.groupby(['Tissue_ID', 'cell_type', 'diversity_feature', 'metric', 'subset'])
deduped_diversity_df_timepoint = deduped_diversity_df_grouped['value'].agg([np.mean, np.std])
deduped_diversity_df_timepoint.reset_index(inplace=True)
deduped_diversity_df_timepoint = deduped_diversity_df_timepoint.merge(harmonized_metadata.drop(['fov', 'MIBI_data_generated'], axis=1).drop_duplicates(), on='Tissue_ID')

# save timepoint df
deduped_diversity_df_timepoint.to_csv(os.path.join(output_dir, 'diversity_df_per_timepoint_filtered_deduped.csv'), index=False)


# process linear distance dfs
total_df_distance = create_long_df_by_functional(func_table=cell_table_distances_broad,
                                             result_name='cluster_broad_freq',
                                             cluster_col_name='cell_cluster_broad',
                                             drop_cols=['label'],
                                             normalize=True,
                                             subset_col='tumor_region')

# create combined df
total_df_distance = total_df_distance.merge(harmonized_metadata, on='fov', how='inner')
total_df_distance = total_df_distance.rename(columns={'functional_marker': 'linear_distance'})
total_df_distance.dropna(inplace=True)

# save df
total_df_distance.to_csv(os.path.join(output_dir, 'distance_df_per_core.csv'), index=False)

# filter distance scores to only include FOVs with at least the specified number of cells
total_df = pd.read_csv(os.path.join(output_dir, 'cluster_df_per_core.csv'))
min_cells = 5

filtered_dfs = []
metrics = [['cluster_broad_count', 'cluster_broad_freq']]

for metric in metrics:
    # subset count df to include cells at the relevant clustering resolution
    for compartment in ['cancer_core', 'cancer_border', 'stroma_core',
                        'stroma_border', 'immune_agg', 'all']:
        count_df = total_df[total_df.metric == metric[0]]
        count_df = count_df[count_df.subset == compartment]

        # subset distance df to only include distance metrics at this resolution
        distance_df = total_df_distance[total_df_distance.metric == metric[1]]
        distance_df = distance_df[distance_df.subset == compartment]

        # for each cell type, determine which FOVs have high enough counts to be included
        for cell_type in count_df.cell_type.unique():
            keep_df = count_df[count_df.cell_type == cell_type]
            keep_df = keep_df[keep_df.value >= min_cells]
            keep_fovs = keep_df.fov.unique()

            # subset morphology df to only include FOVs with high enough counts
            keep_features = distance_df[distance_df.cell_type == cell_type]
            keep_features = keep_features[keep_features.fov.isin(keep_fovs)]

            # append to list of filtered dfs
            filtered_dfs.append(keep_features)

filtered_distance_df = pd.concat(filtered_dfs)

# save filtered df
filtered_distance_df.to_csv(os.path.join(output_dir, 'distance_df_per_core_filtered.csv'), index=False)

#
# Remove distances that are correlated with abundance of cell type. This code can be skipped after
# running the first time and the output file can be loaded instead
#

# load data
# total_df = pd.read_csv(os.path.join(output_dir, 'cluster_df_per_core.csv'))
# density_df = total_df.loc[(total_df.metric == 'cluster_broad_density') & (total_df.subset == 'all')]
# filtered_distance_df = filtered_distance_df.loc[(filtered_distance_df.metric == 'cluster_broad_freq') & (filtered_distance_df.subset == 'all')]
#
# # remove images without tumor cells
# density_df = density_df.loc[density_df.Timepoint != 'lymphnode_neg', :]
# filtered_distance_df = filtered_distance_df.loc[filtered_distance_df.fov.isin(density_df.fov.unique())]
# cell_types = filtered_distance_df.cell_type.unique()
#
# # calculate which pairings to keep
# keep_cells, keep_features = [], []
#
# for cell_type in cell_types:
#     density_subset = density_df.loc[density_df.cell_type == cell_type]
#     distance_subset = filtered_distance_df.loc[filtered_distance_df.linear_distance == cell_type]
#     distance_wide = distance_subset.pivot(index='fov', columns='cell_type', values='value')
#     distance_wide.reset_index(inplace=True)
#     distance_wide = pd.merge(distance_wide, density_subset[['fov', 'value']], on='fov', how='inner')
#
#     # get correlations
#     corr_df = distance_wide.corr(method='spearman', numeric_only=True)
#
#     # determine which features to keep
#     corr_vals = corr_df.loc['value', :].abs()
#     corr_vals = corr_vals[corr_vals < 0.7]
#
#     # add to list of features to keep
#     keep_cells.extend(corr_vals.index)
#     keep_features.extend([cell_type] * len(corr_vals.index))
#
# keep_df = pd.DataFrame({'cell_type': keep_cells, 'feature_name': keep_features})
#
# keep_df.to_csv(os.path.join(intermediate_dir, 'post_processing', 'distance_df_keep_features.csv'), index=False)

# filter distance df to only include features with low correlation with abundance
keep_df = pd.read_csv(os.path.join(intermediate_dir, 'post_processing', 'distance_df_keep_features.csv'))


deduped_dfs = []
for cell_type in keep_df.cell_type.unique():
    keep_features = keep_df.loc[keep_df.cell_type == cell_type, 'feature_name'].unique()
    if len(keep_features) > 0:
        keep_df_subset = filtered_distance_df.loc[filtered_distance_df.cell_type == cell_type]
        keep_df_subset = keep_df_subset.loc[keep_df_subset.linear_distance.isin(keep_features)]
        deduped_dfs.append(keep_df_subset)

deduped_distance_df = pd.concat(deduped_dfs)

# save filtered df
deduped_distance_df.to_csv(os.path.join(output_dir, 'distance_df_per_core_deduped.csv'), index=False)

# create version aggregated by timepoint
deduped_distance_df_grouped = deduped_distance_df.groupby(['Tissue_ID', 'cell_type', 'linear_distance', 'metric', 'subset'])
deduped_distance_df_timepoint = deduped_distance_df_grouped['value'].agg([np.mean, np.std])
deduped_distance_df_timepoint.reset_index(inplace=True)
deduped_distance_df_timepoint = deduped_distance_df_timepoint.merge(harmonized_metadata.drop(['fov', 'MIBI_data_generated'], axis=1).drop_duplicates(), on='Tissue_ID')

# save timepoint df
deduped_distance_df_timepoint.to_csv(os.path.join(output_dir, 'distance_df_per_timepoint_deduped.csv'), index=False)


# fiber analysis
fiber_stats.columns = fiber_stats.columns.str.replace('avg_', '')
fiber_stats.columns = fiber_stats.columns.str.replace('fiber_', '')
fiber_stats = fiber_stats.loc[:, ~fiber_stats.columns.isin(['label', 'centroid-0', 'centroid-1'])]

fiber_stats = fiber_stats.merge(harmonized_metadata[['Tissue_ID', 'fov']], on=['fov'], how='left')

fiber_df_long = pd.melt(fiber_stats, id_vars=['Tissue_ID', 'fov'], var_name='fiber_metric', value_name='value')
fiber_df_long['fiber_metric'] = 'fiber_' + fiber_df_long['fiber_metric']

fiber_df_long.to_csv(os.path.join(output_dir, 'fiber_df_per_core.csv'), index=False)

# create version aggregated by timepoint
fiber_df_grouped = fiber_df_long.groupby(['Tissue_ID', 'fiber_metric'])
fiber_df_timepoint = fiber_df_grouped['value'].agg([np.mean, np.std])
fiber_df_timepoint.reset_index(inplace=True)
fiber_df_timepoint = fiber_df_timepoint.merge(harmonized_metadata.drop(['fov', 'MIBI_data_generated'], axis=1).drop_duplicates(), on='Tissue_ID')

# save timepoint df
fiber_df_timepoint.to_csv(os.path.join(output_dir, 'fiber_df_per_timepoint.csv'), index=False)

# for tiles, get max per image
fiber_tile_df = fiber_tile_df.dropna()
fiber_tile_df = fiber_tile_df.loc[:, ~fiber_tile_df.columns.isin(['pixel_density', 'tile_y', 'tile_x'])]
fiber_tile_df.columns = fiber_tile_df.columns.str.replace('avg_', '')
fiber_tile_df.columns = fiber_tile_df.columns.str.replace('fiber_', '')
fiber_tile_df = fiber_tile_df.merge(harmonized_metadata[['Tissue_ID', 'fov']], on=['fov'], how='left')

# group by fov
fiber_tile_df_means = fiber_tile_df.groupby(['Tissue_ID', 'fov']).agg(np.max)
fiber_tile_df_means.reset_index(inplace=True)

fiber_tile_df_long = pd.melt(fiber_tile_df_means, id_vars=['Tissue_ID', 'fov'], var_name='fiber_metric', value_name='value')
fiber_tile_df_long['fiber_metric'] = 'max_fiber_' + fiber_tile_df_long['fiber_metric']

fiber_tile_df_long.to_csv(os.path.join(output_dir, 'fiber_df_per_tile.csv'), index=False)

# create version aggregated by timepoint
fiber_tile_df_grouped = fiber_tile_df_long.groupby(['Tissue_ID', 'fiber_metric'])
fiber_tile_df_timepoint = fiber_tile_df_grouped['value'].agg([np.mean, np.std])
fiber_tile_df_timepoint.reset_index(inplace=True)
fiber_tile_df_timepoint = fiber_tile_df_timepoint.merge(harmonized_metadata.drop(['fov', 'MIBI_data_generated'], axis=1).drop_duplicates(), on='Tissue_ID')

# save timepoint df
fiber_tile_df_timepoint.to_csv(os.path.join(output_dir, 'fiber_df_per_tile_timepoint.csv'), index=False)


# kmeans neighborhood proportions
# image level proportions
fov_cell_sum = kmeans_cell_table[['fov', 'kmeans_neighborhood']].groupby(by=['fov']).count().reset_index()
fov_cell_sum = fov_cell_sum.rename(columns={'kmeans_neighborhood': 'cells_in_image'})

# create df with all fovs and all kmeans rows
kmeans_cluster_num = len(np.unique(kmeans_cell_table.kmeans_neighborhood.dropna()))
all_fovs_df = []
for fov in np.unique(kmeans_cell_table.fov):
    df = pd.DataFrame({
        'fov': [fov] * kmeans_cluster_num,
        'kmeans_neighborhood': list(range(1, kmeans_cluster_num+1))
    })

    all_fovs_df.append(df)
all_fovs_df = pd.concat(all_fovs_df)

# get kmeans cluster counts per image, merge with all cluster df, replace nan with zero
cluster_prop = kmeans_cell_table[['fov', 'kmeans_neighborhood', 'label']].groupby(
    by=['fov', 'kmeans_neighborhood']).count().reset_index()

cluster_prop = all_fovs_df.merge(cluster_prop, on=['fov', 'kmeans_neighborhood'], how='left')
cluster_prop.fillna(0, inplace=True)

# calculate proportions
cluster_prop = cluster_prop.merge(fov_cell_sum, on=['fov'])
cluster_prop = cluster_prop.rename(columns={'label': 'cells_in_cluster'})
cluster_prop['proportion'] = cluster_prop.cells_in_cluster / cluster_prop.cells_in_image

cluster_prop.to_csv(os.path.join(output_dir, 'neighborhood_image_proportions.csv'), index=False)


# stroma and cancer compartment proportions
kmeans_cells = kmeans_cell_table[['fov', 'kmeans_neighborhood', 'label']]
compartment_data = annotations_by_mask.merge(kmeans_cells, on=['fov', 'label'])

all_compartments_df = []
for fov in np.unique(kmeans_cell_table.fov):
    df = pd.DataFrame({
        'fov': [fov] * 5 * kmeans_cluster_num,
        'mask_name': ['cancer_border'] * kmeans_cluster_num + ['cancer_core'] * kmeans_cluster_num +
        ['stroma_border'] * kmeans_cluster_num + ['stroma_core'] * kmeans_cluster_num + ['immune_agg'] * kmeans_cluster_num,
        'kmeans_neighborhood': list(range(1, kmeans_cluster_num+1)) * 5,
    })

    all_compartments_df.append(df)
all_compartments_df = pd.concat(all_compartments_df)

# get kmeans cluster counts per compartment in each image, merge with all cluster df, replace nan with zero
compartment_data = compartment_data.groupby(by=['fov', 'mask_name', 'kmeans_neighborhood']).count().reset_index()

all_data = all_compartments_df.merge(compartment_data, on=['fov', 'mask_name', 'kmeans_neighborhood'], how='left')
all_data.fillna(0, inplace=True)
all_data = all_data.rename(columns={'label': 'cells_in_cluster'})

# get compartment cell counts
compartment_cell_sum = all_data[['fov', 'mask_name', 'cells_in_cluster']].groupby(
    by=['fov', 'mask_name']).sum().reset_index()
compartment_cell_sum = compartment_cell_sum.rename(columns={'cells_in_cluster': 'total_cells'})

# calculate proportions
df = all_data.merge(compartment_cell_sum, on=['fov', 'mask_name'])
df['proportion'] = df.cells_in_cluster / df.total_cells
df = df.dropna().sort_values(by=['fov', 'mask_name'])

df.to_csv(os.path.join(output_dir, 'neighborhood_compartment_proportions.csv'), index=False)
