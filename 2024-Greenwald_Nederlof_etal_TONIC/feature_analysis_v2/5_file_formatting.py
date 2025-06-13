import os
import numpy as np
import pandas as pd
from functools import reduce

BASE_DIR = '/Volumes/Shared/Noah Greenwald/TONIC_Cohort'
ANALYSIS_DIR = os.path.join(BASE_DIR, 'analysis_files')
INTERMEDIATE_DIR = os.path.join(BASE_DIR, 'intermediate_files')
FORMATTED_DIR = os.path.join(INTERMEDIATE_DIR, 'formatted_files')

############# SKIP IF NO CHANGES TO INTERMEDIATE FILES / COMPARTMENT ASSIGNMENTS

### KMEANS FORMATTING ###
kmeans_dir = os.path.join(INTERMEDIATE_DIR, 'spatial_analysis/neighborhood_analysis/cell_cluster_radius100_frequency_12')
kmeans_cell_table = pd.read_csv(os.path.join(kmeans_dir, 'cell_table_clusters.csv'))
cell_table = pd.read_csv(os.path.join(ANALYSIS_DIR, 'combined_cell_table_normalized_cell_labels_updated.csv'))

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

# make wide df
output_df = cluster_prop.pivot(index='fov', columns='kmeans_neighborhood', values='proportion')
for i in output_df.columns:
    if i != 'fov':
        output_df = output_df.rename(columns={i: 'cluster_' + str(i) + '__proportion'})

output_df.reset_index(inplace=True)
output_df = output_df[output_df.fov.isin(cell_table.fov.unique())]
output_df.to_csv(os.path.join(FORMATTED_DIR, 'neighborhood_image_proportions.csv'), index=False)


# stroma and cancer compartment proportions
cell_table = cell_table.rename(columns={'compartment_example': 'compartment'})
annotations_by_mask = cell_table[['fov', 'label', 'compartment']]

kmeans_cells = kmeans_cell_table[['fov', 'kmeans_neighborhood', 'label']]
compartment_data = annotations_by_mask.merge(kmeans_cells, on=['fov', 'label'])

all_compartments_df = []
for fov in np.unique(kmeans_cell_table.fov):
    df = pd.DataFrame({
        'fov': [fov] * 5 * kmeans_cluster_num,
        'compartment': ['cancer_border'] * kmeans_cluster_num + ['cancer_core'] * kmeans_cluster_num +
        ['stroma_border'] * kmeans_cluster_num + ['stroma_core'] * kmeans_cluster_num + ['immune_agg'] * kmeans_cluster_num,
        'kmeans_neighborhood': list(range(1, kmeans_cluster_num+1)) * 5,
    })

    all_compartments_df.append(df)
all_compartments_df = pd.concat(all_compartments_df)

# get kmeans cluster counts per compartment in each image, merge with all cluster df, replace nan with zero
compartment_data = compartment_data.groupby(by=['fov', 'compartment', 'kmeans_neighborhood']).count().reset_index()

all_data = all_compartments_df.merge(compartment_data, on=['fov', 'compartment', 'kmeans_neighborhood'], how='left')
all_data.fillna(0, inplace=True)
all_data = all_data.rename(columns={'label': 'cells_in_cluster'})

# get compartment cell counts
compartment_cell_sum = all_data[['fov', 'compartment', 'cells_in_cluster']].groupby(
    by=['fov', 'compartment']).sum().reset_index()
compartment_cell_sum = compartment_cell_sum.rename(columns={'cells_in_cluster': 'total_cells'})

# calculate proportions
df = all_data.merge(compartment_cell_sum, on=['fov', 'compartment'])
df['proportion'] = df.cells_in_cluster / df.total_cells
df = df.dropna().sort_values(by=['fov', 'compartment'])

# make wide df
output_df = df.pivot(index='fov', columns=['compartment', 'kmeans_neighborhood'], values='proportion')
output_df.columns = output_df.columns.map(lambda index: f'cluster_{index[1]}__proportion__{index[0]}')
output_df.reset_index(inplace=True)
output_df.to_csv(os.path.join(FORMATTED_DIR, 'neighborhood_compartment_proportions.csv'), index=False)


### MIXING SCORE FORMATTING ###
mixing_score_dir = os.path.join(INTERMEDIATE_DIR, 'spatial_analysis/mixing_score')
mixing_score_df = pd.read_csv(os.path.join(mixing_score_dir, 'homogeneous_mixing_scores.csv'))

keep_cols = [col for col in mixing_score_df.columns if 'mixing_score' in col]
mixing_score_df_sub = mixing_score_df[['fov'] + keep_cols]
mixing_score_df_sub.to_csv(os.path.join(FORMATTED_DIR, 'formatted_mixing_scores.csv'), index=False)


### CELL DISTANCES FORMATTING ###
cell_dist_dir = os.path.join(INTERMEDIATE_DIR, 'spatial_analysis/cell_neighbor_analysis')
distances_df = pd.read_csv(os.path.join(cell_dist_dir, 'cell_cluster_broad_avg_dists-nearest_1.csv'))
distances_df = distances_df[distances_df.fov.isin(cell_table.fov.unique())].reset_index(drop=True)
keep_df = pd.read_csv(os.path.join(cell_dist_dir, 'distance_df_keep_features.csv'))

for col in sorted(distances_df.cell_cluster_broad.unique()):
    for cell_type in sorted(distances_df.cell_cluster_broad.unique()):
        df = keep_df[np.logical_and(keep_df.feature_name == col, keep_df.cell_type == cell_type)]
        if len(df) == 0:
            distances_df.loc[distances_df['cell_cluster_broad'] == cell_type, col] = np.nan

for col in distances_df.columns:
    if not col in ['fov', 'label', 'cell_cluster_broad']:
        distances_df = distances_df.rename(columns={col: f'distance_to__{col}'})
distances_df.to_csv(os.path.join(FORMATTED_DIR, 'cell_distances_filtered.csv'), index=False)


### FIBER STATS FORMATTING ###
fiber_stats = pd.read_csv(os.path.join(INTERMEDIATE_DIR, 'fiber_segmentation_processed_data/fiber_stats_table.csv'))
fiber_stats.columns = fiber_stats.columns.str.replace('avg_', '')
fiber_stats.columns = fiber_stats.columns.str.replace('fiber_', '')
for i in fiber_stats.columns:
    if i != 'fov':
        fiber_stats = fiber_stats.rename(columns={i: f'fiber_{i}'})
fiber_stats.to_csv(os.path.join(FORMATTED_DIR, 'fiber_stats_table.csv'), index=False)

# for tiles, get max per image
fiber_tile_df = pd.read_csv(os.path.join(INTERMEDIATE_DIR, 'fiber_segmentation_processed_data/fiber_stats_table-tile_512.csv'))
fiber_tile_df = fiber_tile_df.dropna()
fiber_tile_df = fiber_tile_df.loc[:, ~fiber_tile_df.columns.isin(['pixel_density', 'tile_y', 'tile_x'])]
fiber_tile_df.columns = fiber_tile_df.columns.str.replace('avg_', '')
fiber_tile_df.columns = fiber_tile_df.columns.str.replace('fiber_', '')

fiber_tile_df_max = fiber_tile_df.groupby(['fov']).agg('max')
fiber_tile_df_max.reset_index(inplace=True)

for i in fiber_tile_df_max.columns:
    if i != 'fov':
        fiber_tile_df_max = fiber_tile_df_max.rename(columns={i: f'max_fiber_{i}'})

fiber_tile_df_max.reset_index(inplace=True, drop=True)
fiber_tile_df_max.to_csv(os.path.join(FORMATTED_DIR, 'fiber_stats_per_tile.csv'), index=False)


### ECM STATS FORMATTING ###
# ecm clustering results
ecm_df = pd.read_csv(os.path.join(INTERMEDIATE_DIR, 'ecm/fov_cluster_counts.csv'))
ecm_frac = 1 - ecm_df.No_ECM.values
ecm_df['Cold_Coll_Norm'] = ecm_df.Cold_Coll.values / ecm_frac
ecm_df['Hot_Coll_Norm'] = ecm_df.Hot_Coll.values / ecm_frac
ecm_df.loc[ecm_frac == 0, ['Cold_Coll_Norm', 'Hot_Coll_Norm']] = 0

ecm_df = ecm_df[[col for col in ecm_df.columns if col != 'fov_cluster']]
for col in ecm_df.columns:
    if col != 'fov':
        ecm_df = ecm_df.rename(columns={col: col + '__proportion'})

# compute ECM fraction for each image
# ecm_mask_dir = os.path.join(intermediate_dir, 'ecm/masks')
# fovs = list_folders(ecm_mask_dir)
# fov_name, fov_frac = [], []
#
# for fov in fovs:
#     mask = io.imread(os.path.join(ecm_mask_dir, fov, 'total_ecm.tiff'))
#     fov_frac.append(mask.sum() / mask.size)
#     fov_name.append(fov)
#
# ecm_frac_df = pd.DataFrame({'fov': fov_name, 'ecm_frac': fov_frac})
# ecm_frac_df.to_csv(os.path.join(FORMATTED_DIR, 'ecm/ecm_fraction_fov.csv'), index=False)

# ecm fraction
ecm_frac_df = pd.read_csv(os.path.join(INTERMEDIATE_DIR, 'ecm/ecm_fraction_fov.csv'))
ecm_frac_df = ecm_frac_df.rename(columns={'ecm_frac': 'ecm_fraction'})

# ecm pixel cluster density
ecm_clusters = pd.read_csv(os.path.join(INTERMEDIATE_DIR, 'ecm_pixel_clustering/fov_pixel_cluster_counts.csv'))
area_df = pd.read_csv(os.path.join(INTERMEDIATE_DIR, 'mask_dir/fov_annotation_mask_area.csv'))
area_df = area_df.loc[area_df.compartment == 'all', ['fov', 'area']]
ecm_clusters_density = pd.merge(ecm_clusters, area_df, on='fov')
ecm_clusters_density['density'] = ecm_clusters_density['counts'] / ecm_clusters_density['area']

ecm_clusters_density = ecm_clusters_density.pivot(index='fov', columns=['pixel_meta_cluster_rename'], values='density')
ecm_clusters_density.reset_index(inplace=True)

for col in ecm_clusters_density.columns:
    if col != 'fov':
        ecm_clusters_density = ecm_clusters_density.rename(columns={col: 'Pixie__cluster__' + col + '__density'})

# ecm pixel cluster proportion
ecm_clusters_wide = pd.pivot(ecm_clusters, index='fov', columns='pixel_meta_cluster_rename', values='counts')
ecm_clusters_wide = ecm_clusters_wide.apply(lambda x: x / x.sum(), axis=1)
ecm_clusters_wide = ecm_clusters_wide.reset_index()

for col in ecm_clusters_wide.columns:
    if col != 'fov':
        ecm_clusters_wide = ecm_clusters_wide.rename(columns={col: 'Pixie__cluster__' + col + '__proportion'})

# ecm neighborhood density
ecm_neighborhoods = pd.read_csv(os.path.join(INTERMEDIATE_DIR, 'ecm_pixel_clustering/neighborhood/fov_neighborhood_counts.csv'))
ecm_neighborhoods = ecm_neighborhoods.rename(
    columns={'Cluster1': 'Collagen', 'Cluster2': 'SMA', 'Cluster3': 'Collagen_Vim', 'Cluster4': 'Vim_FAP',
             'Cluster5': 'Fibronectin'})
ecm_neighborhood_density = pd.melt(ecm_neighborhoods.iloc[:, :-1], id_vars='fov', var_name='ecm_neighborhood',
                                   value_name='counts')
ecm_neighborhood_density = pd.merge(ecm_neighborhood_density, area_df, on='fov')
ecm_neighborhood_density['density'] = ecm_neighborhood_density['counts'] / ecm_neighborhood_density['area']

ecm_neighborhood_density = pd.pivot(ecm_neighborhood_density, index='fov', columns='ecm_neighborhood', values='density')
ecm_neighborhood_density = ecm_neighborhood_density.reset_index()
for col in ecm_neighborhood_density.columns:
    if col != 'fov':
        ecm_neighborhood_density = ecm_neighborhood_density.rename(
            columns={col: 'Pixie__neighborhood__' + col + '__density'})

# ecm neighborhood proportion
ecm_neighborhoods_prop = ecm_neighborhoods.copy()
for col in ecm_neighborhoods_prop.columns[1:-1]:
    ecm_neighborhoods_prop[col] = ecm_neighborhoods_prop[col] / ecm_neighborhoods_prop['total']
    ecm_neighborhoods_prop = ecm_neighborhoods_prop.rename(
        columns={col: 'Pixie__neighborhood__' + col + '__proportion'})
ecm_neighborhoods_prop = ecm_neighborhoods_prop[[col for col in ecm_neighborhoods_prop.columns if col != 'total']]

# ecm shape axis ratio
ecm_object_ratio = pd.read_csv(os.path.join(INTERMEDIATE_DIR, 'ecm_pixel_clustering/shape_analysis/fov_object_mean_ratio.csv'))
ecm_object_ratio = pd.pivot(ecm_object_ratio, index='fov', columns='cluster', values='axis_ratio')
for col in ecm_object_ratio.columns:
    if col != 'fov':
        ecm_object_ratio = ecm_object_ratio.rename(columns={col: 'Pixie__major_minor_ratio__' + col})
ecm_object_ratio.reset_index(inplace=True)

# ecm shape normalized difference
ecm_object_diff = pd.read_csv(os.path.join(INTERMEDIATE_DIR, 'ecm_pixel_clustering/shape_analysis/fov_object_mean_diff_norm.csv'))
ecm_object_diff = pd.pivot(ecm_object_diff, index='fov', columns='cluster', values='axis_diff_norm')
for col in ecm_object_diff.columns:
    if col != 'fov':
        ecm_object_diff = ecm_object_diff.rename(columns={col: 'Pixie__major_minor_diff__' + col})
ecm_object_diff.reset_index(inplace=True)

data_frames = [ecm_clusters_density, ecm_clusters_wide, ecm_neighborhood_density, ecm_neighborhoods_prop, ecm_object_ratio, ecm_object_diff]
pixie_ecm = reduce(lambda left, right: pd.merge(left,right,on=['fov'], how='outer'), data_frames)
pixie_ecm.to_csv(os.path.join(FORMATTED_DIR, 'pixie_ecm_stats.csv'), index=False)

ecm_frac_df = ecm_frac_df[ecm_frac_df.fov.isin(np.unique(pixie_ecm.fov))]
ecm_frac_df.to_csv(os.path.join(FORMATTED_DIR, 'ecm_fraction_stats.csv'), index=False)

ecm_df = ecm_df[ecm_df.fov.isin(np.unique(pixie_ecm.fov))]
ecm_df.to_csv(os.path.join(FORMATTED_DIR, 'ecm_cluster_stats.csv'), index=False)
