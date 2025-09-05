import os
import anndata
import pandas as pd
import numpy as np

from postprocessing_utils import TIMEPOINT_COLUMNS, combine_features, generate_feature_rankings, prediction_preprocessing

from SpaceCat.features import SpaceCat


BASE_DIR = '/Volumes/Shared/Noah Greenwald/TONIC_Cohort'
ANALYSIS_DIR = os.path.join(BASE_DIR, 'analysis_files')
INTERMEDIATE_DIR = os.path.join(BASE_DIR, 'intermediate_files')
FORMATTED_DIR = os.path.join(INTERMEDIATE_DIR, 'formatted_files')


############# SKIP ANNDATA GENERATION IF NO CHANGES TO CELL TABLE

# read in data
cell_table = pd.read_csv(os.path.join(ANALYSIS_DIR, 'combined_cell_table_normalized_cell_labels_updated.csv'))
func_table = pd.read_csv(os.path.join(ANALYSIS_DIR, 'cell_table_func_all.csv'))

for col in func_table.columns:
    if col not in ['fov', 'label', 'cell_cluster_broad', 'cell_cluster', 'cell_meta_cluster']:
        if col not in ['H3K9ac_H3K27me3_ratio', 'CD45RO_CD45RB_ratio']:
            func_table[col] = func_table[col].astype(int)
        func_table = func_table.rename(columns={col: f'{col}+'})

cell_table = cell_table.merge(func_table, on=['fov', 'label', 'cell_cluster_broad', 'cell_cluster', 'cell_meta_cluster'])

# remove complex morphology stats for non-Cancer cells
complex_morph_cols = ['convex_hull_resid', 'centroid_dif', 'centroid_dif_nuclear', 'eccentricity', 'num_concavities',
                      'perim_square_over_area', 'convex_hull_resid_nuclear', 'eccentricity_nuclear', 'num_concavities_nuclear']
for col in complex_morph_cols:
    cell_table.loc[cell_table.cell_cluster_broad != 'Cancer', col] = np.nan

# add linear distances
linear_dists = pd.read_csv(os.path.join(INTERMEDIATE_DIR, 'spatial_analysis/cell_neighbor_analysis/cell_distances_filtered.csv'))
cell_table = cell_table.merge(linear_dists, on=['fov', 'label', 'cell_cluster_broad'], how='left')

# add compartment annotations
compartment_annotations = pd.read_csv(os.path.join(INTERMEDIATE_DIR, 'mask_dir/cell_annotation_mask.csv'))
compartment_areas = pd.read_csv(os.path.join(INTERMEDIATE_DIR, 'mask_dir/fov_annotation_mask_area.csv'))
compartment_annotations = compartment_annotations.rename(columns={'mask_name': 'compartment'})
compartment_areas = compartment_areas.rename(columns={'area': 'compartment_area'})

cell_table = cell_table.merge(compartment_annotations, on=['fov', 'label'])
cell_table = cell_table.merge(compartment_areas, on=['fov', 'compartment'])


# define column groupings
markers = ['CD11c', 'CD14', 'CD163', 'CD20', 'CD3', 'CD31', 'CD38', 'CD4', 'CD45', 'CD45RB', 'CD45RO', 'CD56', 'CD57',
           'CD68', 'CD69', 'CD8', 'CK17', 'Calprotectin', 'ChyTr', 'Collagen1', 'ECAD', 'FAP', 'FOXP3', 'Fe', 'Fibronectin',
           'GLUT1', 'H3K27me3', 'H3K9ac', 'HLA1', 'HLADR', 'IDO', 'Ki67', 'LAG3', 'PD1', 'PDL1', 'SMA', 'TBET',
           'TCF1', 'TIM3', 'Vim'] + [col for col in cell_table.columns if '+' in col]
centroid_cols = ['centroid-0', 'centroid-1']
cell_data_cols = ['fov', 'label', 'cell_meta_cluster', 'cell_cluster', 'cell_cluster_broad',
                  'compartment', 'compartment_area',
                  'area', 'area_nuclear', 'nc_ratio', 'cell_size', 'cell_size_nuclear', 'centroid-0', 'centroid-1',
                  'centroid-0_nuclear', 'centroid-1_nuclear', 'major_axis_length',
                  'distance_to__B', 'distance_to__Cancer', 'distance_to__Granulocyte', 'distance_to__Mono_Mac',
                  'distance_to__NK', 'distance_to__Other', 'distance_to__Structural', 'distance_to__T'] + complex_morph_cols

# create anndata from table subsetted for marker info, which will be stored in adata.X
adata = anndata.AnnData(cell_table.loc[:, markers])

# store all other cell data in adata.obs
adata.obs = cell_table.loc[:, cell_data_cols, ]
adata.obs_names = [str(i) for i in adata.obs_names]

# store cell centroid data in adata.obsm
adata.obsm['spatial'] = cell_table.loc[:, centroid_cols].values

# save the anndata object
adata.write_h5ad(os.path.join(ANALYSIS_DIR, 'adata.h5ad'))

#############


# run SpaceCat
adata = anndata.read_h5ad(os.path.join(ANALYSIS_DIR, 'adata.h5ad'))

# read in image level dataframes
fiber_df = pd.read_csv(os.path.join(FORMATTED_DIR, 'fiber_stats_table.csv'))
fiber_tile_df = pd.read_csv(os.path.join(FORMATTED_DIR, 'fiber_stats_per_tile.csv'))
mixing_df = pd.read_csv(os.path.join(FORMATTED_DIR, 'formatted_mixing_scores.csv'))
kmeans_img_proportions = pd.read_csv(os.path.join(FORMATTED_DIR, 'neighborhood_image_proportions.csv'))
kmeans_compartment_proportions = pd.read_csv(os.path.join(FORMATTED_DIR, 'neighborhood_compartment_proportions.csv'))
pixie_ecm = pd.read_csv(os.path.join(FORMATTED_DIR, 'pixie_ecm_stats.csv'))
ecm_frac = pd.read_csv(os.path.join(FORMATTED_DIR, 'ecm_fraction_stats.csv'))
ecm_clusters = pd.read_csv(os.path.join(FORMATTED_DIR, 'ecm_cluster_stats.csv'))

# specify cell type pairs to compute a ratio for
ratio_pairings = [('CD8T', 'CD4T'), ('CD4T', 'Treg'), ('CD8T', 'Treg'), ('CD68_Mac', 'CD163_Mac')]

# specify additional per cell and per image stats
per_cell_stats = [
    ['morphology', 'cell_cluster', ['area', 'area_nuclear', 'nc_ratio', 'convex_hull_resid', 'centroid_dif', 'centroid_dif_nuclear', 'eccentricity', 'num_concavities',
                                    'perim_square_over_area', 'convex_hull_resid_nuclear', 'eccentricity_nuclear', 'num_concavities_nuclear']],
    ['linear_distance', 'cell_cluster_broad', ['distance_to__B', 'distance_to__Cancer', 'distance_to__Granulocyte', 'distance_to__Mono_Mac',
                                               'distance_to__NK', 'distance_to__Other', 'distance_to__Structural', 'distance_to__T']]

]

per_img_stats = [
    ['fiber', fiber_df],
    ['fiber', fiber_tile_df],
    ['mixing_score', mixing_df],
    ['kmeans_cluster', kmeans_img_proportions],
    ['kmeans_cluster', kmeans_compartment_proportions],
    ['pixie_ecm', pixie_ecm],
    ['ecm_fraction', ecm_frac],
    ['ecm_cluster', ecm_clusters]
]

features = SpaceCat(adata, image_key='fov', seg_label_key='label', cell_area_key='area',
                    cluster_key=['cell_cluster', 'cell_cluster_broad'],
                    compartment_key='compartment', compartment_area_key='compartment_area')

# Generate features and save anndata
adata_processed = features.run_spacecat(functional_feature_level='cell_cluster', diversity_feature_level='cell_cluster', pixel_radius=50,
                                        specified_ratios_cluster_key='cell_cluster', specified_ratios=ratio_pairings,
                                        per_cell_stats=per_cell_stats, per_img_stats=per_img_stats)

adata_processed.write_h5ad(os.path.join(ANALYSIS_DIR, 'adata_processed.h5ad'))

# Save finalized tables to csv
adata_processed.uns['combined_feature_data'].to_csv(os.path.join(ANALYSIS_DIR, 'combined_feature_data.csv'), index=False)
adata_processed.uns['combined_feature_data_filtered'].to_csv(os.path.join(ANALYSIS_DIR, 'combined_feature_data_filtered.csv'), index=False)
adata_processed.uns['feature_metadata'].to_csv(os.path.join(ANALYSIS_DIR, 'feature_metadata.csv'), index=False)
adata_processed.uns['excluded_features'].to_csv(os.path.join(ANALYSIS_DIR, 'excluded_features.csv'), index=False)

# clean features and filter out immune_agg compartment
for file in ['combined_feature_data', 'combined_feature_data_filtered', 'feature_metadata', 'excluded_features']:
    df = pd.read_csv(os.path.join(ANALYSIS_DIR, file + '.csv'))
    if file == 'excluded_features':
        df = df[~df.feature_name_unique.str.contains('immune_agg')]
        df.to_csv(os.path.join(ANALYSIS_DIR, file + '.csv'), index=False)
        continue
    df = df[df.compartment != 'immune_agg']
    df = df[~df.feature_name_unique.str.contains('immune_agg')]
    df = df[~np.logical_and(df.feature_type == 'linear_distance', df.compartment != 'all')]
    df = df[~np.logical_and(df.feature_type == 'morphology', df.compartment != 'all')]
    df.loc[df.cell_pop_level.isna(), 'cell_pop_level'] = 'all'

    for comp in ['cancer_core', 'cancer_border', 'stroma_core', 'stroma_border']:
        df.loc[df.feature_name == comp + '__proportion', 'compartment'] = comp
        df.loc[df.feature_name == comp + '__proportion', 'compartment'] = comp
    df.to_csv(os.path.join(ANALYSIS_DIR, file + '.csv'), index=False)

# group by timepoint
harmonized_metadata = pd.read_csv(os.path.join(ANALYSIS_DIR, 'harmonized_metadata.csv'))

fov_data_df = pd.read_csv(os.path.join(ANALYSIS_DIR, 'combined_feature_data.csv'))
fov_data_df = pd.merge(fov_data_df, harmonized_metadata[['Tissue_ID', 'fov']], on='fov', how='left')
grouped = fov_data_df.groupby(['Tissue_ID', 'feature_name', 'feature_name_unique', 'compartment',
                               'cell_pop_level', 'feature_type']).agg({'raw_value': ['mean', 'std'],
                                                                       'normalized_value': ['mean', 'std']})
grouped.columns = ['raw_mean', 'raw_std', 'normalized_mean', 'normalized_std']
grouped = grouped.reset_index()
grouped.to_csv(os.path.join(ANALYSIS_DIR, 'timepoint_features.csv'), index=False)

fov_data_df_filtered = pd.read_csv(os.path.join(ANALYSIS_DIR, 'combined_feature_data_filtered.csv'))
fov_data_df_filtered = pd.merge(fov_data_df_filtered, harmonized_metadata[['Tissue_ID', 'fov']], on='fov', how='left')
grouped = fov_data_df_filtered.groupby(['Tissue_ID', 'feature_name', 'feature_name_unique', 'compartment',
                                 'cell_pop_level', 'feature_type']).agg({'raw_value': ['mean', 'std'],
                                                                            'normalized_value': ['mean', 'std']})

grouped.columns = ['raw_mean', 'raw_std', 'normalized_mean', 'normalized_std']
grouped = grouped.reset_index()
grouped.to_csv(os.path.join(ANALYSIS_DIR, 'timepoint_features_filtered.csv'), index=False)


## 7_create_evolution_df.py converted
study_name = 'TONIC'

timepoint_features = pd.read_csv(os.path.join(ANALYSIS_DIR, 'timepoint_features_filtered.csv'))
timepoint_features_agg = timepoint_features.merge(
    harmonized_metadata[['Tissue_ID', 'Timepoint', 'Patient_ID'] + TIMEPOINT_COLUMNS].drop_duplicates(), on='Tissue_ID',
    how='left')
patient_metadata = pd.read_csv(os.path.join(INTERMEDIATE_DIR, f'metadata/{study_name}_data_per_patient.csv'))

# add evolution features to get finalized features specified by timepoint
combine_features(ANALYSIS_DIR, harmonized_metadata, timepoint_features, timepoint_features_agg, patient_metadata,
                 timepoint_columns=TIMEPOINT_COLUMNS)


## genomics_nivo_outcomes.py converted
patient_metadata = pd.read_csv(os.path.join(INTERMEDIATE_DIR, 'metadata/TONIC_data_per_patient.csv'))
feature_metadata = pd.read_csv(os.path.join(ANALYSIS_DIR, 'feature_metadata.csv'))

#
# To generate the feature rankings, you must have downloaded the patient outcome data.
#
outcome_data = pd.read_csv(os.path.join(INTERMEDIATE_DIR, 'metadata/patient_clinical_data.csv'))

# load previously computed results
combined_df = pd.read_csv(os.path.join(ANALYSIS_DIR, 'timepoint_combined_features.csv'))
combined_df = combined_df.merge(outcome_data, on='Patient_ID')
combined_df = combined_df.loc[combined_df.Clinical_benefit.isin(['Yes', 'No']), :]
combined_df.to_csv(os.path.join(ANALYSIS_DIR, 'timepoint_combined_features_outcome_labels.csv'), index=False)

# generate  pvalues and feature ranking
generate_feature_rankings(ANALYSIS_DIR, combined_df, feature_metadata)


# preprocess feature sets for modeling
df_feature = pd.read_csv(os.path.join(ANALYSIS_DIR, f'timepoint_combined_features_outcome_labels.csv'))
prediction_dir = os.path.join(ANALYSIS_DIR, 'prediction_model')
os.makedirs(prediction_dir, exist_ok=True)
df_feature.to_csv(os.path.join(prediction_dir, 'timepoint_combined_features_outcome_labels.csv'), index=False)

prediction_preprocessing(df_feature, prediction_dir)
os.makedirs(os.path.join(prediction_dir, 'patient_outcomes'), exist_ok=True)
