# script to generate summary stats for each fov
import os

import itertools
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import spearmanr
import skimage.io as io

from alpineer.io_utils import list_folders

# This script takes the many distinct dataframes generated from the previous script and
# combines them together into a single unified dataframe for downstream analysis


base_dir = '/Volumes/Shared/Noah Greenwald/TONIC_Cohort/'
intermediate_dir = os.path.join(base_dir, 'intermediate_files')
output_dir = os.path.join(base_dir, 'output_files')
analysis_dir =  os.path.join(base_dir, 'analysis_files')

TIMEPOINT_NAMES = ['primary', 'baseline', 'pre_nivo', 'on_nivo']
study_name = 'TONIC'

# load datasets
cluster_df_core = pd.read_csv(os.path.join(output_dir, 'cluster_df_per_core.csv'))
metadata_df_core = pd.read_csv(os.path.join(intermediate_dir, f'metadata/{study_name}_data_per_core.csv'))
functional_df_core = pd.read_csv(os.path.join(output_dir, 'functional_df_per_core_filtered_deduped.csv'))
morph_df_core = pd.read_csv(os.path.join(output_dir, 'morph_df_per_core_filtered_deduped.csv'))
mixing_df = pd.read_csv(os.path.join(output_dir, 'formatted_mixing_scores.csv'))
diversity_df = pd.read_csv(os.path.join(output_dir, 'diversity_df_per_core_filtered_deduped.csv'))
distance_df = pd.read_csv(os.path.join(output_dir, 'distance_df_per_core_deduped.csv'))
fiber_df = pd.read_csv(os.path.join(output_dir, 'fiber_df_per_core.csv'))
fiber_tile_df = pd.read_csv(os.path.join(output_dir, 'fiber_df_per_tile.csv'))
ecm_df = pd.read_csv(os.path.join(intermediate_dir, 'ecm/fov_cluster_counts.csv'))
ecm_clusters = pd.read_csv(os.path.join(intermediate_dir, 'ecm_pixel_clustering/fov_pixel_cluster_counts.csv'))
ecm_object_ratio = pd.read_csv(os.path.join(intermediate_dir, 'ecm_pixel_clustering/shape_analysis/fov_object_mean_ratio.csv'))
ecm_object_diff = pd.read_csv(os.path.join(intermediate_dir, 'ecm_pixel_clustering/shape_analysis/fov_object_mean_diff_norm.csv'))
ecm_neighborhoods = pd.read_csv(os.path.join(intermediate_dir, 'ecm_pixel_clustering/neighborhood/fov_neighborhood_counts.csv'))
image_clusters = pd.read_csv(os.path.join(output_dir, 'neighborhood_image_proportions.csv'))
compartment_clusters = pd.read_csv(os.path.join(output_dir, 'neighborhood_compartment_proportions.csv'))

# load metadata
harmonized_metadata_df = pd.read_csv(os.path.join(analysis_dir, 'harmonized_metadata.csv'))
compartment_area = pd.read_csv(os.path.join(intermediate_dir, 'mask_dir/fov_annotation_mask_area.csv'))


# compute shannon diversity from list of proportions
def shannon_diversity(proportions):
    proportions = [prop for prop in proportions if prop > 0]
    return -np.sum(proportions * np.log2(proportions))


# create lookup table for mapping individual cell types to broader categories
broad_to_narrow = {'Cancer': ['Cancer_CD56', 'Cancer_CK17', 'Cancer_Ecad', 'Cancer_SMA',
                              'Cancer_Vim', 'Cancer_Other', 'Cancer_Mono', 'Cancer_1', 'Cancer_2', 'Cancer_3', 'Cancer'],
                   'Immune': ['CD68', 'CD163', 'CD68_CD163_DP', 'CD4_Mono', 'CD14', 'CD11c_HLADR',
                              'CD20', 'CD56', 'Neutrophil', 'Mast', 'CD4T', 'CD4T_HLADR', 'CD8T', 'Treg', 'CD3_DN',
                              'CD4T_CD8T_DP', 'Immune_Other', 'CD68_Mac', 'CD163_Mac', 'Mac_Other',
                              'Monocyte', 'APC', 'B', 'CD4T', 'CD8T', 'Treg', 'T_Other',
                              'Neutrophil', 'Mast', 'NK', 'Mono_Mac', 'Granulocyte', 'T'],
                   'Structural': ['Stroma_Collagen', 'Stroma_Fibronectin', 'VIM', 'CD31', 'CD31_VIM',
                              'FAP', 'FAP_SMA', 'SMA', 'Other', 'Endothelium', 'Fibroblast',
                              'Structural', 'CAF', 'Smooth_Muscle'],
                   }
# reverse lookup table
narrow_to_broad = {}
for broad, narrow in broad_to_narrow.items():
    for n in narrow:
        narrow_to_broad[n] = broad

fov_data = []

# compute diversity of different levels of granularity
diversity_features = [['cluster_broad_freq', 'cluster_broad_diversity'],
                      ['immune_freq', 'immune_diversity'],
                      ['cancer_freq', 'cancer_diversity'],
                      ['structural_freq', 'structural_diversity']]

# diversity_features = [['meta_cluster_freq', 'meta_cluster_diversity']]

for cluster_name, feature_name in diversity_features:
    input_df = cluster_df_core[cluster_df_core['metric'].isin([cluster_name])]
    for compartment in ['cancer_core', 'cancer_border', 'stroma_core', 'stroma_border', 'all']:
        compartment_df = input_df[input_df.subset == compartment].copy()
        wide_df = pd.pivot(compartment_df, index='fov', columns=['cell_type'], values='value')
        wide_df['value'] = wide_df.apply(shannon_diversity, axis=1)
        wide_df.reset_index(inplace=True)
        wide_df['feature_name'] = feature_name
        wide_df['compartment'] = compartment

        if compartment == 'all':
            wide_df['feature_name_unique'] = feature_name
        else:
            wide_df['feature_name_unique'] = feature_name + '_' + compartment

        if cluster_name == 'cluster_broad_freq':
            cell_pop = 'all'
            cell_pop_level = 'broad'
        else:
            cell_pop = cluster_name.split('_')[0]
            cell_pop = cell_pop[0].upper() + cell_pop[1:]
            cell_pop_level = 'broad'

        wide_df['cell_pop'] = cell_pop
        wide_df['cell_pop_level'] = cell_pop_level
        wide_df['feature_type'] = 'region_diversity'
        wide_df['feature_type_detail'] = 'region_diversity'
        wide_df['feature_type_detail_2'] = ''
        wide_df = wide_df[['fov', 'value', 'feature_name', 'feature_name_unique', 'compartment', 'cell_pop', 'cell_pop_level',
                           'feature_type', 'feature_type_detail', 'feature_type_detail_2']]
        fov_data.append(wide_df)


# compute abundance of cell types
abundance_features = [['cluster_density', 'cluster_density', 'med'],
                      ['total_cell_density', 'total_density', 'broad'],
                      ['cluster_broad_density', 'cluster_broad_density', 'broad']]

# abundance_features = [['meta_cluster_density', 'meta_cluster_density', 'narrow']]

for cluster_name, feature_name, cell_pop_level in abundance_features:
    input_df = cluster_df_core[cluster_df_core['metric'].isin([cluster_name])]
    if cluster_name == 'cluster_density':
        # B and NK are the same as cluster_broad, keep just cluster broad
        input_df = input_df[~input_df.cell_type.isin(['B', 'NK'])]

    for compartment in ['cancer_core', 'cancer_border', 'stroma_core', 'stroma_border', 'all']:
        compartment_df = input_df[input_df.subset == compartment].copy()
        compartment_df['feature_name'] = compartment_df.cell_type + '__' + feature_name
        compartment_df = compartment_df.rename(columns={'subset': 'compartment'})
        if cluster_name != 'total_cell_density':
            compartment_df['cell_pop'] = compartment_df.cell_type.apply(lambda x: narrow_to_broad[x])
        else:
            compartment_df['cell_pop'] = 'all'

        if compartment == 'all':
            compartment_df['feature_name_unique'] = compartment_df.cell_type + '__' + feature_name
        else:
            compartment_df['feature_name_unique'] = compartment_df.cell_type + '__' + feature_name + '__' + compartment

        compartment_df['cell_pop_level'] = cell_pop_level
        compartment_df['feature_type'] = cluster_name.split('_')[-1]
        compartment_df['feature_type_detail'] = compartment_df.cell_type
        compartment_df['feature_type_detail_2'] = ''
        compartment_df = compartment_df[['fov', 'value', 'feature_name', 'feature_name_unique', 'compartment', 'cell_pop',
                                         'cell_pop_level', 'feature_type', 'feature_type_detail', 'feature_type_detail_2']]
        fov_data.append(compartment_df)

# compute ratio of broad cell type abundances
input_df = cluster_df_core[cluster_df_core['metric'].isin(['cluster_broad_density'])]
for compartment in ['cancer_core', 'cancer_border', 'stroma_core', 'stroma_border', 'all']:
    compartment_df = input_df[input_df.subset == compartment].copy()
    cell_types = compartment_df.cell_type.unique()

    # put cancer and structural cells last
    cell_types = [cell_type for cell_type in cell_types if cell_type not in ['Cancer', 'Structural']]
    cell_types = cell_types + ['Cancer', 'Structural']

    for cell_type1, cell_type2 in itertools.combinations(cell_types, 2):
        cell_type1_df = compartment_df[compartment_df.cell_type == cell_type1].copy()
        cell_type2_df = compartment_df[compartment_df.cell_type == cell_type2].copy()

        # only keep FOVS with at least one cell type over the minimum density
        minimum_density = 0.0005
        cell_type1_mask = cell_type1_df.value > minimum_density
        cell_type2_mask = cell_type2_df.value > minimum_density
        cell_mask = cell_type1_mask.values | cell_type2_mask.values

        cell_type1_df = cell_type1_df[cell_mask]
        cell_type2_df = cell_type2_df[cell_mask]

        # add minimum density to avoid log2(0)
        cell_type1_df['ratio'] = np.log2((cell_type1_df.value.values + minimum_density) /
                                         (cell_type2_df.value.values + minimum_density))

        cell_type1_df['value'] = cell_type1_df.ratio.values

        cell_type1_df['feature_name'] = cell_type1 + '__' + cell_type2 + '__ratio'
        cell_type1_df['compartment'] = compartment

        if compartment == 'all':
            cell_type1_df['feature_name_unique'] = cell_type1 + '__' + cell_type2 + '__ratio'
        else:
            cell_type1_df['feature_name_unique'] = cell_type1 + '__' + cell_type2 + '__ratio__' + compartment

        cell_type1_df['cell_pop'] = 'multiple'
        cell_type1_df['cell_pop_level'] = 'broad'
        cell_type1_df['feature_type'] = 'density_ratio'
        cell_type1_df['feature_type_detail'] = cell_type1
        cell_type1_df['feature_type_detail_2'] = cell_type2
        cell_type1_df = cell_type1_df[['fov', 'value', 'feature_name', 'feature_name_unique', 'compartment', 'cell_pop',
                                       'cell_pop_level', 'feature_type', 'feature_type_detail', 'feature_type_detail_2']]
        fov_data.append(cell_type1_df)

# compute ratio of specific cell type abundances
cell_ratios = [('CD8T', 'CD4T'), ('CD4T', 'Treg'), ('CD8T', 'Treg'), ('CD68_Mac', 'CD163_Mac')]
input_df = cluster_df_core[cluster_df_core.metric == 'cluster_density'].copy()
for compartment in ['cancer_core', 'cancer_border', 'stroma_core', 'stroma_border', 'all']:
    compartment_df = input_df[input_df.subset == compartment].copy()
    for cell_type1, cell_type2 in cell_ratios:
        cell_type1_df = compartment_df[compartment_df.cell_type == cell_type1].copy()
        cell_type2_df = compartment_df[compartment_df.cell_type == cell_type2].copy()

        # only keep FOVS with at least one cell type over the minimum density
        minimum_density = 0.0005
        cell_type1_mask = cell_type1_df.value > minimum_density
        cell_type2_mask = cell_type2_df.value > minimum_density
        cell_mask = cell_type1_mask.values | cell_type2_mask.values

        cell_type1_df = cell_type1_df[cell_mask]
        cell_type2_df = cell_type2_df[cell_mask]

        # add minimum density to avoid log2(0)
        cell_type1_df['ratio'] = np.log2((cell_type1_df.value.values + minimum_density) /
                                         (cell_type2_df.value.values + minimum_density))

        cell_type1_df['value'] = cell_type1_df.ratio.values
        cell_type1_df['feature_name'] = cell_type1 + '__' + cell_type2 + '__ratio'

        if compartment == 'all':
            cell_type1_df['feature_name_unique'] = cell_type1 + '__' + cell_type2 + '__ratio'
        else:
            cell_type1_df['feature_name_unique'] = cell_type1 + '__' + cell_type2 + '__ratio__' + compartment

        cell_type1_df['compartment'] = compartment
        cell_type1_df['cell_pop'] = 'Immune'
        cell_type1_df['cell_pop_level'] = 'med'
        cell_type1_df['feature_type'] = 'density_ratio'
        cell_type1_df['feature_type_detail'] = cell_type1
        cell_type1_df['feature_type_detail_2'] = cell_type2
        cell_type1_df = cell_type1_df[['fov', 'value', 'feature_name', 'feature_name_unique','compartment', 'cell_pop',
                                       'cell_pop_level',  'feature_type', 'feature_type_detail', 'feature_type_detail_2']]
        fov_data.append(cell_type1_df)


# compute proportion of cells in a given cell type
cell_groups = {'Cancer': ['Cancer_1', 'Cancer_2', 'Cancer_3'],
               'Mono_Mac': ['CD68_Mac', 'CD163_Mac', 'Mac_Other', 'Monocyte', 'APC'],
               'T': ['CD4T', 'CD8T', 'Treg', 'T_Other'],
               'Granulocyte': ['Neutrophil', 'Mast'],
               'Structural': ['Endothelium', 'CAF', 'Fibroblast', 'Smooth_Muscle']}

input_df = cluster_df_core[cluster_df_core.metric == 'cluster_density'].copy()
for compartment in ['cancer_core', 'cancer_border', 'stroma_core', 'stroma_border', 'all']:
    compartment_df = input_df[input_df.subset == compartment].copy()
    for broad_cell_type, cell_types in cell_groups.items():
        # get the total for all cell types
        cell_type_df = compartment_df[compartment_df.cell_type.isin(cell_types)].copy()
        grouped_df = cell_type_df[['fov', 'value']].groupby('fov').sum().reset_index()
        grouped_df.columns = ['fov', 'fov_sum']

        # normalize each cell type by the total
        cell_type_df = cell_type_df.merge(grouped_df, on='fov')
        idx_nonzero = np.where(cell_type_df.fov_sum != 0)[0]
        cell_type_df = cell_type_df.iloc[idx_nonzero, :].copy()
        cell_type_df['value'] = cell_type_df.value / cell_type_df.fov_sum

        cell_type_df['feature_name'] = cell_type_df.cell_type + '__proportion_of__' + broad_cell_type

        if compartment == 'all':
            cell_type_df['feature_name_unique'] = cell_type_df.cell_type + '__proportion_of__' + broad_cell_type
        else:
            cell_type_df['feature_name_unique'] = cell_type_df.cell_type + '__proportion_of__' + broad_cell_type + '__' + compartment

        cell_type_df['compartment'] = compartment
        cell_type_df['cell_pop'] = broad_cell_type
        cell_type_df['cell_pop_level'] = 'med'
        cell_type_df['feature_type'] = 'density_proportion'
        cell_type_df['feature_type_detail'] = cell_type_df.cell_type
        cell_type_df['feature_type_detail_2'] = broad_cell_type
        cell_type_df = cell_type_df[['fov', 'value', 'feature_name', 'feature_name_unique', 'compartment', 'cell_pop',
                                     'cell_pop_level', 'feature_type', 'feature_type_detail', 'feature_type_detail_2']]
        fov_data.append(cell_type_df)


# compute functional marker positivity for different levels of granularity
functional_features = [['cluster_freq', 'med'],
                       ['total_freq', 'broad']]

# functional_features = [['meta_cluster_freq', 'narrow']]
for functional_name, cell_pop_level in functional_features:
    input_df = functional_df_core[functional_df_core['metric'].isin([functional_name])]
    for compartment in ['cancer_core', 'cancer_border', 'stroma_core', 'stroma_border', 'all']:
        compartment_df = input_df[input_df.subset == compartment].copy()
        compartment_df['feature_name'] = compartment_df.functional_marker + '+__' + compartment_df.cell_type

        if compartment == 'all':
            compartment_df['feature_name_unique'] = compartment_df.functional_marker + '+__' + compartment_df.cell_type
        else:
            compartment_df['feature_name_unique'] = compartment_df.functional_marker + '+__' + compartment_df.cell_type + '__' + compartment

        compartment_df = compartment_df.rename(columns={'subset': 'compartment'})

        if functional_name != 'total_freq':
            compartment_df['cell_pop'] = compartment_df.cell_type.apply(lambda x: narrow_to_broad[x])
        else:
            compartment_df['cell_pop'] = 'all'

        compartment_df['cell_pop_level'] = cell_pop_level
        compartment_df['feature_type'] = 'functional_marker'
        compartment_df['feature_type_detail'] = compartment_df.functional_marker
        compartment_df['feature_type_detail_2'] = ''
        compartment_df = compartment_df[['fov', 'value', 'feature_name','feature_name_unique','compartment', 'cell_pop',
                                         'cell_pop_level', 'feature_type', 'feature_type_detail', 'feature_type_detail_2']]
        fov_data.append(compartment_df)

# compute morphology features for different levels of granularity
morphology_features = [['cluster_freq', 'med'],
                          ['total_freq', 'broad']]

# morphology_features = [['meta_cluster_freq', 'narrow']]
for morphology_name, cell_pop_level in morphology_features:
    input_df = morph_df_core[morph_df_core['metric'].isin([morphology_name])]
    for compartment in ['all']:
        compartment_df = input_df[input_df.subset == compartment].copy()
        compartment_df['feature_name'] = compartment_df.morphology_feature + '__' + compartment_df.cell_type

        if compartment == 'all':
            compartment_df['feature_name_unique'] = compartment_df.morphology_feature + '__' + compartment_df.cell_type
        else:
            compartment_df[ 'feature_name_unique'] = compartment_df.morphology_feature + '__' + compartment_df.cell_type + '__' + compartment

        compartment_df = compartment_df.rename(columns={'subset': 'compartment'})

        if morphology_name != 'total_freq':
            compartment_df['cell_pop'] = compartment_df.cell_type.apply(lambda x: narrow_to_broad[x])
        else:
            compartment_df['cell_pop'] = 'all'

        compartment_df['cell_pop_level'] = cell_pop_level
        compartment_df['feature_type'] = 'morphology'
        compartment_df['feature_type_detail'] = compartment_df.morphology_feature
        compartment_df['feature_type_detail_2'] = ''
        compartment_df = compartment_df[['fov', 'value', 'feature_name','feature_name_unique','compartment', 'cell_pop',
                                         'cell_pop_level', 'feature_type', 'feature_type_detail', 'feature_type_detail_2']]
        fov_data.append(compartment_df)


# compute diversity features
input_df = diversity_df[diversity_df['metric'] == 'cluster_freq']
input_df = input_df[input_df.diversity_feature == 'diversity_cell_cluster']
for compartment in ['cancer_core', 'cancer_border', 'stroma_core', 'stroma_border', 'all']:
    compartment_df = input_df[input_df.subset == compartment].copy()
    compartment_df['feature_name'] = compartment_df.diversity_feature + '__' + compartment_df.cell_type

    if compartment == 'all':
        compartment_df['feature_name_unique'] = compartment_df.diversity_feature + '__' + compartment_df.cell_type
    else:
        compartment_df['feature_name_unique'] = compartment_df.diversity_feature + '__' + compartment_df.cell_type + '__' + compartment

    compartment_df = compartment_df.rename(columns={'subset': 'compartment'})

    compartment_df['cell_pop'] = compartment_df.cell_type.apply(lambda x: narrow_to_broad[x])
    compartment_df['cell_pop_level'] = 'med'
    compartment_df['feature_type'] = 'cell_diversity'
    compartment_df['feature_type_detail'] = compartment_df.cell_type
    compartment_df['feature_type_detail_2'] = compartment_df.diversity_feature
    compartment_df = compartment_df[['fov', 'value', 'feature_name','feature_name_unique','compartment', 'cell_pop',
                                     'cell_pop_level', 'feature_type', 'feature_type_detail', 'feature_type_detail_2']]
    fov_data.append(compartment_df)

# compute distance features
input_df = distance_df[distance_df['metric'].isin(['cluster_broad_freq'])]
for compartment in ['cancer_core', 'cancer_border', 'stroma_core', 'stroma_border', 'all']:
    compartment_df = input_df[input_df.subset == compartment].copy()
    compartment_df['feature_name'] = compartment_df.cell_type + '__distance_to__' + compartment_df.linear_distance

    if compartment == 'all':
        compartment_df['feature_name_unique'] = compartment_df.cell_type + '__distance_to__' + compartment_df.linear_distance
    else:
        compartment_df[ 'feature_name_unique'] = compartment_df.cell_type + '__distance_to__' + compartment_df.linear_distance + '__' + compartment

    compartment_df = compartment_df.rename(columns={'subset': 'compartment'})

    compartment_df['cell_pop'] = compartment_df.cell_type.apply(lambda x: narrow_to_broad[x])
    compartment_df['cell_pop_level'] = 'broad'
    compartment_df['feature_type'] = 'linear_distance'
    compartment_df['feature_type_detail'] = compartment_df.cell_type
    compartment_df['feature_type_detail_2'] = compartment_df.linear_distance
    compartment_df = compartment_df[['fov', 'value', 'feature_name','feature_name_unique','compartment', 'cell_pop',
                                        'cell_pop_level', 'feature_type', 'feature_type_detail', 'feature_type_detail_2']]

    fov_data.append(compartment_df)

# compute compartment abundance and ratios
compartments = ['cancer_core', 'cancer_border', 'stroma_core', 'stroma_border']
for idx, compartment in enumerate(compartments):
    compartment_df = compartment_area[compartment_area.compartment == compartment].copy()
    total_area = compartment_area[compartment_area.compartment == 'all']
    compartment_df['value'] = compartment_df.area.values / total_area.area.values
    compartment_df['feature_name'] = compartment + '__proportion'
    compartment_df['feature_name_unique'] = compartment + '__proportion'
    compartment_df['compartment'] = compartment
    compartment_df['cell_pop'] = 'all'
    compartment_df['cell_pop_level'] = 'broad'
    compartment_df['feature_type'] = 'compartment_area'
    compartment_df['feature_type_detail'] = compartment
    compartment_df['feature_type_detail_2'] = ''
    compartment_df = compartment_df[['fov', 'value', 'feature_name', 'feature_name_unique','compartment', 'cell_pop',
                                     'cell_pop_level', 'feature_type', 'feature_type_detail', 'feature_type_detail_2']]
    fov_data.append(compartment_df)

    # now look at combinations of compartments, except for rare ones
    if idx > 2:
        continue
    compartment2 = compartments[idx + 1]
    compartment2_df = compartment_area[compartment_area.compartment == compartment2].copy()
    compartment2_df['value'] = compartment2_df.area.values / total_area.area.values

    minimum_abundance = 0.01
    compartment1_keep_mask = compartment_df.value > minimum_abundance
    compartment2_keep_mask = compartment2_df.value > minimum_abundance
    keep_mask = compartment1_keep_mask.values | compartment2_keep_mask.values
    compartment_df = compartment_df[keep_mask]
    compartment2_df = compartment2_df[keep_mask]
    compartment2_df['value'] = np.log2((compartment_df.value.values + minimum_abundance) /
                                       (compartment2_df.value.values + minimum_abundance))

    # add metadata
    compartment2_df['feature_name'] = compartment + '__' + compartment2 + '__log2_ratio'
    compartment2_df['feature_name_unique'] = compartment + '__' + compartment2 + '__log2_ratio'
    compartment2_df['compartment'] = 'all'
    compartment2_df['cell_pop'] = 'all'
    compartment2_df['cell_pop_level'] = 'broad'
    compartment2_df['feature_type'] = 'compartment_area_ratio'
    compartment2_df['feature_type_detail'] = compartment
    compartment2_df['feature_type_detail_2'] = compartment2
    compartment2_df = compartment2_df[['fov', 'value', 'feature_name', 'feature_name_unique', 'compartment', 'cell_pop',
                                       'cell_pop_level', 'feature_type', 'feature_type_detail', 'feature_type_detail_2']]
    fov_data.append(compartment2_df)


# integrate mixing scores
mixing_df = mixing_df.dropna()
mixing_df = mixing_df.rename(columns={'mixing_score': 'feature_name'})
mixing_df['feature_name_unique'] = mixing_df.feature_name
mixing_df['compartment'] = 'all'
mixing_df['cell_pop'] = 'multiple'
mixing_df['cell_pop_level'] = 'broad'
mixing_df['feature_type'] = 'mixing_score'
mixing_df['feature_type_detail'] = 'mixing_score'
mixing_df['feature_type_detail_2'] = ''
mixing_df = mixing_df[['fov', 'value', 'feature_name', 'feature_name_unique', 'compartment',
                       'cell_pop', 'cell_pop_level', 'feature_type', 'feature_type_detail', 'feature_type_detail_2']]
fov_data.append(mixing_df)

# add kmeans cluster proportions
# image proportions
image_clusters = image_clusters.rename(columns={'proportion': 'value'})
image_clusters['feature_name'] = 'cluster_' + image_clusters.kmeans_neighborhood.astype(str) + '__proportion'
image_clusters['feature_name_unique'] = 'cluster_' + image_clusters.kmeans_neighborhood.astype(str) + '__proportion'
image_clusters['compartment'] = 'all'
image_clusters['cell_pop'] = 'kmeans_cluster'
image_clusters['cell_pop_level'] = 'med'
image_clusters['feature_type'] = 'kmeans_cluster'
image_clusters['feature_type_detail'] = 'cluster_' + image_clusters.kmeans_neighborhood.astype(str)
image_clusters['feature_type_detail_2'] = ''
image_clusters = image_clusters[['fov', 'value', 'feature_name', 'feature_name_unique',
                                 'compartment', 'cell_pop', 'cell_pop_level', 'feature_type', 'feature_type_detail', 'feature_type_detail_2']]
fov_data.append(image_clusters)

# compartment proportions
compartment_clusters = compartment_clusters.rename(columns={'proportion': 'value'})
compartment_clusters['feature_name'] = 'cluster_' + compartment_clusters.kmeans_neighborhood.astype(str) + '__proportion'
compartment_clusters['feature_name_unique'] = 'cluster_' + compartment_clusters.kmeans_neighborhood.astype(str) + '__proportion' + '__' + compartment_clusters.mask_name
compartment_clusters['compartment'] = compartment_clusters.mask_name
compartment_clusters['cell_pop'] = 'kmeans_cluster'
compartment_clusters['cell_pop_level'] = 'med'
compartment_clusters['feature_type'] = 'kmeans_cluster'
compartment_clusters['feature_type_detail'] = 'cluster_' + compartment_clusters.kmeans_neighborhood.astype(str)
compartment_clusters['feature_type_detail_2'] = compartment_clusters.mask_name
compartment_clusters = compartment_clusters[['fov', 'value', 'feature_name', 'feature_name_unique',
                                             'compartment', 'cell_pop', 'cell_pop_level', 'feature_type', 'feature_type_detail', 'feature_type_detail_2']]
fov_data.append(compartment_clusters)

# add ecm clustering results
ecm_frac = 1 - ecm_df.No_ECM.values
ecm_df['Cold_Coll_Norm'] = ecm_df.Cold_Coll.values / ecm_frac
ecm_df['Hot_Coll_Norm'] = ecm_df.Hot_Coll.values / ecm_frac
ecm_df.loc[ecm_frac == 0, ['Cold_Coll_Norm', 'Hot_Coll_Norm']] = 0

for col_name in ['Cold_Coll', 'Hot_Coll', 'No_ECM', 'Cold_Coll_Norm', 'Hot_Coll_Norm']:
    ecm_df_subset = ecm_df[['fov', col_name]]
    ecm_df_subset = ecm_df_subset.rename(columns={col_name: 'value'})
    ecm_df_subset['feature_name'] = col_name + '__proportion'
    ecm_df_subset['feature_name_unique'] = col_name + '__proportion'
    ecm_df_subset['compartment'] = 'all'
    ecm_df_subset['cell_pop'] = 'ecm'
    ecm_df_subset['cell_pop_level'] = 'broad'
    ecm_df_subset['feature_type'] = 'ecm_cluster'
    ecm_df_subset['feature_type_detail'] = col_name
    ecm_df_subset['feature_type_detail_2'] = ''
    ecm_df_subset = ecm_df_subset[['fov', 'value', 'feature_name', 'feature_name_unique',
                                   'compartment', 'cell_pop', 'cell_pop_level', 'feature_type', 'feature_type_detail', 'feature_type_detail_2']]
    fov_data.append(ecm_df_subset)

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
# ecm_frac_df.to_csv(os.path.join(intermediate_dir, 'ecm/ecm_fraction_fov.csv'), index=False)

# add ecm fraction
ecm_frac_df = pd.read_csv(os.path.join(intermediate_dir, 'ecm/ecm_fraction_fov.csv'))
ecm_frac_df = ecm_frac_df.rename(columns={'ecm_frac': 'value'})
ecm_frac_df['feature_name'] = 'ecm_fraction'
ecm_frac_df['feature_name_unique'] = 'ecm_fraction'
ecm_frac_df['compartment'] = 'all'
ecm_frac_df['cell_pop'] = 'ecm'
ecm_frac_df['cell_pop_level'] = 'broad'
ecm_frac_df['feature_type'] = 'ecm_fraction'
ecm_frac_df['feature_type_detail'] = 'ecm_fraction'
ecm_frac_df['feature_type_detail_2'] = ''
ecm_frac_df = ecm_frac_df[['fov', 'value', 'feature_name', 'feature_name_unique',
                            'compartment', 'cell_pop', 'cell_pop_level', 'feature_type', 'feature_type_detail', 'feature_type_detail_2']]

fov_data.append(ecm_frac_df)


# add ecm pixel cluster density
area_df = pd.read_csv(os.path.join(intermediate_dir, 'mask_dir/fov_annotation_mask_area.csv'))
area_df = area_df.loc[area_df.compartment == 'all', ['fov', 'area']]
ecm_clusters_density = pd.merge(ecm_clusters, area_df, on='fov')
ecm_clusters_density['density'] = ecm_clusters_density['counts'] / ecm_clusters_density['area']

ecm_clusters_density['feature_name'] = 'Pixie__cluster__' + ecm_clusters_density['pixel_meta_cluster_rename'] + '__density'
ecm_clusters_density['feature_name_unique'] = 'Pixie__cluster__' + ecm_clusters_density['pixel_meta_cluster_rename'] + '__density'
ecm_clusters_density['compartment'] = 'all'
ecm_clusters_density['cell_pop'] = 'ecm'
ecm_clusters_density['cell_pop_level'] = 'broad'
ecm_clusters_density['feature_type'] = 'pixie_ecm'
ecm_clusters_density['feature_type_detail'] = ecm_clusters_density['pixel_meta_cluster_rename']
ecm_clusters_density['feature_type_detail_2'] = ''
ecm_clusters_density = ecm_clusters_density.rename(columns={'density': 'value'})
ecm_clusters_density = ecm_clusters_density[['fov', 'value', 'feature_name', 'feature_name_unique',
                                'compartment', 'cell_pop', 'cell_pop_level', 'feature_type', 'feature_type_detail', 'feature_type_detail_2']]

fov_data.append(ecm_clusters_density)

# add ecm pixel cluster proportion
ecm_clusters_wide = pd.pivot(ecm_clusters, index='fov', columns='pixel_meta_cluster_rename', values='counts')
ecm_clusters_wide = ecm_clusters_wide.apply(lambda x: x / x.sum(), axis=1)
ecm_clusters_wide = ecm_clusters_wide.reset_index()

ecm_clusters = pd.melt(ecm_clusters_wide, id_vars='fov', var_name='ecm_cluster_name', value_name='value')
ecm_clusters['feature_name'] = 'Pixie__cluster__' + ecm_clusters['ecm_cluster_name'] + '__proportion'
ecm_clusters['feature_name_unique'] = 'Pixie__cluster__' + ecm_clusters['ecm_cluster_name'] + '__proportion'
ecm_clusters['compartment'] = 'all'
ecm_clusters['cell_pop'] = 'ecm'
ecm_clusters['cell_pop_level'] = 'broad'
ecm_clusters['feature_type'] = 'pixie_ecm'
ecm_clusters['feature_type_detail'] = ecm_clusters['ecm_cluster_name']
ecm_clusters['feature_type_detail_2'] = ''
ecm_clusters = ecm_clusters[['fov', 'value', 'feature_name', 'feature_name_unique',
                                'compartment', 'cell_pop', 'cell_pop_level', 'feature_type', 'feature_type_detail', 'feature_type_detail_2']]

fov_data.append(ecm_clusters)

# add ecm neighborhood density
ecm_neighborhoods = ecm_neighborhoods.rename(columns={'Cluster1': 'Collagen', 'Cluster2': 'SMA', 'Cluster3': 'Collagen_Vim', 'Cluster4': 'Vim_FAP', 'Cluster5': 'Fibronectin'})
ecm_neighborhood_density = pd.melt(ecm_neighborhoods.iloc[:, :-1], id_vars='fov', var_name='ecm_neighborhood', value_name='counts')
ecm_neighborhood_density = pd.merge(ecm_neighborhood_density, area_df, on='fov')
ecm_neighborhood_density['density'] = ecm_neighborhood_density['counts'] / ecm_neighborhood_density['area']

ecm_neighborhood_density['feature_name'] = 'Pixie__neighborhood__' + ecm_neighborhood_density['ecm_neighborhood'] + '__density'
ecm_neighborhood_density['feature_name_unique'] = 'Pixie__neighborhood__' + ecm_neighborhood_density['ecm_neighborhood'] + '__density'
ecm_neighborhood_density['compartment'] = 'all'
ecm_neighborhood_density['cell_pop'] = 'ecm'
ecm_neighborhood_density['cell_pop_level'] = 'broad'
ecm_neighborhood_density['feature_type'] = 'pixie_ecm'
ecm_neighborhood_density['feature_type_detail'] = ecm_neighborhood_density['ecm_neighborhood']
ecm_neighborhood_density['feature_type_detail_2'] = ''
ecm_neighborhood_density = ecm_neighborhood_density.rename(columns={'density': 'value'})
ecm_neighborhood_density = ecm_neighborhood_density[['fov', 'value', 'feature_name', 'feature_name_unique',
                                'compartment', 'cell_pop', 'cell_pop_level', 'feature_type', 'feature_type_detail', 'feature_type_detail_2']]

fov_data.append(ecm_neighborhood_density)

# add ecm neighborhood proportion
ecm_neighborhoods_prop = ecm_neighborhoods.copy()
for col in ecm_neighborhoods_prop.columns[1:-1]:
    ecm_neighborhoods_prop[col] = ecm_neighborhoods_prop[col] / ecm_neighborhoods_prop['total']

ecm_neighborhoods_prop = pd.melt(ecm_neighborhoods_prop.iloc[:, :-1], id_vars='fov', var_name='ecm_neighborhood', value_name='value')
ecm_neighborhoods_prop['feature_name'] = 'Pixie__neighborhood__' + ecm_neighborhoods_prop['ecm_neighborhood'] + '__proportion'
ecm_neighborhoods_prop['feature_name_unique'] = 'Pixie__neighborhood__' + ecm_neighborhoods_prop['ecm_neighborhood'] + '__proportion'
ecm_neighborhoods_prop['compartment'] = 'all'
ecm_neighborhoods_prop['cell_pop'] = 'ecm'
ecm_neighborhoods_prop['cell_pop_level'] = 'broad'
ecm_neighborhoods_prop['feature_type'] = 'pixie_ecm'
ecm_neighborhoods_prop['feature_type_detail'] = ecm_neighborhoods_prop['ecm_neighborhood']
ecm_neighborhoods_prop['feature_type_detail_2'] = ''
ecm_neighborhoods_prop = ecm_neighborhoods_prop[['fov', 'value', 'feature_name', 'feature_name_unique',
                                'compartment', 'cell_pop', 'cell_pop_level', 'feature_type', 'feature_type_detail', 'feature_type_detail_2']]
fov_data.append(ecm_neighborhoods_prop)

# add ecm shape axis ratio
ecm_object_ratio['feature_name'] = 'Pixie__major_minor_ratio__' + ecm_object_ratio['cluster']
ecm_object_ratio['feature_name_unique'] = 'Pixie__major_minor_ratio__' + ecm_object_ratio['cluster']
ecm_object_ratio['compartment'] = 'all'
ecm_object_ratio['cell_pop'] = 'ecm'
ecm_object_ratio['cell_pop_level'] = 'broad'
ecm_object_ratio['feature_type'] = 'pixie_ecm'
ecm_object_ratio['feature_type_detail'] = 'major_minor_ratio'
ecm_object_ratio['feature_type_detail_2'] = ''
ecm_object_ratio = ecm_object_ratio.rename(columns={'axis_ratio': 'value'})
ecm_object_ratio = ecm_object_ratio[['fov', 'value', 'feature_name', 'feature_name_unique',
                                'compartment', 'cell_pop', 'cell_pop_level', 'feature_type', 'feature_type_detail', 'feature_type_detail_2']]

fov_data.append(ecm_object_ratio)

# add ecm shape normalized difference
ecm_object_diff['feature_name'] = 'Pixie__major_minor_diff__' + ecm_object_diff['cluster']
ecm_object_diff['feature_name_unique'] = 'Pixie__major_minor_diff__' + ecm_object_diff['cluster']
ecm_object_diff['compartment'] = 'all'
ecm_object_diff['cell_pop'] = 'ecm'
ecm_object_diff['cell_pop_level'] = 'broad'
ecm_object_diff['feature_type'] = 'pixie_ecm'
ecm_object_diff['feature_type_detail'] = 'major_minor_diff'
ecm_object_diff['feature_type_detail_2'] = ''
ecm_object_diff = ecm_object_diff.rename(columns={'axis_diff_norm': 'value'})
ecm_object_diff = ecm_object_diff[['fov', 'value', 'feature_name', 'feature_name_unique',
                                'compartment', 'cell_pop', 'cell_pop_level', 'feature_type', 'feature_type_detail', 'feature_type_detail_2']]


# add fiber stats
fiber_df = fiber_df.rename(columns={'fiber_metric': 'feature_name'})

# drop rows with NAs or Inf
fiber_df = fiber_df.dropna()
fiber_df = fiber_df[~fiber_df.isin([np.nan, np.inf, -np.inf]).any(axis=1)]

fiber_df['feature_name_unique'] = fiber_df['feature_name']
fiber_df['compartment'] = 'all'
fiber_df['cell_pop'] = 'ecm'
fiber_df['cell_pop_level'] = 'broad'
fiber_df['feature_type'] = 'fiber'
fiber_df['feature_type_detail'] = fiber_df['feature_name']
fiber_df['feature_type_detail_2'] = ''
fiber_df = fiber_df[['fov', 'value', 'feature_name', 'feature_name_unique',
                        'compartment', 'cell_pop', 'cell_pop_level', 'feature_type', 'feature_type_detail', 'feature_type_detail_2']]
fov_data.append(fiber_df)

# add fiber tile stats
fiber_tile_df = fiber_tile_df.rename(columns={'fiber_metric': 'feature_name'})
fiber_tile_df['feature_name_unique'] = fiber_tile_df['feature_name']
fiber_tile_df['compartment'] = 'all'
fiber_tile_df['cell_pop'] = 'ecm'
fiber_tile_df['cell_pop_level'] = 'broad'
fiber_tile_df['feature_type'] = 'fiber'
fiber_tile_df['feature_type_detail'] = fiber_tile_df['feature_name']
fiber_tile_df['feature_type_detail_2'] = ''
fiber_tile_df = fiber_tile_df[['fov', 'value', 'feature_name', 'feature_name_unique',
                        'compartment', 'cell_pop', 'cell_pop_level', 'feature_type', 'feature_type_detail', 'feature_type_detail_2']]
fov_data.append(fiber_tile_df)

#
# Once individual files have been generated, we combine and standardize them
#

# compute z-scores for each feature
fov_data_df = pd.concat(fov_data)
fov_data_df = fov_data_df.rename(columns={'value': 'raw_value'})
fov_data_wide = fov_data_df.pivot(index='fov', columns='feature_name_unique', values='raw_value')
zscore_df = (fov_data_wide - fov_data_wide.mean()) / fov_data_wide.std()

# add z-scores to fov_data_df
zscore_df = zscore_df.reset_index()
zscore_df_long = pd.melt(zscore_df, id_vars='fov', var_name='feature_name_unique', value_name='normalized_value')
fov_data_df = pd.merge(fov_data_df, zscore_df_long, on=['fov', 'feature_name_unique'], how='left')

# add metadata
fov_data_df = pd.merge(fov_data_df, harmonized_metadata_df[['Tissue_ID', 'fov']], on='fov', how='left')

# rearrange columns
fov_data_df = fov_data_df[['Tissue_ID', 'fov', 'raw_value', 'normalized_value', 'feature_name', 'feature_name_unique',
                            'compartment', 'cell_pop', 'cell_pop_level', 'feature_type', 'feature_type_detail', 'feature_type_detail_2']]

fov_data_df.to_csv(os.path.join(analysis_dir, 'fov_features.csv'), index=False)


# create timepoint-level stats file
grouped = fov_data_df.groupby(['Tissue_ID', 'feature_name', 'feature_name_unique', 'compartment', 'cell_pop',
                               'cell_pop_level', 'feature_type']).agg({'raw_value': ['mean', 'std'],
                                                                       'normalized_value': ['mean', 'std']})
grouped.columns = ['raw_mean', 'raw_std', 'normalized_mean', 'normalized_std']
grouped = grouped.reset_index()

grouped.to_csv(os.path.join(analysis_dir, 'timepoint_features.csv'), index=False)

#
# The code below only needs to be run once to identify correlated features, after which
# it can be skipped and the output file can be read in
#

# filter FOV features based on correlation in compartments
#
# study_fovs = harmonized_metadata_df.loc[harmonized_metadata_df.Timepoint.isin(TIMEPOINT_NAMES), 'fov'].unique()
# fov_data_df = pd.read_csv(os.path.join(analysis_dir, 'fov_features.csv'))
#
# # filter out features that are highly correlated in compartments
# feature_names = fov_data_df.feature_name.unique()
# exclude_list = []
#
# # set minimum number of FOVs for compartment feature
# min_fovs = 100
#
# for feature_name in feature_names:
#     fov_data_feature = fov_data_df.loc[fov_data_df.feature_name == feature_name, :]
#     fov_data_feature = fov_data_feature.loc[fov_data_feature.fov.isin(study_fovs), :]
#
#     # get the compartments present for this feature
#     compartments = fov_data_feature.compartment.unique()
#
#     # if only one compartment, skip
#     if len(compartments) == 1:
#         continue
#
#     fov_data_wide = fov_data_feature.pivot(index='fov', columns='compartment', values='raw_value')
#
#     # filter out features that are nans or mostly zeros
#     for compartment in compartments:
#         nan_count = fov_data_wide[compartment].isna().sum()
#         zero_count = (fov_data_wide[compartment] == 0).sum()
#
#         if len(fov_data_wide) - nan_count - zero_count < min_fovs:
#             exclude_list.append(feature_name + '__' + compartment)
#             fov_data_wide = fov_data_wide.drop(columns=compartment)
#
#     # compute correlations
#     compartments = fov_data_wide.columns
#     compartments = compartments[compartments != 'all']
#
#     for compartment in compartments:
#         corr, _ = spearmanr(fov_data_wide['all'].values, fov_data_wide[compartment].values, nan_policy='omit')
#         if corr > 0.7:
#             exclude_list.append(feature_name + '__' + compartment)
#
# exclude_df = pd.DataFrame({'feature_name_unique': exclude_list})
# exclude_df.to_csv(os.path.join(intermediate_dir, 'post_processing', 'exclude_features_compartment_correlation.csv'), index=False)

# use pre-defined list of features to exclude
exclude_df = pd.read_csv(os.path.join(intermediate_dir, 'post_processing', 'exclude_features_compartment_correlation.csv'))

fov_data_df_filtered = fov_data_df.loc[~fov_data_df.feature_name_unique.isin(exclude_df.feature_name_unique.values), :]
fov_data_df_filtered.to_csv(os.path.join(analysis_dir, 'fov_features_filtered.csv'), index=False)

# group by timepoint
grouped = fov_data_df_filtered.groupby(['Tissue_ID', 'feature_name', 'feature_name_unique', 'compartment', 'cell_pop',
                                 'cell_pop_level', 'feature_type']).agg({'raw_value': ['mean', 'std'],
                                                                            'normalized_value': ['mean', 'std']})

grouped.columns = ['raw_mean', 'raw_std', 'normalized_mean', 'normalized_std']
grouped = grouped.reset_index()

grouped.to_csv(os.path.join(analysis_dir, 'timepoint_features_filtered.csv'), index=False)


# get feature metadata
feature_metadata = fov_data_df_filtered[['feature_name', 'feature_name_unique', 'compartment', 'cell_pop', 'cell_pop_level', 'feature_type', 'feature_type_detail', 'feature_type_detail_2']]
feature_metadata = feature_metadata.drop_duplicates()

feature_metadata.to_csv(os.path.join(analysis_dir, 'feature_metadata.csv'), index=False)
