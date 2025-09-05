import os
import pandas as pd
import matplotlib
import itertools
import shutil
import matplotlib.lines as mlines

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from venny4py.venny4py import venny4py
from matplotlib_venn import venn2
import anndata

BASE_DIR = "/Volumes/Shared/Noah Greenwald/TONIC_Cohort/"
SUPPLEMENTARY_FIG_DIR = os.path.join(BASE_DIR, "supplementary_figs")
REVIEW_FIG_DIR = os.path.join(BASE_DIR, "supplementary_figs/review_figures")
ANALYSIS_DIR = os.path.join(BASE_DIR, "analysis_files")


'''
## GENERATE SPACECAT FEATURES FOR NT DATA ##
## RUN SPACECAT ON NT DATA
NT_DIR = '/Volumes/Shared/Noah Greenwald/NTPublic'
harmonized_metadata = pd.read_csv(os.path.join(NT_DIR, 'analysis_files', 'harmonized_metadata.csv'))

## SPACECAT FEATURES ONLY
adata_processed = anndata.read_h5ad(os.path.join(NT_DIR, 'adata', 'adata_preprocessed.h5ad'))

features = SpaceCat(adata_processed, image_key='fov', seg_label_key='label', cell_area_key='area',
                    cluster_key=['cell_cluster', 'cell_cluster_broad'],
                    compartment_key='compartment', compartment_area_key='compartment_area')

mixing_df = pd.read_csv(os.path.join(NT_DIR, 'intermediate_files/spatial_analysis/mixing_score/homogeneous_mixing_scores-ratio5_count200.csv'))
cols = mixing_df.columns.tolist()
keep_cols = [col for col in cols if 'mixing_score' in col]
mixing_df = mixing_df[['fov'] + keep_cols]
mixing_df = mixing_df[mixing_df.fov.isin(adata_processed.obs.fov.unique())]

# specify cell type pairs to compute a ratio for
ratio_pairings = [('CD8T', 'CD4T'), ('CD4T', 'Treg'), ('CD8T', 'Treg')]
# specify addtional per cell and per image stats
per_cell_stats = [
    ['morphology', 'cell_cluster', ['cell_size']]
]
per_img_stats = [
    ['mixing_score', mixing_df],
]

import time
start = time.time()
print('Start SpaceCat')
# Generate features and save anndata
adata_processed = features.run_spacecat(functional_feature_level='cell_cluster', diversity_feature_level='cell_cluster',
                                        pixel_radius=25,
                                        specified_ratios_cluster_key='cell_cluster', specified_ratios=ratio_pairings,
                                        per_cell_stats=per_cell_stats, per_img_stats=per_img_stats,
                                        correlation_filtering_thresh=0.95)

end = time.time()
length = end - start
print(f'End SpaceCat - {length/60}')

folder = 'SpaceCat'
out_dir = os.path.join(NT_DIR, folder)
os.makedirs(out_dir, exist_ok=True)
adata_processed.write_h5ad(os.path.join(out_dir, 'adata_processed.h5ad'))

# Save finalized tables to csv
adata_processed.uns['combined_feature_data'].to_csv(os.path.join(out_dir, 'combined_feature_data.csv'), index=False)
adata_processed.uns['combined_feature_data_filtered'].to_csv(os.path.join(out_dir, 'combined_feature_data_filtered.csv'), index=False)
adata_processed.uns['feature_metadata'].to_csv(os.path.join(out_dir, 'feature_metadata.csv'), index=False)
adata_processed.uns['excluded_features'].to_csv(os.path.join(out_dir, 'excluded_features.csv'), index=False)

df = pd.read_csv(os.path.join(out_dir, 'combined_feature_data_filtered.csv'))
feature_data = df.merge(harmonized_metadata[['fov', 'Tissue_ID']], on='fov', how='left')
grouped = feature_data.groupby(['Tissue_ID', 'feature_name', 'feature_name_unique', 'compartment',
                                 'cell_pop_level', 'feature_type']).agg({'raw_value': ['mean', 'std'],
                                                                            'normalized_value': ['mean', 'std']})
grouped.columns = ['raw_mean', 'raw_std', 'normalized_mean', 'normalized_std']
grouped = grouped.reset_index()

output_dir = os.path.join(out_dir, 'output_files')
analysis_dir = os.path.join(out_dir, 'analysis_files')

os.makedirs(analysis_dir, exist_ok=True)
os.makedirs(output_dir, exist_ok=True)
os.makedirs(os.path.join(out_dir, 'prediction_model'), exist_ok=True)
grouped.to_csv(os.path.join(analysis_dir, 'timepoint_features_filtered.csv'), index=False)

timepoint_features = pd.read_csv(os.path.join(analysis_dir, f'timepoint_features_filtered.csv'))
tissue_treatment = harmonized_metadata[['Tissue_ID', 'Arm']].drop_duplicates()
timepoint_features_subset = timepoint_features.merge(tissue_treatment)

timepoint_features_C = timepoint_features_subset[timepoint_features_subset.Arm=='C']
timepoint_features_C = timepoint_features_C.drop(columns='Arm')
timepoint_features_C.to_csv(os.path.join(analysis_dir, 'timepoint_features_filtered_chemotherapy.csv'), index=False)
timepoint_features_C_I = timepoint_features_subset[timepoint_features_subset.Arm=='C&I']
timepoint_features_C_I = timepoint_features_C_I.drop(columns='Arm')
timepoint_features_C_I.to_csv(os.path.join(analysis_dir, 'timepoint_features_filtered_immunotherapy+chemotherapy.csv'), index=False)


## CALCULATE INTERACTION FEATURES
adata = anndata.read_h5ad(os.path.join(NT_DIR, 'SpaceCat', 'adata_processed.h5ad'))
neighborhood_mat = adata.obsm['neighbors_counts_cell_cluster_radius25']
neighborhood_mat = pd.concat([adata.obs.loc[:, ['fov', 'label', 'cell_cluster']], neighborhood_mat], axis=1)
interactions_df = neighborhood_mat[['fov', 'cell_cluster']]

# define epithelial and TME cell clusters
ep_cells = ['Epithelial_1', 'Epithelial_3', 'Epithelial_2']
tme_cells = ['APC', 'B', 'CD4T', 'CD8T', 'DC', 'Endothelial',
       'Fibroblast', 'M2_Mac', 'NK', 'Neutrophil', 'PDPN', 'Plasma',
       'Treg']
interactions_df['Epi'] = neighborhood_mat[ep_cells].sum(axis=1)
interactions_df['TME'] = neighborhood_mat[tme_cells].sum(axis=1)
interactions_df = interactions_df.groupby(by=['fov', 'cell_cluster']).sum(numeric_only=True).reset_index()

interactions_df = interactions_df.merge(adata.obs[['fov', 'cell_size']].groupby(by='fov').count().reset_index(), on='fov')
interactions_df = interactions_df.rename(columns={'cell_size': 'total_cells'})
interactions_df['Epi'] = interactions_df['Epi'] / interactions_df['total_cells']
interactions_df['TME'] = interactions_df['TME'] / interactions_df['total_cells']
interactions_df['isEpi'] = interactions_df['cell_cluster'].isin(ep_cells)
interactions_df['isEpi'] = ['Epi' if check else 'TME' for check in interactions_df['isEpi']]

# sum total homotypic and heterotypic interactions per cell type
interactions_df = pd.melt(interactions_df, id_vars=['fov', 'cell_cluster', 'isEpi'], value_vars=['Epi', 'TME'])
interactions_df['feature_type'] = [cell + 'Hom' if cell == interactions_df['variable'][i] else cell + 'Het' for i, cell in enumerate(interactions_df.isEpi)]
interactions_df['feature'] = interactions_df["cell_cluster"].astype(str) + '__' + interactions_df["feature_type"].astype(str)
interactions_df = interactions_df.pivot(index='fov', columns='feature', values='value').reset_index()

folder = 'SpaceCat_NT_combined'
out_dir = os.path.join(NT_DIR, folder)
os.makedirs(out_dir, exist_ok=True)
interactions_df.to_csv(os.path.join(out_dir, 'Epi_TME_interactions.csv'), index=False)


## SPACECAT + WANG INTERACTION FEATURES
adata_processed = anndata.read_h5ad(os.path.join(NT_DIR, 'adata', 'adata_preprocessed.h5ad'))

features = SpaceCat(adata_processed, image_key='fov', seg_label_key='label', cell_area_key='area',
                    cluster_key=['cell_cluster', 'cell_cluster_broad'],
                    compartment_key='compartment', compartment_area_key='compartment_area')

mixing_df = pd.read_csv(os.path.join(NT_DIR, 'intermediate_files/spatial_analysis/mixing_score/homogeneous_mixing_scores-ratio5_count200.csv'))
cols = mixing_df.columns.tolist()
keep_cols = [col for col in cols if 'mixing_score' in col]
mixing_df = mixing_df[['fov'] + keep_cols]
mixing_df = mixing_df[mixing_df.fov.isin(adata_processed.obs.fov.unique())]
cell_interactions = pd.read_csv(os.path.join(out_dir, 'Epi_TME_interactions.csv'))

# specify cell type pairs to compute a ratio for
ratio_pairings = [('CD8T', 'CD4T'), ('CD4T', 'Treg'), ('CD8T', 'Treg')]
# specify addtional per cell and per image stats
per_cell_stats = [
    ['morphology', 'cell_cluster', ['cell_size']]
]
per_img_stats = [
    ['mixing_score', mixing_df],
    ['cell_interactions', cell_interactions]
]

# Generate features and save anndata
adata_processed = features.run_spacecat(functional_feature_level='cell_cluster', diversity_feature_level='cell_cluster',
                                        pixel_radius=25,
                                        specified_ratios_cluster_key='cell_cluster', specified_ratios=ratio_pairings,
                                        per_cell_stats=per_cell_stats, per_img_stats=per_img_stats,
                                        correlation_filtering_thresh=0.95)
adata_processed.write_h5ad(os.path.join(out_dir, 'adata_processed.h5ad'))

# Save finalized tables to csv
adata_processed.uns['combined_feature_data'].to_csv(os.path.join(out_dir, 'combined_feature_data.csv'), index=False)
adata_processed.uns['combined_feature_data_filtered'].to_csv(os.path.join(out_dir, 'combined_feature_data_filtered.csv'), index=False)
adata_processed.uns['feature_metadata'].to_csv(os.path.join(out_dir, 'feature_metadata.csv'), index=False)
adata_processed.uns['excluded_features'].to_csv(os.path.join(out_dir, 'excluded_features.csv'), index=False)

df = pd.read_csv(os.path.join(out_dir, 'combined_feature_data_filtered.csv'))
feature_data = df.merge(harmonized_metadata[['fov', 'Tissue_ID']], on='fov', how='left')
grouped = feature_data.groupby(['Tissue_ID', 'feature_name', 'feature_name_unique', 'compartment',
                                 'cell_pop_level', 'feature_type']).agg({'raw_value': ['mean', 'std'],
                                                                            'normalized_value': ['mean', 'std']})
grouped.columns = ['raw_mean', 'raw_std', 'normalized_mean', 'normalized_std']
grouped = grouped.reset_index()

output_dir = os.path.join(out_dir, 'output_files')
analysis_dir = os.path.join(out_dir, 'analysis_files')

os.makedirs(analysis_dir, exist_ok=True)
os.makedirs(output_dir, exist_ok=True)
os.makedirs(os.path.join(out_dir, 'prediction_model'), exist_ok=True)
grouped.to_csv(os.path.join(analysis_dir, 'timepoint_features_filtered.csv'), index=False)

timepoint_features = pd.read_csv(os.path.join(analysis_dir, f'timepoint_features_filtered.csv'))
tissue_treatment = harmonized_metadata[['Tissue_ID', 'Arm']].drop_duplicates()
timepoint_features_subset = timepoint_features.merge(tissue_treatment)

timepoint_features_C = timepoint_features_subset[timepoint_features_subset.Arm=='C']
timepoint_features_C = timepoint_features_C.drop(columns='Arm')
timepoint_features_C.to_csv(os.path.join(analysis_dir, 'timepoint_features_filtered_chemotherapy.csv'), index=False)
timepoint_features_C_I = timepoint_features_subset[timepoint_features_subset.Arm=='C&I']
timepoint_features_C_I = timepoint_features_C_I.drop(columns='Arm')
timepoint_features_C_I.to_csv(os.path.join(analysis_dir, 'timepoint_features_filtered_immunotherapy+chemotherapy.csv'), index=False)


## 7_create_evolution_df.py converted
TIMEPOINT_NAMES = ['Baseline', 'On-treatment']

# define the timepoint pairs to use, these index into harmonized_metadata
TIMEPOINT_PAIRS = [
    ("Baseline", "On-treatment")
]

# define the timepoint columns as they appear in timepoint_df
TIMEPOINT_COLUMNS = [
    "Baseline__On-treatment",
]

for analysis_method, treatment in itertools.product(['SpaceCat', 'SpaceCat_NT_combined'], ['_chemotherapy', '_immunotherapy+chemotherapy']):
    analysis_dir = os.path.join(NT_DIR, analysis_method, 'analysis_files')
    harmonized_metadata = pd.read_csv(os.path.join(os.path.join(NT_DIR, 'analysis_files', 'harmonized_metadata.csv')))
    timepoint_features = pd.read_csv(os.path.join(analysis_dir, f'timepoint_features_filtered{treatment}.csv'))
    timepoint_features_agg = timepoint_features.merge(harmonized_metadata[['Tissue_ID', 'Timepoint', 'Patient_ID'] + TIMEPOINT_COLUMNS].drop_duplicates(), on='Tissue_ID', how='left')
    patient_metadata = pd.read_csv(os.path.join(os.path.join(NT_DIR, 'intermediate_files', 'metadata/NTPublic_data_per_patient.csv')))


    # add evolution features to get finalized features specified by timepoint
    combine_features(analysis_dir, harmonized_metadata, timepoint_features, timepoint_features_agg, 
                    patient_metadata, metadata_cols=['Patient_ID', 'pCR', 'Arm'], 
                    timepoint_columns=TIMEPOINT_COLUMNS, timepoint_names=TIMEPOINT_NAMES, timepoint_pairs = TIMEPOINT_PAIRS,
                    file_suffix=treatment)

    # generate  pvalues and feature ranking
    feature_metadata = pd.read_csv(os.path.join(NT_DIR, analysis_method, 'feature_metadata.csv'))
    combined_df = pd.read_csv(os.path.join(analysis_dir, f'timepoint_combined_features{treatment}.csv'))
    generate_feature_rankings(analysis_dir, combined_df, feature_metadata, pop_col='pCR', 
                              populations=['RD', 'pCR'], timepoint_names=TIMEPOINT_NAMES, file_suffix=treatment)
                              
## PREDICTION PREPROCESSING
for analysis_method in ['SpaceCat', 'SpaceCat_NT_combined']:
    sub_dir = os.path.join(NT_DIR, analysis_method)
    top_dir = os.path.join(sub_dir ,'prediction_model')
    analysis_dir = os.path.join(NT_DIR, analysis_method, 'analysis_files')

    os.makedirs(os.path.join(top_dir, 'chemotherapy'), exist_ok=True)
    os.makedirs(os.path.join(top_dir, 'immunotherapy+chemotherapy'), exist_ok=True)

    for treatment in ['_chemotherapy', '_immunotherapy+chemotherapy']:

        # read data
        prediction_dir = os.path.join(top_dir, treatment[1:])
        shutil.copy(os.path.join(analysis_dir, f'timepoint_combined_features{treatment}.csv'), os.path.join(prediction_dir, f'timepoint_combined_features{treatment}.csv'))
        
        df_feature = pd.read_csv(os.path.join(prediction_dir, f'timepoint_combined_features{treatment}.csv'))
        prediction_preprocessing(df_feature, prediction_dir)
        
        os.makedirs(os.path.join(prediction_dir, 'patient_outcomes'), exist_ok=True)
        
    ## RUN all_timepoints.R

'''

# plot top baseline features
NT_DIR = '/Volumes/Shared/Noah Greenwald/NTPublic'
combined_df = pd.read_csv(
    os.path.join(NT_DIR, 'SpaceCat/analysis_files/timepoint_combined_features_immunotherapy+chemotherapy.csv'))

for feature in combined_df.feature_name_unique.unique():
    if 'Epithelial' in feature:
        feature_new = feature.replace('Epithelial', 'Cancer')
        combined_df = combined_df.replace({feature: feature_new})

for timepoint in ['Baseline', 'On-treatment']:
    for feature, lims in zip(['Ki67+__Cancer_1', 'Ki67+__CD8T', 'Cancer_1__cell_cluster_density'],
                             [[0, 1], [0, 0.6], [0, 1]]):
        plot_df = combined_df.loc[(combined_df.feature_name_unique == feature) &
                                  (combined_df.Timepoint == timepoint), :]

        fig, ax = plt.subplots(1, 1, figsize=(2, 4))
        sns.stripplot(data=plot_df, x='pCR', y='raw_mean', order=['pCR', 'RD'], color='black', ax=ax)
        sns.boxplot(data=plot_df, x='pCR', y='raw_mean', order=['pCR', 'RD'], color='grey', ax=ax, showfliers=False,
                    width=0.3)
        ax.set_title(timepoint)
        ax.set_ylim(lims)
        sns.despine()
        plt.tight_layout()
        os.makedirs(os.path.join(SUPPLEMENTARY_FIG_DIR, 'ED_figure_10'), exist_ok=True)
        plt.savefig(os.path.join(SUPPLEMENTARY_FIG_DIR, 'ED_figure_10', '{}_{}.pdf'.format(feature, timepoint)))
        plt.close()


# WANG DATA COMPARISON
pred_chemo = pd.read_csv(os.path.join(NT_DIR, 'SpaceCat/prediction_model/chemotherapy/patient_outcomes/all_timepoints_results.csv'))
pred_chemo_immuno = pd.read_csv(os.path.join(NT_DIR, 'SpaceCat/prediction_model/immunotherapy+chemotherapy/patient_outcomes/all_timepoints_results.csv'))
pred_chemo = pred_chemo.rename(columns={'auc_baseline_list': 'baseline_C', 'auc_on_treatment_list': 'on_treatment_C',
                                        'auc_baseline__on_treatment_list': 'both_C'})
pred_chemo_immuno = pred_chemo_immuno.rename(columns={'auc_baseline_list': 'baseline_C&I', 'auc_on_treatment_list': 'on_treatment_C&I',
                                                      'auc_baseline__on_treatment_list': 'both_C&I'})
NT_preds = pd.concat([pred_chemo, pred_chemo_immuno], axis=1)
NT_preds = NT_preds[['baseline_C', 'baseline_C&I', 'on_treatment_C', 'on_treatment_C&I', 'both_C', 'both_C&I']]
NT_preds = NT_preds.rename(columns={'baseline_C': 'Baseline (C)', 'baseline_C&I': 'Baseline (C&I)', 'on_treatment_C': 'On-treatment (C)',
                                    'on_treatment_C&I': 'On-treatment (C&I)', 'both_C': 'Both (C)', 'both_C&I': 'Both (C&I)'})
NT_preds = NT_preds.melt()
NT_preds['cancer_revised'] = 1

pred_chemo = pd.read_csv(os.path.join(NT_DIR, 'SpaceCat_NT_combined/prediction_model/chemotherapy/patient_outcomes/all_timepoints_results.csv'))
pred_chemo_immuno = pd.read_csv(os.path.join(NT_DIR, 'SpaceCat_NT_combined/prediction_model/immunotherapy+chemotherapy/patient_outcomes/all_timepoints_results.csv'))
pred_chemo = pred_chemo.rename(columns={'auc_baseline_list': 'baseline_C', 'auc_on_treatment_list': 'on_treatment_C',
                                        'auc_baseline__on_treatment_list': 'both_C'})
pred_chemo_immuno = pred_chemo_immuno.rename(columns={'auc_baseline_list': 'baseline_C&I', 'auc_on_treatment_list': 'on_treatment_C&I',
                                                      'auc_baseline__on_treatment_list': 'both_C&I'})
combo_preds = pd.concat([pred_chemo, pred_chemo_immuno], axis=1)
combo_preds = combo_preds[['baseline_C', 'baseline_C&I', 'on_treatment_C', 'on_treatment_C&I', 'both_C', 'both_C&I']]
combo_preds = combo_preds.rename(columns={'baseline_C': 'Baseline (C)', 'baseline_C&I': 'Baseline (C&I)', 'on_treatment_C': 'On-treatment (C)',
                                          'on_treatment_C&I': 'On-treatment (C&I)', 'both_C': 'Both (C)', 'both_C&I': 'Both (C&I)'})
combo_preds = combo_preds.melt()
combo_preds['cancer_revised'] = 2

og_preds = pd.read_csv(os.path.join(NT_DIR, 'NT_preds.csv'))
og_preds = og_preds.replace('Base&On', 'Both')
og_preds['variable'] = og_preds['Timepoint'] + ' (' + og_preds['Arm'] + ')'
og_preds = og_preds[['Fold', 'LassoAUC', 'variable']]
og_preds = og_preds.pivot(index='Fold', columns='variable')
og_preds.columns = og_preds.columns.droplevel(0)
og_preds = og_preds.melt()
og_preds['cancer_revised'] = 0
all_preds = pd.concat([NT_preds, og_preds, combo_preds])

fig, ax = plt.subplots()
sns.boxplot(data=all_preds, x='variable', y='value', ax=ax, hue='cancer_revised',
            palette=sns.color_palette(["gold", "#1f77b4", "darkseagreen"]), showfliers=False)
sns.stripplot(data=all_preds, x='variable', y='value', ax=ax, hue='cancer_revised',
              palette=sns.color_palette(["gold", "#1f77b4", "darkseagreen"]), dodge=True)
fig.set_figheight(4)
fig.set_figwidth(8)
plt.xticks(rotation=45)
plt.title('Wang et al. dataset predictions')
plt.ylabel('AUC')
plt.xlabel('')
plt.ylim((0, 1))

# Add the custom legend
yellow_line = mlines.Line2D([], [], color="gold", marker="o", label="Wang et al. predictions", linestyle='None')
blue_line = mlines.Line2D([], [], color="#1f77b4", marker="o", label="SpaceCat predictions", linestyle='None')
green_line = mlines.Line2D([], [], color="darkseagreen", marker="o", label="SpaceCat & Wang et al. features", linestyle='None')
plt.legend(handles=[yellow_line, blue_line, green_line], loc='lower right')
sns.despine()
plt.savefig(os.path.join(SUPPLEMENTARY_FIG_DIR, 'ED_figure_10f.pdf'), bbox_inches='tight', dpi=300)

# TONIC DATA COMPARISON
SpaceCat_dir = os.path.join(BASE_DIR, 'TONIC_SpaceCat')
preds = pd.read_csv(os.path.join(SpaceCat_dir, 'SpaceCat/prediction_model/patient_outcomes/all_timepoints_results_MIBI.csv'))
preds = preds[['auc_primary_list', 'auc_baseline_list', 'auc_post_induction_list', 'auc_on_nivo_list']]
preds = preds.rename(columns={'auc_primary_list': 'Primary', 'auc_baseline_list': 'Baseline',
                              'auc_post_induction_list': 'Pre nivo', 'auc_on_nivo_list': 'On nivo'})
preds = preds.melt()
preds['Analysis'] = 0

adj_preds = pd.read_csv(os.path.join(SpaceCat_dir, 'SpaceCat_NT_combined/prediction_model/patient_outcomes/all_timepoints_results_MIBI.csv'))
adj_preds = adj_preds[['auc_primary_list', 'auc_baseline_list', 'auc_post_induction_list', 'auc_on_nivo_list']]
adj_preds = adj_preds.rename(columns={'auc_primary_list': 'Primary', 'auc_baseline_list': 'Baseline',
                                      'auc_post_induction_list': 'Pre nivo', 'auc_on_nivo_list': 'On nivo'})
adj_preds = adj_preds.melt()
adj_preds['Analysis'] = 2

nt_feats_preds = pd.read_csv(os.path.join(SpaceCat_dir, 'NT_features_only/prediction_model/patient_outcomes/all_timepoints_results_MIBI.csv'))
nt_feats_preds = nt_feats_preds[['auc_primary_list', 'auc_baseline_list', 'auc_post_induction_list', 'auc_on_nivo_list']]
nt_feats_preds = nt_feats_preds.rename(columns={'auc_primary_list': 'Primary', 'auc_baseline_list': 'Baseline',
                                                'auc_post_induction_list': 'Pre nivo', 'auc_on_nivo_list': 'On nivo'})
nt_feats_preds = nt_feats_preds.melt()
nt_feats_preds['Analysis'] = 1
all_preds = pd.concat([preds, adj_preds, nt_feats_preds])

fig, ax = plt.subplots()
sns.boxplot(data=all_preds, x='variable', y='value', ax=ax, hue='Analysis',
            palette=sns.color_palette(["#1f77b4", 'gold', "darkseagreen"]), showfliers=False)
sns.stripplot(data=all_preds, x='variable', y='value', ax=ax, hue='Analysis',
              palette=sns.color_palette(["#1f77b4", 'gold', "darkseagreen"]), dodge=True)

fig.set_figheight(4)
fig.set_figwidth(8)
plt.xticks(rotation=45)
plt.title('TONIC dataset prediction')
plt.ylabel('AUC')
plt.xlabel('')
plt.ylim((0, 1))
plt.legend(loc='lower right').set_title('')
sns.despine()

blue_line = mlines.Line2D([], [], color="#1f77b4", marker="o", label="SpaceCat features", linestyle='None')
yellow_line = mlines.Line2D([], [], color="gold", marker="o", label="Wang et al. features", linestyle='None')
green_line = mlines.Line2D([], [], color="darkseagreen", marker="o", label="SpaceCat & Wang et al. features", linestyle='None')
plt.legend(handles=[blue_line, yellow_line, green_line], loc='lower right')
plt.savefig(os.path.join(SUPPLEMENTARY_FIG_DIR, 'ED_figure_10g.pdf'), bbox_inches='tight', dpi=300)


# NT feature enrichment on TONIC data
harmonized_metadata = pd.read_csv(os.path.join(ANALYSIS_DIR, 'harmonized_metadata.csv'))
ranked_features_all = pd.read_csv(os.path.join(BASE_DIR, 'TONIC_SpaceCat/NT_features_only/analysis_files/feature_ranking.csv'))
ranked_features = ranked_features_all.loc[ranked_features_all.comparison.isin(['primary', 'baseline', 'pre_nivo', 'on_nivo'])]
top_features = ranked_features.loc[ranked_features.comparison.isin(['primary', 'baseline', 'pre_nivo', 'on_nivo']), :]
top_features = top_features.iloc[:100, :]

# summarize distribution of top features
top_features_by_comparison = top_features[['feature_name_unique', 'comparison']].groupby('comparison').count().reset_index()
top_features_by_comparison.columns = ['comparison', 'num_features']
top_features_by_comparison = top_features_by_comparison.sort_values('num_features', ascending=False)

fig, ax = plt.subplots(figsize=(4, 4))
sns.barplot(data=top_features_by_comparison, x='comparison', y='num_features', color='grey', ax=ax)
plt.xticks(rotation=90)
plt.tight_layout()
sns.despine()
plt.savefig(os.path.join(SUPPLEMENTARY_FIG_DIR, 'ED_figure_10j.pdf'))
plt.close()

## 2.8 / 4.8  Pre-treatment and On-treatment NT vs TONIC comparisons ##
# Original NT features
file_path = os.path.join(NT_DIR, 'data/41586_2023_6498_MOESM3_ESM.xlsx')
NT_features = pd.read_excel(file_path, sheet_name=None)
cell_table = pd.read_csv(os.path.join(NT_DIR, 'analysis_files/cell_table.csv'))
cell_table = cell_table.replace({'{': '', '}': ''}, regex=True)
cell_table['cell_meta_cluster'] = [cell.replace('^', '')for cell in cell_table['cell_meta_cluster']]
cell_dict = dict(zip(cell_table.cell_meta_cluster, cell_table.cell_cluster))

density_pvals = NT_features['Table 8 Densities p values'][['Time point', 'Cell phenotype', 'p.value', 'Arm']]
density_pvals = density_pvals.rename(columns={'Cell phenotype': 'feature_name_unique'})
density_pvals = density_pvals.replace(cell_dict)
density_pvals['feature_name_unique'] = density_pvals['feature_name_unique'] + '__cell_cluster_density'

ki67_pvals = NT_features['Table 11 Ki67 p values'][['Time point', 'Cell phenotype', 'p.value', 'Unnamed: 10']]
ki67_pvals = ki67_pvals.rename(columns={'Unnamed: 10': 'Arm'})
ki67_pvals = ki67_pvals.rename(columns={'Cell phenotype': 'feature_name_unique'})
ki67_pvals = ki67_pvals.replace(cell_dict)
ki67_pvals['feature_name_unique'] = 'Ki67+__' + ki67_pvals['feature_name_unique']

interaction_pvals = NT_features['Table 10 Interaction p values'][['Time point', 'from', 'to cell phenotype', 'p.value', 'Arm']]
ep_cells = cell_table[cell_table.isEpithelial == 1].cell_meta_cluster.unique()
tme_cells = cell_table[cell_table.isEpithelial == 0].cell_meta_cluster.unique()
interaction_pvals['from'] = interaction_pvals['from'].replace('Epithelial', 'Epi')
interaction_pvals['isEpi'] = ['Epi' if cell in ep_cells else 'TME' for cell in interaction_pvals['to cell phenotype']]
interaction_pvals['interaction'] = ['Hom' if f == Epi else 'Het' for f, Epi in zip(interaction_pvals['from'], interaction_pvals['isEpi'])]
interaction_pvals = interaction_pvals.replace(cell_dict)
interaction_pvals['feature_name_unique'] = interaction_pvals['to cell phenotype'] + '__' + interaction_pvals['from'] + interaction_pvals['interaction']
interaction_pvals = interaction_pvals[['Time point', 'feature_name_unique', 'p.value', 'Arm']]

features = pd.concat([density_pvals, ki67_pvals, interaction_pvals])
features = features[features['p.value']<=0.05]
sig_features = features.replace({'Cancer_1__TMEHet': 'Cancer_Immune_mixing_score',
                                 'Cancer_2__TMEHet': 'Cancer_Immune_mixing_score',
                                 'Cancer_3__TMEHet': 'Cancer_Immune_mixing_score',
                                 'CD4T__TMEHom': 'Structural_T_mixing_score',
                                 'CD8T__TMEHom': 'Structural_T_mixing_score',
                                 'Treg__TMEHom': 'Structural_T_mixing_score'})
pre_treatment_features = sig_features[sig_features['Time point'] == 'Baseline']
on_treatment_features = sig_features[sig_features['Time point'] == 'On-treatment']

tonic_features = pd.read_csv(os.path.join(ANALYSIS_DIR, 'feature_ranking.csv'))
tonic_features = tonic_features[tonic_features['fdr_pval'] <= 0.05]
tonic__sig_features = tonic_features[tonic_features.compartment == 'all']
tonic__sig_features = tonic__sig_features[~tonic__sig_features.feature_name_unique.str.contains('core')]
tonic__sig_features = tonic__sig_features[~tonic__sig_features.feature_name_unique.str.contains('border')]
tonic_pre_treatment_features = tonic__sig_features[tonic__sig_features.comparison.isin(['baseline', 'pre_nivo'])]
tonic_pre_treatment_features = tonic_pre_treatment_features[['feature_name_unique', 'pval', 'comparison']]
tonic_on_treatment_features = tonic__sig_features[tonic__sig_features.comparison == 'on_nivo']
tonic_on_treatment_features = tonic_on_treatment_features[['feature_name_unique', 'pval', 'comparison']]

NT_feats = set(pre_treatment_features.feature_name_unique.unique())
TONIC_feats = set(tonic_pre_treatment_features.feature_name_unique.unique())
sets = {'Wang et al.': NT_feats, 'TONIC': TONIC_feats}
venny4py(sets=sets, colors="yb")
plt.title("Pre-treatment Features")
plt.savefig(os.path.join(SUPPLEMENTARY_FIG_DIR, 'ED_figure_10a.pdf'), bbox_inches='tight', dpi=300)

NT_feats = set(on_treatment_features.feature_name_unique.unique())
TONIC_feats = set(tonic_on_treatment_features.feature_name_unique.unique())
sets = {'Wang et al.': NT_feats, 'TONIC': TONIC_feats}
venny4py(sets=sets, colors="yb")
plt.title("On-treatment Features")
plt.savefig(os.path.join(SUPPLEMENTARY_FIG_DIR, 'ED_figure_10h.pdf'), bbox_inches='tight', dpi=300)

# compare SpaceCat features
NT_features = pd.read_csv(os.path.join(NT_DIR, 'SpaceCat/analysis_files/feature_ranking_immunotherapy+chemotherapy.csv'))
NT_features = NT_features[NT_features['fdr_pval'] <= 0.05]
NT__sig_features = NT_features[NT_features.compartment == 'all']
for feature in NT__sig_features.feature_name_unique.unique():
    if 'Epithelial' in feature:
        feature_new = feature.replace('Epithelial', 'Cancer')
        NT__sig_features = NT__sig_features.replace({feature: feature_new})
NT__sig_features = NT__sig_features[~NT__sig_features.feature_name_unique.str.contains('core')]
NT__sig_features = NT__sig_features[~NT__sig_features.feature_name_unique.str.contains('border')]
NT_pre_treatment_features = NT__sig_features[NT__sig_features.comparison == 'Baseline']
NT_pre_treatment_features = NT_pre_treatment_features[['feature_name_unique', 'pval', 'comparison']]
NT_on_treatment_features = NT__sig_features[NT__sig_features.comparison == 'On-treatment']
NT_on_treatment_features = NT_on_treatment_features[['feature_name_unique', 'pval', 'comparison']]

NT_feats = set(NT_pre_treatment_features.feature_name_unique.unique())
TONIC_feats = set(tonic_pre_treatment_features.feature_name_unique.unique())
sets = {'Wang et al.': NT_feats, 'TONIC': TONIC_feats}
venny4py(sets=sets, colors="yb")
plt.title("Pre-treatment Features")
plt.savefig(os.path.join(SUPPLEMENTARY_FIG_DIR, 'ED_figure_10b.pdf'), bbox_inches='tight', dpi=300)

NT_feats = set(NT_on_treatment_features.feature_name_unique.unique())
TONIC_feats = set(tonic_on_treatment_features.feature_name_unique.unique())
sets = {'Wang et al.': NT_feats, 'TONIC': TONIC_feats}
venny4py(sets=sets, colors="yb")
plt.title("On-treatment Features")
plt.savefig(os.path.join(SUPPLEMENTARY_FIG_DIR, 'ED_figure_10i.pdf'), bbox_inches='tight', dpi=300)


colors_dict = {'original':'#C9C9C9', 'updated': '#78CE8B'}

predicted_features_new = pd.read_csv(os.path.join(REVIEW_FIG_DIR, 'Cancer_reclustering', 'prediction_model', 'all_timepoints_results_MIBI.csv'))
predicted_features_new.columns = ['on nivo', 'baseline', 'pre nivo', 'primary']
predicted_features_new = predicted_features_new.loc[:, ['primary', 'baseline', 'pre nivo', 'on nivo']]
predicted_features_new = predicted_features_new.melt()
predicted_features_new['annotation'] = 'updated'

predicted_features_old = pd.read_csv(os.path.join(BASE_DIR, 'prediction_model', 'patient_outcomes', 'all_timepoints_results_MIBI.csv'))
predicted_features_old.columns = ['on nivo', 'baseline', 'pre nivo', 'primary']
predicted_features_old = predicted_features_old.loc[:, ['primary', 'baseline', 'pre nivo', 'on nivo']]
predicted_features_old = predicted_features_old.melt()
predicted_features_old['annotation'] = 'original'
predicted_features_joint = pd.concat([predicted_features_old, predicted_features_new], axis=0)

_, axes = plt.subplots(1, 1, figsize=(4.5, 4.5), gridspec_kw={'hspace': 0.65, 'wspace': 0.3, 'bottom': 0.15})
g = sns.boxplot(x='variable', y='value', hue='annotation', data=predicted_features_joint, linewidth=1, fliersize=0, width=0.6, ax=axes, palette=colors_dict)
g = sns.stripplot(x='variable', y='value', hue='annotation', data=predicted_features_joint, linewidth=0.8,  size=5, edgecolor="black", ax=axes, palette=colors_dict, legend=False, dodge=True, jitter=True)

g.tick_params(labelsize=10)
g.set_xlabel('Timepoint', fontsize=10)
g.set_ylabel('AUC', fontsize=10)
g.set_ylim(0.4, 1.0)
plt.savefig(os.path.join('auc.pdf'), bbox_inches='tight')
plt.savefig(os.path.join(SUPPLEMENTARY_FIG_DIR, 'ED_figure_10l.pdf'), bbox_inches='tight', dpi=300)

feature_ranking = pd.read_csv(os.path.join(REVIEW_FIG_DIR, 'Cancer_reclustering', 'SpaceCat', 'feature_ranking.csv'))
feature_ranking_df = feature_ranking[np.isin(feature_ranking['comparison'], ['primary', 'baseline', 'pre_nivo', 'on_nivo'])]
feature_ranking_df['feature_name_unique'] = [i.replace('cell_cluster_revised', 'cell_cluster') for i in feature_ranking_df['feature_name_unique'] ]
joint_features_new = feature_ranking_df['feature_name_unique'] + '_time_' + feature_ranking_df['comparison']
joint_features_new = joint_features_new[:100].values
joint_features_new = [i.replace('cell_cluster_revised', 'cell_cluster') for i in joint_features_new]

feature_ranking_old = pd.read_csv(os.path.join(REVIEW_FIG_DIR, 'SpaceCat_original', '', 'feature_ranking.csv'))
feature_ranking_df_old = feature_ranking_old[np.isin(feature_ranking_old['comparison'], ['primary', 'baseline', 'pre_nivo', 'on_nivo'])]
joint_features_old = feature_ranking_df_old['feature_name_unique'] + '_time_' + feature_ranking_df_old['comparison']
joint_features_old = joint_features_old[:100].values

new_not_old = list(set(joint_features_new).difference(set(joint_features_old)))
old_not_new = list(set(joint_features_old).difference(set(joint_features_new)))

original_set = set(joint_features_old)
updated_set = set(joint_features_new)
plt.figure(figsize=(4, 4))

v = venn2([original_set, updated_set], set_labels=['original', 'updated'], set_colors=['#C9C9C9', '#78CE8B'])

for region_id in ('10', '01', '11'):
    label = v.get_label_by_id(region_id)
    if label is not None:
        label.set_fontsize(14)

for label in v.set_labels:
    if label is not None:
        label.set_fontsize(16)

plt.savefig(os.path.join(SUPPLEMENTARY_FIG_DIR, 'ED_figure_10k.pdf'), bbox_inches='tight')
