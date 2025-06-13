import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

import os
import pandas as pd
import numpy as np

from sklearn import metrics
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold


BASE_DIR = "/Volumes/Shared/Noah Greenwald/TONIC_Cohort/"
ANALYSIS_DIR = os.path.join(BASE_DIR, 'analysis_files')

# supplementary tables
save_dir = os.path.join(BASE_DIR, 'supplementary_figs/supplementary_tables')
if not os.path.exists(save_dir):
    os.makedirs(save_dir)


# sample summary
harmonized_metadata = pd.read_csv(os.path.join(BASE_DIR, 'analysis_files/harmonized_metadata.csv'))
wes_metadata = pd.read_csv(os.path.join(BASE_DIR, 'sequencing_data/preprocessing/TONIC_WES_meta_table.tsv'), sep='\t')
rna_metadata = pd.read_csv(os.path.join(BASE_DIR, 'sequencing_data/preprocessing/TONIC_tissue_rna_id.tsv'), sep='\t')
rna_metadata = rna_metadata.merge(harmonized_metadata[['Patient_ID', 'Tissue_ID', 'Timepoint']].drop_duplicates(), on='Tissue_ID', how='left')
clinical_data = pd.read_csv(os.path.join(BASE_DIR, 'intermediate_files/metadata/patient_clinical_data.csv'))

harmonized_metadata = harmonized_metadata.loc[harmonized_metadata.MIBI_data_generated, :]
harmonized_metadata = harmonized_metadata[harmonized_metadata.Timepoint.isin(['primary', 'baseline', 'pre_nivo', 'on_nivo'])]
harmonized_metadata = harmonized_metadata[['Patient_ID', 'Timepoint']].merge(clinical_data, on='Patient_ID')
harmonized_metadata = harmonized_metadata[harmonized_metadata.Clinical_benefit.isin(['Yes', 'No'])]
wes_metadata = wes_metadata.rename(columns={'Individual.ID': 'Patient_ID', 'timepoint': 'Timepoint'})
# wes_metadata = wes_metadata[['Patient_ID', 'Timepoint']].merge(clinical_data, on='Patient_ID').drop_duplicates()
wes_metadata = wes_metadata[wes_metadata.Clinical_benefit.isin(['Yes', 'No'])]
rna_metadata = rna_metadata[['Patient_ID', 'Timepoint']].merge(clinical_data, on='Patient_ID').drop_duplicates()
rna_metadata = rna_metadata[rna_metadata.Clinical_benefit.isin(['Yes', 'No'])]

modality = ['MIBI'] * 4 + ['RNA'] * 3 + ['DNA'] * 1
timepoint = ['primary', 'baseline', 'pre_nivo', 'on_nivo'] + ['baseline', 'pre_nivo', 'on_nivo'] + ['baseline']

sample_summary_df = pd.DataFrame({'modality': modality, 'timepoint': timepoint, 'sample_num': [0] * 8, 'patient_num': [0] * 8, 'responder_num': [0] * 8, 'nonresponder_num': [0] * 8})

# populate dataframe
for idx, row in sample_summary_df.iterrows():
    if row.modality == 'MIBI':
        sample_summary_df.loc[idx, 'sample_num'] = len(harmonized_metadata.loc[harmonized_metadata.Timepoint == row.timepoint, :])
        sample_summary_df.loc[idx, 'patient_num'] = len(harmonized_metadata.loc[harmonized_metadata.Timepoint == row.timepoint, 'Patient_ID'].unique())
        sample_summary_df.loc[idx, 'responder_num'] = len(harmonized_metadata.loc[np.logical_and(harmonized_metadata.Timepoint == row.timepoint, harmonized_metadata.Clinical_benefit == 'Yes')])
        sample_summary_df.loc[idx, 'nonresponder_num'] = len(harmonized_metadata.loc[np.logical_and(harmonized_metadata.Timepoint == row.timepoint, harmonized_metadata.Clinical_benefit == 'No')])
    elif row.modality == 'RNA':
        sample_summary_df.loc[idx, 'sample_num'] = len(rna_metadata.loc[rna_metadata.Timepoint == row.timepoint, :])
        sample_summary_df.loc[idx, 'patient_num'] = len(rna_metadata.loc[rna_metadata.Timepoint == row.timepoint, 'Patient_ID'].unique())
        sample_summary_df.loc[idx, 'responder_num'] = len(rna_metadata.loc[np.logical_and(rna_metadata.Timepoint == row.timepoint, rna_metadata.Clinical_benefit == 'Yes')])
        sample_summary_df.loc[idx, 'nonresponder_num'] = len(rna_metadata.loc[np.logical_and(rna_metadata.Timepoint == row.timepoint, rna_metadata.Clinical_benefit == 'No')])
    elif row.modality == 'DNA':
        sample_summary_df.loc[idx, 'sample_num'] = len(wes_metadata.loc[wes_metadata.Timepoint == row.timepoint, :])
        sample_summary_df.loc[idx, 'patient_num'] = len(wes_metadata.loc[wes_metadata.Timepoint == row.timepoint, 'Patient_ID'].unique())
        sample_summary_df.loc[idx, 'responder_num'] = len(wes_metadata.loc[np.logical_and(wes_metadata.Timepoint == row.timepoint, wes_metadata.Clinical_benefit == 'Yes')])
        sample_summary_df.loc[idx, 'nonresponder_num'] = len(wes_metadata.loc[np.logical_and(wes_metadata.Timepoint == row.timepoint, wes_metadata.Clinical_benefit == 'No')])

sample_summary_df.to_csv(os.path.join(save_dir, 'Supplementary_Table_3.csv'), index=False)

# FOV counts per patient and timepoint
clinical_data = pd.read_csv(os.path.join(BASE_DIR, 'intermediate_files/metadata/patient_clinical_data.csv'))

harmonized_metadata = pd.read_csv(os.path.join(ANALYSIS_DIR, 'harmonized_metadata.csv'))
harmonized_metadata = harmonized_metadata[harmonized_metadata.Timepoint.isin(['primary', 'baseline', 'pre_nivo', 'on_nivo'])]
harmonized_metadata = harmonized_metadata[['Patient_ID', 'Timepoint', 'MIBI_data_generated', 'fov', 'rna_seq_sample_id']].merge(clinical_data, on='Patient_ID')
harmonized_metadata = harmonized_metadata[harmonized_metadata.Clinical_benefit.isin(['Yes', 'No'])]
mibi_metadata = harmonized_metadata[harmonized_metadata.MIBI_data_generated]

wes_metadata = pd.read_csv(os.path.join(BASE_DIR, 'sequencing_data/preprocessing/TONIC_WES_meta_table.tsv'), sep='\t')
wes_metadata = wes_metadata[wes_metadata.Clinical_benefit.isin(['Yes', 'No'])]
dna_counts = wes_metadata[['Individual.ID', 'timepoint', 'Experiment.System.ID']].drop_duplicates().groupby(['Individual.ID', 'timepoint'])['Experiment.System.ID'].count().unstack(fill_value=0).stack().reset_index()
dna_counts = dna_counts.rename(columns={0: 'DNA samples', 'Individual.ID': 'Patient_ID', 'timepoint': 'Timepoint'})

mibi_counts = mibi_metadata[['Patient_ID', 'Timepoint', 'fov']].groupby(['Patient_ID', 'Timepoint'])['fov'].count().unstack(fill_value=0).stack().reset_index()
mibi_counts = mibi_counts.rename(columns={0: 'MIBI fovs'})

rna_counts = harmonized_metadata[['Patient_ID', 'Timepoint', 'rna_seq_sample_id']].drop_duplicates().groupby(['Patient_ID', 'Timepoint'])['rna_seq_sample_id'].count().unstack(fill_value=0).stack().reset_index()
rna_counts = rna_counts.rename(columns={0: 'RNA samples'})

all_counts = mibi_counts.merge(rna_counts, on=['Patient_ID', 'Timepoint'], how='outer').fillna(0)
all_counts = all_counts.merge(dna_counts, on=['Patient_ID', 'Timepoint'], how='outer').fillna(0)
for col in ['MIBI fovs', 'RNA samples', 'DNA samples']:
    all_counts[col] = all_counts[col].astype(int)
all_counts.sort_values(by=['Patient_ID', 'Timepoint']).to_csv(os.path.join(save_dir, 'Supplementary_Table_4.csv'), index=False)

# feature metadata
feature_metadata = pd.read_csv(os.path.join(BASE_DIR, 'analysis_files/feature_metadata.csv'))
feature_metadata.columns = ['Feature name', 'Feature name including compartment', 'Compartment the feature is calculated in',
                            'Level of clustering granularity for cell types', 'Type of feature']

correlation_feature_order = pd.read_csv(os.path.join(BASE_DIR, 'supplementary_figs/review_figures/Correlation clustermap/clustermap_feature_order.csv'))
feature_metadata = feature_metadata.merge(correlation_feature_order, on='Feature name including compartment')
feature_metadata.to_csv(os.path.join(save_dir, 'Supplementary_Table_5.csv'), index=False)

# get overlap between static and evolution top features
ranked_features = pd.read_csv(os.path.join(BASE_DIR, 'analysis_files/feature_ranking.csv'))

overlap_type_dict = {'global': [['primary', 'baseline', 'pre_nivo', 'on_nivo'],
                                ['primary__baseline', 'baseline__pre_nivo', 'baseline__on_nivo', 'pre_nivo__on_nivo']],
                     'primary': [['primary'], ['primary__baseline']],
                     'baseline': [['baseline'], ['primary__baseline', 'baseline__pre_nivo', 'baseline__on_nivo']],
                     'pre_nivo': [['pre_nivo'], ['baseline__pre_nivo', 'pre_nivo__on_nivo']],
                     'on_nivo': [['on_nivo'], ['baseline__on_nivo', 'pre_nivo__on_nivo']]}
overlap_results = {}
for overlap_type, comparisons in overlap_type_dict.items():
    static_comparisons, evolution_comparisons = comparisons

    overlap_top_features = ranked_features.copy()
    overlap_top_features = overlap_top_features.loc[overlap_top_features.comparison.isin(static_comparisons + evolution_comparisons)]
    overlap_top_features.loc[overlap_top_features.comparison.isin(static_comparisons), 'comparison'] = 'static'
    overlap_top_features.loc[overlap_top_features.comparison.isin(evolution_comparisons), 'comparison'] = 'evolution'
    overlap_top_features = overlap_top_features[['feature_name_unique', 'comparison']].drop_duplicates()
    overlap_top_features = overlap_top_features.iloc[:100, :]
    static_ids = overlap_top_features.loc[
        overlap_top_features.comparison == 'static', 'feature_name_unique'].unique()
    evolution_ids = overlap_top_features.loc[
        overlap_top_features.comparison == 'evolution', 'feature_name_unique'].unique()

    overlap_results[overlap_type] = {'static_ids': static_ids, 'evolution_ids': evolution_ids}

# get counts of features in each category
static_ids = overlap_results['global']['static_ids']
evolution_ids = overlap_results['global']['evolution_ids']
overlap_features = list(set(static_ids).intersection(set(evolution_ids)))
evolution_ids = [feature for feature in evolution_ids if feature not in overlap_features]
static_ids = [feature for feature in static_ids if feature not in overlap_features]

static_features = pd.DataFrame({'Feature name': static_ids, 'Feature type': 'static'})
evolution_features = pd.DataFrame({'Feature name': evolution_ids, 'Feature type': 'evolution'})
shared_features = pd.DataFrame({'Feature name': overlap_features, 'Feature type': 'shared'})
all_features = pd.concat([static_features, evolution_features, shared_features])
all_features.to_csv(os.path.join(save_dir, 'Supplementary_Table_7.csv'), index=False)

# NT pre-treatment feature comparison
tonic_features = pd.read_csv(os.path.join(BASE_DIR, 'analysis_files/feature_ranking.csv'))
tonic_features = tonic_features[tonic_features['fdr_pval'] <= 0.05]
tonic__sig_features = tonic_features[tonic_features.compartment == 'all']
tonic__sig_features = tonic__sig_features[~tonic__sig_features.feature_name_unique.str.contains('core')]
tonic__sig_features = tonic__sig_features[~tonic__sig_features.feature_name_unique.str.contains('border')]
tonic_pre_treatment_features = tonic__sig_features[tonic__sig_features.comparison.isin(['baseline', 'pre_nivo'])]
tonic_pre_treatment_features = tonic_pre_treatment_features[['feature_name_unique', 'pval', 'comparison']]
tonic_on_treatment_features = tonic__sig_features[tonic__sig_features.comparison == 'on_nivo']
tonic_on_treatment_features = tonic_on_treatment_features[['feature_name_unique', 'pval', 'comparison']]

NT_DIR = '/Volumes/Shared/Noah Greenwald/NTPublic'
file_path = os.path.join(NT_DIR, 'data/41586_2023_6498_MOESM3_ESM.xlsx')
NT_features = pd.read_excel(file_path, sheet_name=None)
cell_table = pd.read_csv(os.path.join(NT_DIR, 'analysis_files/cell_table.csv'))
cell_table = cell_table.replace({'{': '', '}': ''}, regex=True)
cell_table['cell_meta_cluster'] = [cell.replace('^', '') for cell in cell_table['cell_meta_cluster']]

cell_dict = dict(zip(cell_table.cell_meta_cluster, cell_table.cell_cluster))
cell_dict['CD8+PD1+TEx'] = 'CD8T'
cell_dict['M2 Mac'] = 'CD163_Mac'
for key in cell_dict.keys():
    if 'Epithelial' in cell_dict[key]:
        cell_dict[key] = cell_dict[key].replace('Epithelial', 'Cancer')

density_pvals = NT_features['Table 8 Densities p values'][['Time point', 'Cell phenotype', 'p.value', 'Arm']]
interaction_pvals = NT_features['Table 10 Interaction p values'][
    ['Time point', 'from', 'to cell phenotype', 'p.value', 'Arm']]
ki67_pvals = NT_features['Table 11 Ki67 p values'][['Time point', 'Cell phenotype', 'p.value', 'Unnamed: 10']]
ki67_pvals = ki67_pvals.rename(columns={'Unnamed: 10': 'Arm'})

density_pvals = density_pvals.rename(columns={'Cell phenotype': 'feature_name_unique'})
density_pvals = density_pvals.replace(cell_dict)
density_pvals['feature_name_unique'] = density_pvals['feature_name_unique'] + '__cell_cluster_density'
ki67_pvals = ki67_pvals.rename(columns={'Cell phenotype': 'feature_name_unique'})
ki67_pvals = ki67_pvals.replace(cell_dict)
ki67_pvals['feature_name_unique'] = 'Ki67+__' + ki67_pvals['feature_name_unique']

ep_cells = cell_table[cell_table.isEpithelial == 1].cell_meta_cluster.unique()
tme_cells = cell_table[cell_table.isEpithelial == 0].cell_meta_cluster.unique()
interaction_pvals['from'] = interaction_pvals['from'].replace('Epithelial', 'Epi')
interaction_pvals['isEpi'] = ['Epi' if cell in ep_cells else 'TME' for cell in interaction_pvals['to cell phenotype']]
interaction_pvals['interaction'] = ['Hom' if f == Epi else 'Het' for f, Epi in zip(interaction_pvals['from'], interaction_pvals['isEpi'])]
interaction_pvals = interaction_pvals.replace(cell_dict)
interaction_pvals['feature_name_unique'] = interaction_pvals['to cell phenotype'] + '__' + interaction_pvals['from'] + interaction_pvals['interaction']
interaction_pvals = interaction_pvals[['Time point', 'feature_name_unique', 'p.value', 'Arm']]

features = pd.concat([density_pvals, ki67_pvals, interaction_pvals])
features = features[features['p.value'] <= 0.05]
sig_features = features.replace({'Cancer_1__TMEHet': 'Cancer_Immune_mixing_score',
                                 'Cancer_2__TMEHet': 'Cancer_Immune_mixing_score',
                                 'Cancer_3__TMEHet': 'Cancer_Immune_mixing_score',
                                 'CD4T__TMEHom': 'Structural_T_mixing_score',
                                 'CD8T__TMEHom': 'Structural_T_mixing_score',
                                 'Treg__TMEHom': 'Structural_T_mixing_score'})
Wang_pre_treatment_features = sig_features[sig_features['Time point'] == 'Baseline']

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

TONIC_feats = list(set(tonic_pre_treatment_features.feature_name_unique.unique()))
TONIC_feats.sort()
WANG_feats = list(set(Wang_pre_treatment_features.feature_name_unique.unique()))
WANG_feats.sort()
NT_SPACECAT_feats = list(set(NT_pre_treatment_features.feature_name_unique.unique()))
NT_SPACECAT_feats.sort()

pre_treatment_table = pd.DataFrame({'TONIC features': pd.Series(TONIC_feats), 'Wang et al features': pd.Series(WANG_feats),
                                    'NT SpaceCat features': pd.Series(NT_SPACECAT_feats)})
pre_treatment_table.to_csv(os.path.join(save_dir, 'Supplementary_Table_8.csv'), index=False)

# feature ranking
feature_rank = pd.read_csv(os.path.join(ANALYSIS_DIR, 'feature_ranking.csv'))
sub_columns = ['feature_name_unique', 'comparison', 'pval', 'fdr_pval', 'med_diff', 'pval_rank', 'cor_rank',
               'combined_rank', 'importance_score', 'signed_importance_score',
               'feature_name', 'compartment', 'cell_pop_level', 'feature_type', 'feature_type_broad']
feature_rank_sub = feature_rank[sub_columns]
feature_rank_sub = feature_rank_sub[feature_rank_sub.comparison.isin(['primary', 'baseline', 'pre_nivo', 'on_nivo'])]

feature_values = pd.read_csv(os.path.join(ANALYSIS_DIR, 'timepoint_combined_features_outcome_labels.csv'))
feature_values = feature_values[feature_values.Timepoint.isin(['primary', 'baseline', 'pre_nivo', 'on_nivo'])]
feature_tp_df = feature_values[['feature_name_unique', 'Timepoint']].drop_duplicates()

feature_tp_df['AUC'] = np.nan
for _, row in feature_tp_df.iloc[:, :2].iterrows():
    feature, tp = row
    data = feature_values.loc[
        np.logical_and(feature_values.feature_name_unique == feature, feature_values.Timepoint == tp)]
    data.reset_index(inplace=True)
    X = data.raw_mean.values.reshape(-1, 1)
    y = data.Clinical_benefit

    # Perform stratified k-fold cross-validation
    kfold_avg_aucs = []
    for i in range(10):
        # Define the number of folds
        k = 3
        skf = StratifiedKFold(n_splits=k, shuffle=True, random_state=i)
        model = LogisticRegression()
        aucs = []
        for train_index, test_index in skf.split(X, y):
            try:
                X_train, X_test = X[train_index], X[test_index]
                y_train, y_test = y[train_index], y[test_index]

                model.fit(X_train, y_train)
                y_pred = model.predict(X_test)
                probabilities = model.predict_proba(X_test)[::, 1]
                auc = metrics.roc_auc_score(y_test, probabilities)
                aucs.append(auc)
            except ValueError:
                continue
        kfold_avg_aucs.append(np.array(aucs).mean())
    mean_auc = np.array(kfold_avg_aucs).mean()
    feature_tp_df.loc[np.logical_and(feature_values.feature_name_unique == feature, feature_values.Timepoint == tp), 'AUC'] = mean_auc

feature_tp_df = feature_tp_df.rename(columns={'Timepoint': 'comparison'})
feature_rank__auc = feature_rank_sub.merge(feature_tp_df, on=['feature_name_unique', 'comparison'], how='left')
feature_rank__auc.to_csv(os.path.join(save_dir, 'Supplementary_Table_9.csv'), index=False)

# top features for multivariate modeling
all_model_rankings = pd.read_csv(os.path.join(BASE_DIR, 'multivariate_lasso/intermediate_results/all_model_rankings.csv'))
top_model_features = all_model_rankings[all_model_rankings.top_ranked]
top_model_features = top_model_features[['timepoint', 'modality', 'feature_name_unique', 'importance_score', 'coef_norm']]
top_model_features.to_csv(os.path.join(save_dir, 'Supplementary_Table_10.csv'), index=False)

# sequencing features
sequencing_features = pd.read_csv(os.path.join(BASE_DIR, 'sequencing_data/processed_genomics_features.csv'))
sequencing_features = sequencing_features[['feature_name', 'data_type', 'feature_type']].drop_duplicates()
sequencing_features = sequencing_features.loc[sequencing_features.feature_type != 'gene_rna', :]

sequencing_features.to_csv(os.path.join(save_dir, 'Supplementary_Table_11.csv'), index=False)
