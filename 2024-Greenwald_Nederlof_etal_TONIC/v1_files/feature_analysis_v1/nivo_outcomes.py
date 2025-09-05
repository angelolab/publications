import os

import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


from python_files.utils import find_conserved_features, compare_timepoints, compare_populations
from python_files.utils import summarize_population_enrichment, summarize_timepoint_enrichment, compute_feature_enrichment

from statsmodels.stats.multitest import multipletests


plot_dir = '/Users/noahgreenwald/Documents/Grad_School/Lab/TNBC/plots/'
base_dir = '/Volumes/Shared/Noah Greenwald/TONIC_Cohort/'


harmonized_metadata = pd.read_csv(os.path.join(base_dir, 'intermediate_files/metadata/harmonized_metadata.csv'))
patient_metadata = pd.read_csv(os.path.join(base_dir, 'intermediate_files/metadata/TONIC_data_per_patient.csv'))
feature_metadata = pd.read_csv(os.path.join(base_dir, 'analysis_files/feature_metadata.csv'))

#
# To generate the feature rankings, you must have downloaded the patient outcome data.
#
outcome_data = pd.read_csv(os.path.join(base_dir, 'intermediate_files/metadata/patient_clinical_data.csv'))

# load previously computed results
combined_df = pd.read_csv(os.path.join(base_dir, 'analysis_files/timepoint_combined_features.csv'))
combined_df = combined_df.merge(outcome_data, on='Patient_ID')
combined_df = combined_df.loc[combined_df.Clinical_benefit.isin(['Yes', 'No']), :]
combined_df.to_csv(os.path.join(base_dir, 'analysis_files/timepoint_combined_features_outcome_labels.csv'), index=False)

# generate a single set of top hits across all comparisons

# settings for generating hits
plot_hits = False
method = 'ttest'

# placeholder for all values
total_dfs = []

for comparison in combined_df.Timepoint.unique():
    population_df = compare_populations(feature_df=combined_df, pop_col='Clinical_benefit',
                                        timepoints=[comparison], pop_1='No', pop_2='Yes', method=method)

    if plot_hits:
        current_plot_dir = os.path.join(plot_dir, 'responders_nonresponders_{}'.format(comparison))
        if not os.path.exists(current_plot_dir):
            os.makedirs(current_plot_dir)
        summarize_population_enrichment(input_df=population_df, feature_df=combined_df, timepoints=[comparison],
                                        pop_col='Clinical_benefit', output_dir=current_plot_dir, sort_by='med_diff')

    if np.sum(~population_df.log_pval.isna()) == 0:
        continue
    long_df = population_df[['feature_name_unique', 'log_pval', 'mean_diff', 'med_diff']]
    long_df['comparison'] = comparison
    long_df = long_df.dropna()
    long_df['pval'] = 10 ** (-long_df.log_pval)
    long_df['fdr_pval'] = multipletests(long_df.pval, method='fdr_bh')[1]
    total_dfs.append(long_df)


# total_dfs_continuous = []
# for comparison in ['baseline', 'pre_nivo', 'on_nivo', 'baseline__pre_nivo', 'baseline__on_nivo', 'pre_nivo__on_nivo']:
#     input_df = combined_df[combined_df.Timepoint == comparison]
#     continuous_df = compare_continuous(feature_df=input_df, variable_col='Time_to_progression_weeks_RECIST1.1')
#
#     if plot_hits:
#         current_plot_dir = os.path.join(plot_dir, 'responders_nonresponders_continuous_{}'.format(comparison))
#         if not os.path.exists(current_plot_dir):
#             os.makedirs(current_plot_dir)
#         summarize_continuous_enrichment(input_df=continuous_df, feature_df=combined_df, timepoint=comparison,
#                                         variable_col='Time_to_progression_weeks_RECIST1.1', output_dir=current_plot_dir, min_score=0.95)

# summarize hits from all comparisons
ranked_features_df = pd.concat(total_dfs)
ranked_features_df['log10_qval'] = -np.log10(ranked_features_df.fdr_pval)

# get ranking of each row by FDR p-value and correlation
ranked_features_df['pval_rank'] = ranked_features_df.fdr_pval.rank(ascending=True)
ranked_features_df['cor_rank'] = ranked_features_df.med_diff.abs().rank(ascending=False)
ranked_features_df['combined_rank'] = (ranked_features_df.pval_rank.values + ranked_features_df.cor_rank.values) / 2

# generate importance score
max_rank = len(~ranked_features_df.med_diff.isna())
normalized_rank = ranked_features_df.combined_rank / max_rank
ranked_features_df['importance_score'] = 1 - normalized_rank

ranked_features_df = ranked_features_df.sort_values('importance_score', ascending=False)

# generate signed version of score
ranked_features_df['signed_importance_score'] = ranked_features_df.importance_score * np.sign(ranked_features_df.med_diff)

# add feature type
ranked_features_df = ranked_features_df.merge(feature_metadata, on='feature_name_unique', how='left')

feature_type_dict = {'functional_marker': 'phenotype', 'linear_distance': 'interactions',
                     'density': 'density', 'cell_diversity': 'diversity', 'density_ratio': 'density',
                     'mixing_score': 'interactions', 'region_diversity': 'diversity',
                     'compartment_area_ratio': 'compartment', 'density_proportion': 'density',
                      'morphology': 'phenotype', 'pixie_ecm': 'ecm', 'fiber': 'ecm', 'ecm_cluster': 'ecm',
                        'compartment_area': 'compartment', 'ecm_fraction': 'ecm'}
ranked_features_df['feature_type_broad'] = ranked_features_df.feature_type.map(feature_type_dict)

# get ranking of each feature
ranked_features_df['feature_rank_global_evolution'] = ranked_features_df.importance_score.rank(ascending=False)

# get ranking of non-evolution features
ranked_features_no_evo = ranked_features_df.loc[ranked_features_df.comparison.isin(['primary', 'baseline', 'pre_nivo', 'on_nivo']), :]
ranked_features_no_evo['feature_rank_global'] = ranked_features_no_evo.importance_score.rank(ascending=False)
ranked_features_df = ranked_features_df.merge(ranked_features_no_evo.loc[:, ['feature_name_unique', 'comparison', 'feature_rank_global']], on=['feature_name_unique', 'comparison'], how='left')

# get ranking for each comparison
ranked_features_df['feature_rank_comparison'] = np.nan
for comparison in ranked_features_df.comparison.unique():
    # get subset of features from given comparison
    ranked_features_comp = ranked_features_df.loc[ranked_features_df.comparison == comparison, :]
    ranked_features_comp['temp_comparison'] = ranked_features_comp.importance_score.rank(ascending=False)

    # merge with placeholder column
    ranked_features_df = ranked_features_df.merge(ranked_features_comp.loc[:, ['feature_name_unique', 'comparison', 'temp_comparison']], on=['feature_name_unique', 'comparison'], how='left')

    # replace with values from placeholder, then delete
    ranked_features_df['feature_rank_comparison'] = ranked_features_df['temp_comparison'].fillna(ranked_features_df['feature_rank_comparison'])
    ranked_features_df.drop(columns='temp_comparison', inplace=True)

# saved formatted df
ranked_features_df.to_csv(os.path.join(base_dir, 'analysis_files/feature_ranking.csv'), index=False)


# same thing for genomics features
sequence_dir = os.path.join(base_dir, 'sequencing_data')
genomics_df = pd.read_csv(os.path.join(sequence_dir, 'processed_genomics_features.csv'))

genomics_df = pd.merge(genomics_df, outcome_data, on='Patient_ID')

plot_hits = False
method = 'ttest'

genomics_df = genomics_df.loc[genomics_df.Timepoint != 'on_nivo_1_cycle', :]
genomics_df = genomics_df.rename(columns={'feature_name': 'feature_name_unique'})
genomics_df = genomics_df.loc[genomics_df.feature_type != 'gene_rna', :]

# placeholder for all values
total_dfs = []

for comparison in genomics_df.Timepoint.unique():
    population_df = compare_populations(feature_df=genomics_df, pop_col='Clinical_benefit',
                                        timepoints=[comparison], pop_1='No', pop_2='Yes', method=method,
                                        feature_suff='value')

    if plot_hits:
        current_plot_dir = os.path.join(plot_dir, 'responders_nonresponders_{}'.format(comparison))
        if not os.path.exists(current_plot_dir):
            os.makedirs(current_plot_dir)
        summarize_population_enrichment(input_df=population_df, feature_df=genomics_df, timepoints=[comparison],
                                        pop_col='Clinical_benefit', output_dir=current_plot_dir, sort_by='med_diff')

    if np.sum(~population_df.log_pval.isna()) == 0:
        continue
    long_df = population_df[['feature_name_unique', 'log_pval', 'mean_diff', 'med_diff']]
    long_df['comparison'] = comparison
    long_df = long_df.dropna()
    long_df['pval'] = 10 ** (-long_df.log_pval)
    long_df['fdr_pval'] = multipletests(long_df.pval, method='fdr_bh')[1]
    total_dfs.append(long_df)


ranked_genomics_df = pd.concat(total_dfs)
ranked_genomics_df['log10_qval'] = -np.log10(ranked_genomics_df.fdr_pval)

# get ranking of each row by pval and correlation
ranked_genomics_df['pval_rank'] = ranked_genomics_df.log_pval.rank(ascending=False)
ranked_genomics_df['cor_rank'] = ranked_genomics_df.med_diff.abs().rank(ascending=False)
ranked_genomics_df['combined_rank'] = (ranked_genomics_df.pval_rank.values + ranked_genomics_df.cor_rank.values) / 2

# generate importance score
max_rank = len(~ranked_genomics_df.med_diff.isna())
normalized_rank = ranked_genomics_df.combined_rank / max_rank
ranked_genomics_df['importance_score'] = 1 - normalized_rank

ranked_genomics_df = ranked_genomics_df.sort_values('importance_score', ascending=False)

# generate signed version of score
ranked_genomics_df['signed_importance_score'] = ranked_genomics_df.importance_score * np.sign(ranked_genomics_df.med_diff)

# get ranking of each feature
ranked_genomics_df['feature_rank_global'] = ranked_genomics_df.importance_score.rank(ascending=False)

# get ranking for each comparison
ranked_genomics_df['feature_rank_comparison'] = np.nan
for comparison in ranked_genomics_df.comparison.unique():
    # get subset of features from given comparison
    ranked_features_comp = ranked_genomics_df.loc[ranked_genomics_df.comparison == comparison, :]
    ranked_features_comp['temp_comparison'] = ranked_features_comp.importance_score.rank(ascending=False)

    # merge with placeholder column
    ranked_genomics_df = ranked_genomics_df.merge(ranked_features_comp.loc[:, ['feature_name_unique', 'comparison', 'temp_comparison']], on=['feature_name_unique', 'comparison'], how='left')

    # replace with values from placeholder, then delete
    ranked_genomics_df['feature_rank_comparison'] = ranked_genomics_df['temp_comparison'].fillna(ranked_genomics_df['feature_rank_comparison'])
    ranked_genomics_df.drop(columns='temp_comparison', inplace=True)

# saved formatted df
genomics_df = genomics_df.rename(columns={'feature_name': 'feature_name_unique'})
genomics_df = genomics_df[['feature_name_unique', 'feature_type', 'data_type']].drop_duplicates()

ranked_genomics_df = ranked_genomics_df.merge(genomics_df, on='feature_name_unique', how='left')
ranked_genomics_df.to_csv(os.path.join(sequence_dir, 'genomics_outcome_ranking.csv'), index=False)

if plot_hits:
    # plot top X features per comparison
    num_features = 30

    # ranked_features_df = ranked_genomics_df
    # combined_df = genomics_df

    for comparison in ranked_features_df.comparison.unique():
        current_plot_dir = os.path.join(plot_dir, 'top_features_{}'.format(comparison))
        if not os.path.exists(current_plot_dir):
            os.makedirs(current_plot_dir)

        current_df = ranked_features_df.loc[ranked_features_df.comparison == comparison, :]
        current_df = current_df.sort_values('combined_rank', ascending=True)
        current_df = current_df.iloc[:num_features, :]

        # plot results
        for feature_name, rank in zip(current_df.feature_name_unique.values, current_df.combined_rank.values):
            plot_df = combined_df.loc[(combined_df.feature_name_unique == feature_name) &
                                      (combined_df.Timepoint == comparison), :]

            g = sns.catplot(data=plot_df, x='Clinical_benefit', y='raw_value', kind='strip')
            g.fig.suptitle(feature_name)
            g.savefig(os.path.join(current_plot_dir, 'rank_{}_feature_{}.png'.format(rank, feature_name)))
            plt.close()


    # plot top X features overall
    num_features = 30
    output_dir = os.path.join(plot_dir, 'top_features_overall')
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    current_df = ranked_features_df.sort_values('combined_rank', ascending=True)
    current_df = current_df.loc[current_df.comparison.isin(['baseline', 'pre_nivo', 'on_nivo', 'primary']), :]
    current_df = current_df.iloc[:num_features, :]

    # plot results
    for feature_name, comparison, rank in zip(current_df.feature_name_unique.values, current_df.comparison.values, current_df.combined_rank.values):
        plot_df = combined_df.loc[(combined_df.feature_name_unique == feature_name) &
                                  (combined_df.Timepoint == comparison), :]

        g = sns.catplot(data=plot_df, x='Clinical_benefit', y='raw_mean', kind='strip')
        g.fig.suptitle(feature_name)
        g.savefig(os.path.join(output_dir, 'rank_{}_feature_{}.png'.format(rank, feature_name)))
        plt.close()





# compare subsets of features to see effect of dox only

# # create different input dfs
# default_df = timepoint_features.copy()
# control_df = timepoint_features.loc[timepoint_features.Induction_treatment != 'No induction', :]
# dox_df = timepoint_features.loc[timepoint_features.Induction_treatment == 'Doxorubicin', :]
#
# len(dox_df.loc[(dox_df.Timepoint == 'pre_nivo'), 'Patient_ID'].unique())
# induction_pats = default_df.loc[(default_df.Timepoint == 'pre_nivo'), 'Patient_ID'].unique()
# len(induction_pats)
#
# # pick random subset from induction patients
# np.random.seed(42)
# induction_pat_subsets = [np.random.choice(induction_pats, size=29, replace=False) for _ in range(10)]
#
# dfs = [[default_df, 'default'], [control_df, 'control'], [dox_df, 'dox_only']]
# random_dfs = [[default_df.loc[default_df.Patient_ID.isin(vals), :], 'subset_{}'.format(idx)] for idx, vals in enumerate(induction_pat_subsets)]
#
# dfs = dfs + random_dfs
#
# # generate top hits for each dataset
# total_dfs = []
# for input_df, name in dfs:
#     population_df = compare_populations(feature_df=input_df, pop_col='Clinical_benefit',
#                                             timepoints=['pre_nivo'], pop_1='No', pop_2='Yes', method='ttest')
#
#
#     long_df = population_df[['feature_name_unique', 'log_pval', 'mean_diff', 'med_diff']]
#     long_df['comparison'] = name
#     long_df = long_df.dropna()
#
#     long_df['pval_rank'] = long_df.log_pval.rank(ascending=False)
#     long_df['cor_rank'] = long_df.med_diff.abs().rank(ascending=False)
#     long_df['combined_rank'] = (long_df.pval_rank.values + long_df.cor_rank.values) / 2
#
#     # add to overall df
#     total_dfs.append(long_df)
#
# total_dfs = pd.concat(total_dfs)
#
# wide_df = total_dfs.pivot(index='feature_name_unique', columns='comparison', values='combined_rank')
# wide_df = wide_df.reset_index()
# wide_df = wide_df.loc[wide_df.default < 100]
#
# # plot results
# plot_df = wide_df.loc[:, ['default', 'subset_5']]
# plot_df.dropna(axis=0, inplace=True)
# sns.scatterplot(data=plot_df, x='default', y='subset_5')
# plt.title('All samples vs. subset 5')
#
# corr, _ = spearmanr(plot_df.default, plot_df.subset_5)
# plt.text(0.5, 0.5, 'Spearman R = {:.2f}'.format(corr), transform=plt.gca().transAxes)
# plt.savefig(os.path.join(plot_dir, 'Feature_Correlation_subset_5.pdf'))
# plt.close()
#
# corrs = []
# for i in range(10):
#     col = 'subset_{}'.format(i)
#     plot_df = wide_df.loc[:, ['default', col]]
#     plot_df.dropna(axis=0, inplace=True)
#     corr, _ = spearmanr(plot_df.default, plot_df[col])
#     corrs.append(corr)
#
# cor_df = pd.DataFrame(corrs, columns=['correlation'])
#
# plot_df = wide_df.loc[:, ['default', 'dox_only']]
# plot_df.dropna(axis=0, inplace=True)
# dox_cor, _ = spearmanr(plot_df.default, plot_df.dox_only)
#
# cor_df = pd.concat([cor_df, (pd.DataFrame([dox_cor], columns=['correlation']))])
# cor_df['type'] = 'randomized'
# cor_df.iloc[-1, -1] = 'dox_only'
#
# sns.stripplot(data=cor_df, x='type', y='correlation')
# plt.savefig(os.path.join(plot_dir, 'Feature_Correlation_summary.pdf'))
# plt.close()