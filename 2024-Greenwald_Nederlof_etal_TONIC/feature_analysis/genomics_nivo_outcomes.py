import os

import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

import numpy as np
import pandas as pd


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
