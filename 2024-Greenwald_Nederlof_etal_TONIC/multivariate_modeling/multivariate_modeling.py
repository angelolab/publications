import os

import pandas as pd
import matplotlib

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy.stats import ttest_ind


base_dir = '/Volumes/Shared/Noah Greenwald/TONIC_Cohort'
plot_dir = os.path.join(base_dir, 'figures')
multivariate_dir = os.path.join(base_dir, 'multivariate_lasso')

ranked_features_univariate = pd.read_csv(os.path.join(base_dir, 'analysis_files/feature_ranking.csv'))
ranked_features_univariate_genomic = pd.read_csv(os.path.join(base_dir, 'sequencing_data/genomics_outcome_ranking.csv'))


# organize scores from cross validation
cv_scores_mibi = pd.read_csv(os.path.join(multivariate_dir, 'all_timepoints_results_MIBI.csv'))
cv_scores_mibi['fold'] = len(cv_scores_mibi)
cv_scores_mibi.columns = ['on_nivo', 'baseline', 'pre_nivo', 'primary', 'fold']

cv_scores_long = pd.melt(cv_scores_mibi, id_vars=['fold'], value_vars=cv_scores_mibi.columns)
cv_scores_long['assay'] = 'MIBI'

cv_scores_genomic = pd.read_csv(os.path.join(multivariate_dir, 'all_timepoints_results_genomics.csv'))
cv_scores_genomic['fold'] = len(cv_scores_genomic)
cv_scores_genomic.columns = ['baseline', 'pre_nivo', 'on_nivo', 'dna baseline', 'fold']

cv_scores_genomic_long = pd.melt(cv_scores_genomic, id_vars=['fold'], value_vars=cv_scores_genomic.columns)
cv_scores_genomic_long['assay'] = 'RNA'
cv_scores_genomic_long.loc[cv_scores_genomic_long.variable == 'dna baseline', 'assay'] = 'DNA'
cv_scores_genomic_long.loc[cv_scores_genomic_long.variable == 'dna baseline', 'variable'] = 'baseline'

cv_scores_long = pd.concat([cv_scores_long, cv_scores_genomic_long])

cv_scores_long.to_csv(os.path.join(base_dir, 'multivariate_lasso', 'formatted_cv_scores.csv'), index=False)

intermediate_dir = os.path.join(base_dir, 'multivariate_lasso', 'intermediate_results')
if not os.path.exists(intermediate_dir):
    os.makedirs(intermediate_dir)


# model diagnostics per timepoint
for timepoint in ['baseline', 'on_nivo', 'pre_nivo', 'primary']:
    for modality in ['MIBI', 'RNA_signature']:
        file_path = os.path.join(base_dir, 'multivariate_lasso', 'top_features_results_{}_{}.csv'.format(timepoint, modality))
        if not os.path.exists(file_path):
            continue

        # read in top features
        top_features = pd.read_csv(file_path)

        # take columns and append into single df
        top_features_long = pd.DataFrame()

        for i in range(10):
            current_fold = top_features.iloc[:, (i * 2):((i * 2) + 2)].reset_index(drop=True)
            current_fold.columns = ['feature', 'coef']
            current_fold = current_fold.loc[~current_fold.feature.isna(), :]
            current_fold['fold'] = i
            current_fold['rank'] = current_fold.index + 1
            current_fold['rank_norm'] = current_fold['rank'] / current_fold.shape[0]
            current_fold['coef_norm'] = current_fold['coef'].abs() / current_fold['coef'].abs().max()
            top_features_long = pd.concat([top_features_long, current_fold])

        # change remove periods from R script modifications
        top_features_long['feature'] = top_features_long['feature'].apply(lambda x: x.replace('.', '+'))

        # look at top 10 features
        top_10_features = top_features_long.loc[top_features_long['rank'] <= 10, :]
        top_10_grouped = top_10_features.groupby('feature').size().sort_values(ascending=False).head(20)

        # plot barplot of top 10 features grouped
        fig, ax = plt.subplots(1, 1, figsize=(6, 6))
        top_10_grouped.plot(kind='bar', ax=ax)
        ax.set_title('Top 10 features')

        plt.tight_layout()
        plt.savefig(os.path.join(intermediate_dir, 'Figure6_top_10_features_{}_{}.pdf'.format(timepoint, modality)))
        plt.close()

        # plot number of occurrences of each feature
        feature_counts = top_features_long.groupby('feature').size().sort_values(ascending=False)
        feature_counts = feature_counts.head(20)
        feature_counts = feature_counts.reset_index()
        feature_counts.columns = ['feature', 'count']

        fig, ax = plt.subplots(1, 1, figsize=(6, 6))
        feature_counts.plot(kind='bar', x='feature', y='count', ax=ax)
        ax.set_title('Feature occurrences')

        plt.tight_layout()
        plt.savefig(os.path.join(intermediate_dir, 'Figure6_feature_occurrences_{}_{}.pdf'.format(timepoint, modality)))
        plt.close()

        # look at average rank
        ranked_features = top_features_long.groupby('feature').agg({'rank_norm': 'mean', 'coef_norm': 'mean', 'rank': 'mean'}).sort_values('rank_norm', ascending=True)
        ranked_features = ranked_features.reset_index()

        ranked_features = ranked_features.merge(feature_counts, on='feature')
        ranked_features.rename(columns={'feature': 'feature_name_unique'}, inplace=True)

        # compare with univariate
        if modality == 'MIBI':
            ranked_features_univ_sub = ranked_features_univariate.loc[ranked_features_univariate.feature_name_unique.isin(ranked_features.feature_name_unique), :]
            ranked_features_univ_sub = ranked_features_univ_sub.loc[ranked_features_univ_sub.comparison == timepoint, :]

            # merge
            combined_rankings = ranked_features_univ_sub[['feature_name_unique', 'feature_rank_global', 'importance_score']].merge(ranked_features, on='feature_name_unique')
            combined_rankings = combined_rankings.sort_values('rank_norm')

        else:
            ranked_features_univ_sub_genomic = ranked_features_univariate_genomic.loc[ranked_features_univariate_genomic.feature_name_unique.isin(ranked_features.feature_name_unique), :]
            ranked_features_univ_sub_genomic = ranked_features_univ_sub_genomic.loc[ranked_features_univ_sub_genomic.comparison == timepoint, :]

            # merge
            combined_rankings = ranked_features_univ_sub_genomic[['feature_name_unique', 'feature_rank_global', 'feature_rank_comparison', 'importance_score']].merge(ranked_features, on='feature_name_unique')
            combined_rankings = combined_rankings.sort_values('rank_norm')

        # save csv
        combined_rankings['top_ranked'] = np.logical_and(combined_rankings.coef_norm.values >= 0.3,
                                                         combined_rankings['count'] >= 3)
        combined_rankings.to_csv(os.path.join(intermediate_dir, 'combined_feature_rankings_{}_{}.csv'.format(timepoint, modality)), index=False)

        # plot
        sns.scatterplot(data=combined_rankings, x='feature_rank_global', y='rank_norm')
        plt.title('Univariate vs multivariate ranking for all nonzero model features')

        plt.savefig(os.path.join(intermediate_dir, 'Figure6_univariate_vs_multivariate_ranking_{}_{}.pdf'.format(timepoint, modality)))
        plt.close()

        #plot feature score broken down by ranking
        fig, ax = plt.subplots(1, 1, figsize=(3, 4))
        sns.stripplot(data=combined_rankings, x='top_ranked', y='importance_score',
                      color='black', ax=ax)
        sns.boxplot(data=combined_rankings, x='top_ranked', y='importance_score',
                    color='grey', ax=ax, showfliers=False)

        ax.set_title('Feature importance score')
        ax.set_ylim([0, 1.1])
        ax.set_xticklabels(['Top ranked', 'Other'])

        sns.despine()
        plt.tight_layout()
        plt.savefig(os.path.join(intermediate_dir, 'Figure6_feature_importance_by_model_{}_{}.pdf'.format(timepoint, modality)))
        plt.close()

# aggregate combined rankings from each timepoint and modality
all_model_rankings = pd.DataFrame()

for timepoint in ['baseline', 'on_nivo', 'pre_nivo', 'primary']:
    for modality in ['MIBI', 'RNA_signature']:
        file_path = os.path.join(intermediate_dir, 'combined_feature_rankings_{}_{}.csv'.format(timepoint, modality))
        if not os.path.exists(file_path):
            continue
        combined_rankings = pd.read_csv(file_path)
        combined_rankings['timepoint'] = timepoint
        combined_rankings['modality'] = modality.replace('_signature', '')
        combined_rankings = combined_rankings[['timepoint', 'modality', 'feature_name_unique', 'feature_rank_global', 'rank', 'rank_norm', 'importance_score', 'count', 'coef_norm', 'top_ranked']]
        all_model_rankings = pd.concat([all_model_rankings, combined_rankings])

# save
all_model_rankings.to_csv(os.path.join(intermediate_dir, 'all_model_rankings.csv'), index=False)
