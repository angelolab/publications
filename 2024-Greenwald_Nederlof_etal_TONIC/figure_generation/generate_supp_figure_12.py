import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from matplotlib_venn import venn3


BASE_DIR = "/Volumes/Shared/Noah Greenwald/TONIC_Cohort/"
SUPPLEMENTARY_FIG_DIR = os.path.join(BASE_DIR, "supplementary_figs")

all_model_rankings = pd.read_csv(os.path.join(BASE_DIR, 'multivariate_lasso/intermediate_results', 'all_model_rankings.csv'))

# plot top features
all_model_plot = all_model_rankings.loc[all_model_rankings.timepoint != 'primary', :]
sns.stripplot(data=all_model_plot.loc[all_model_plot.top_ranked, :], x='timepoint', y='importance_score', hue='modality',
              order=['baseline', 'pre_nivo', 'on_nivo'])
plt.title('Top ranked features')
plt.ylim([0, 1.05])
plt.savefig(os.path.join(SUPPLEMENTARY_FIG_DIR, 'supp_figure_12d.pdf'))
plt.close()

# plot number of times features are selected
sns.histplot(data=all_model_rankings.loc[all_model_rankings.top_ranked, :], x='count', color='grey', multiple='stack',
             binrange=(1, 10), discrete=True)
plt.title('Number of times features are selected')
plt.savefig(os.path.join(SUPPLEMENTARY_FIG_DIR, 'supp_figure_12b.pdf'))
plt.close()

# plot venn diagram
rna_rankings_top = all_model_rankings.loc[np.logical_and(all_model_rankings.modality == 'RNA', all_model_rankings.top_ranked), :]
rna_baseline = rna_rankings_top.loc[rna_rankings_top.timepoint == 'baseline', 'feature_name_unique'].values
rna_nivo = rna_rankings_top.loc[rna_rankings_top.timepoint == 'on_nivo', 'feature_name_unique'].values
rna_induction = rna_rankings_top.loc[rna_rankings_top.timepoint == 'pre_nivo', 'feature_name_unique'].values

venn3([set(rna_baseline), set(rna_nivo), set(rna_induction)], ('Baseline', 'Nivo', 'Induction'))
plt.title('RNA top ranked features')
plt.savefig(os.path.join(SUPPLEMENTARY_FIG_DIR, 'supp_figure_12f.pdf'))
plt.close()

# top ranked features from each timepoint
mibi_rankings_top = all_model_rankings.loc[np.logical_and(all_model_rankings.modality == 'MIBI', all_model_rankings.top_ranked), :]
mibi_baseline = mibi_rankings_top.loc[mibi_rankings_top.timepoint == 'baseline', 'feature_name_unique'].values
mibi_nivo = mibi_rankings_top.loc[mibi_rankings_top.timepoint == 'on_nivo', 'feature_name_unique'].values
mibi_induction = mibi_rankings_top.loc[mibi_rankings_top.timepoint == 'pre_nivo', 'feature_name_unique'].values

venn3([set(mibi_baseline), set(mibi_nivo), set(mibi_induction)], ('Baseline', 'Nivo', 'Induction'))
plt.title('MIBI top ranked features')
plt.savefig(os.path.join(SUPPLEMENTARY_FIG_DIR, 'supp_figure_12e.pdf'))
plt.close()

# compare correlations between top ranked features
ranked_features_univariate = pd.read_csv(os.path.join(BASE_DIR, 'analysis_files/feature_ranking.csv'))

nivo_features_model = all_model_rankings.loc[np.logical_and(all_model_rankings.timepoint == 'on_nivo', all_model_rankings.top_ranked), :]
nivo_features_model = nivo_features_model.loc[nivo_features_model.modality == 'MIBI', 'feature_name_unique'].values

nivo_features_univariate = ranked_features_univariate.loc[np.logical_and(ranked_features_univariate.comparison == 'on_nivo',
                                                                         ranked_features_univariate.feature_rank_global <= 100), :]

timepoint_features = pd.read_csv(os.path.join(BASE_DIR, 'analysis_files/timepoint_combined_features.csv'))
timepoint_features = timepoint_features.loc[timepoint_features.Timepoint == 'on_nivo', :]

timepoint_features_model = timepoint_features.loc[timepoint_features.feature_name_unique.isin(nivo_features_model), :]
timepoint_features_model = timepoint_features_model[['feature_name_unique', 'normalized_mean', 'Patient_ID']]
timepoint_features_model = timepoint_features_model.pivot(index='Patient_ID', columns='feature_name_unique', values='normalized_mean')

# get values
model_corr = timepoint_features_model.corr()
model_corr = model_corr.where(np.triu(np.ones(model_corr.shape), k=1).astype(np.bool)).values.flatten()
model_corr = model_corr[~np.isnan(model_corr)]

# get univariate features
timepoint_features_univariate = timepoint_features.loc[timepoint_features.feature_name_unique.isin(nivo_features_univariate.feature_name_unique), :]
timepoint_features_univariate = timepoint_features_univariate[['feature_name_unique', 'normalized_mean', 'Patient_ID']]
timepoint_features_univariate = timepoint_features_univariate.pivot(index='Patient_ID', columns='feature_name_unique', values='normalized_mean')

# get values
univariate_corr = timepoint_features_univariate.corr()
univariate_corr = univariate_corr.where(np.triu(np.ones(univariate_corr.shape), k=1).astype(np.bool)).values.flatten()
univariate_corr = univariate_corr[~np.isnan(univariate_corr)]

corr_values = pd.DataFrame({'correlation': np.concatenate([model_corr, univariate_corr]),
                            'model': ['model'] * len(model_corr) + ['univariate'] * len(univariate_corr)})

# plot correlations by model
fig, ax = plt.subplots(1, 1, figsize=(3, 4))
sns.boxplot(data=corr_values, x='model', y='correlation',
            color='grey', ax=ax, showfliers=False)

ax.set_title('Feature correlation')
sns.despine()
plt.tight_layout()
plt.savefig(os.path.join(SUPPLEMENTARY_FIG_DIR, 'supp_figure_12c.pdf'))
plt.close()

# look at interesting features
combined_df = pd.read_csv(os.path.join(BASE_DIR, 'analysis_files/timepoint_combined_features.csv'))

for timepoint in ['primary', 'baseline', 'pre_nivo', 'on_nivo']:

    plot_df = combined_df.loc[(combined_df.feature_name_unique == 'CD38+__all') &
                              (combined_df.Timepoint == timepoint), :]

    fig, ax = plt.subplots(1, 1, figsize=(2, 4))
    sns.stripplot(data=plot_df, x='Clinical_benefit', y='raw_mean', order=['Yes', 'No'],
                    color='black', ax=ax)
    sns.boxplot(data=plot_df, x='Clinical_benefit', y='raw_mean', order=['Yes', 'No'],
                    color='grey', ax=ax, showfliers=False, width=0.3)
    ax.set_title('CD38+ ' + timepoint)
    ax.set_ylim([0, 0.5])
    sns.despine()
    plt.tight_layout()
    plt.savefig(os.path.join(SUPPLEMENTARY_FIG_DIR, 'CD38_positivity in {}.pdf'.format(timepoint)))
    plt.close()
