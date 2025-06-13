import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from matplotlib_venn import venn3
from alpineer.io_utils import list_files


BASE_DIR = "/Volumes/Shared/Noah Greenwald/TONIC_Cohort/"
SUPPLEMENTARY_FIG_DIR = os.path.join(BASE_DIR, "supplementary_figs")
REVIEW_FIG_DIR = os.path.join(BASE_DIR, "supplementary_figs/review_figures")

all_model_rankings = pd.read_csv(os.path.join(BASE_DIR, 'multivariate_lasso/intermediate_results', 'all_model_rankings.csv'))

# plot top features
all_model_plot = all_model_rankings.loc[all_model_rankings.timepoint != 'primary', :]
sns.stripplot(data=all_model_plot.loc[all_model_plot.top_ranked, :], x='timepoint', y='importance_score', hue='modality',
              order=['baseline', 'pre_nivo', 'on_nivo'])
plt.title('Top ranked features')
plt.ylim([0, 1.05])
plt.savefig(os.path.join(SUPPLEMENTARY_FIG_DIR, 'supp_figure_15c.pdf'))
plt.close()

# plot number of times features are selected
sns.histplot(data=all_model_rankings.loc[all_model_rankings.top_ranked, :], x='count', color='grey', multiple='stack',
             binrange=(1, 10), discrete=True)
plt.title('Number of times features are selected')
plt.savefig(os.path.join(SUPPLEMENTARY_FIG_DIR, 'supp_figure_15b.pdf'))
plt.close()

# plot venn diagram
rna_rankings_top = all_model_rankings.loc[np.logical_and(all_model_rankings.modality == 'RNA', all_model_rankings.top_ranked), :]
rna_baseline = rna_rankings_top.loc[rna_rankings_top.timepoint == 'baseline', 'feature_name_unique'].values
rna_nivo = rna_rankings_top.loc[rna_rankings_top.timepoint == 'on_nivo', 'feature_name_unique'].values
rna_induction = rna_rankings_top.loc[rna_rankings_top.timepoint == 'pre_nivo', 'feature_name_unique'].values

venn3([set(rna_baseline), set(rna_nivo), set(rna_induction)], ('Baseline', 'Nivo', 'Induction'))
plt.title('RNA top ranked features')
plt.savefig(os.path.join(SUPPLEMENTARY_FIG_DIR, 'supp_figure_15e.pdf'))
plt.close()

# top ranked features from each timepoint
mibi_rankings_top = all_model_rankings.loc[np.logical_and(all_model_rankings.modality == 'MIBI', all_model_rankings.top_ranked), :]
mibi_baseline = mibi_rankings_top.loc[mibi_rankings_top.timepoint == 'baseline', 'feature_name_unique'].values
mibi_nivo = mibi_rankings_top.loc[mibi_rankings_top.timepoint == 'on_nivo', 'feature_name_unique'].values
mibi_induction = mibi_rankings_top.loc[mibi_rankings_top.timepoint == 'pre_nivo', 'feature_name_unique'].values

venn3([set(mibi_baseline), set(mibi_nivo), set(mibi_induction)], ('Baseline', 'Nivo', 'Induction'))
plt.title('MIBI top ranked features')
plt.savefig(os.path.join(SUPPLEMENTARY_FIG_DIR, 'supp_figure_15d.pdf'))
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

'''
# plot correlations by model
fig, ax = plt.subplots(1, 1, figsize=(3, 4))
sns.boxplot(data=corr_values, x='model', y='correlation',
            color='grey', ax=ax, showfliers=False)

ax.set_title('Feature correlation')
sns.despine()
plt.tight_layout()
plt.savefig(os.path.join(SUPPLEMENTARY_FIG_DIR, 'supp_figure_15c.pdf'))
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
'''

## 4.4 Limit multivariate model features ##
limit_features_dir = os.path.join(REVIEW_FIG_DIR, 'limit_model_features')
prediction_dir = os.path.join(BASE_DIR, 'prediction_model')

'''
for feature_cap_num in [5, 10, 15, 20]:
    feature_cap_dir = os.path.join(limit_features_dir, f'prediction_model_feature_cap_num')
    os.makedirs(os.path.join(feature_cap_dir, 'patient_outcomes'), exist_ok=True)
    
## RUN all_timepoints-feature_cap.R for each feature_cap_num
'''

preds_5 = pd.read_csv(os.path.join(limit_features_dir, 'prediction_model_5/patient_outcomes/all_timepoints_results_MIBI.csv'))
preds_10 = pd.read_csv(os.path.join(limit_features_dir, 'prediction_model_10/patient_outcomes/all_timepoints_results_MIBI.csv'))
preds_15 = pd.read_csv(os.path.join(limit_features_dir, 'prediction_model_15/patient_outcomes/all_timepoints_results_MIBI.csv'))
preds_20 = pd.read_csv(os.path.join(limit_features_dir, 'prediction_model_20/patient_outcomes/all_timepoints_results_MIBI.csv'))
preds_reg = pd.read_csv(os.path.join(prediction_dir, 'patient_outcomes/all_timepoints_results_MIBI.csv'))

df_5 = preds_5.mean()
df_10 = preds_10.mean()
df_15 = preds_15.mean()
df_20 = preds_20.mean()
df_reg = preds_reg.mean()
df = pd.concat([df_reg, df_20, df_15, df_10, df_5], axis=1)
df = df.rename(columns={0: 'All features', 1: 'Top 20 features', 2: 'Top 15 features', 3: 'Top 10 features', 4: 'Top 5 features'})

df = df.reset_index()
df.replace('auc_on_nivo_list', 'On nivo', inplace=True)
df.replace('auc_post_induction_list', 'Pre nivo', inplace=True)
df.replace('auc_primary_list', 'Primary', inplace=True)
df.replace('auc_baseline_list', 'Baseline', inplace=True)
df['order'] = df['index'].replace({'Primary':0, 'Baseline':1, 'Pre nivo':2, 'On nivo': 3})
df = df.sort_values(by='order')
df = df.drop(columns=['order'])
df = df.rename(columns={'index': 'Timepoint'})
df = pd.melt(df, ['Timepoint'])

sns.scatterplot(data=df[df.variable == 'All features'], x='Timepoint', y='value', hue='variable', palette=sns.color_palette(['black']), edgecolors='black')
sns.scatterplot(data=df[df.variable != 'All features'], x='Timepoint', y='value', hue='variable', palette=sns.color_palette(['dimgrey', 'darkgrey', 'lightgrey', 'whitesmoke']), edgecolors='black')
plt.xticks(rotation=30)
plt.ylabel('Mean AUC')
plt.title('Model accuracy by feature amount')
sns.despine()
plt.gca().legend(loc='lower right').set_title('')
plt.savefig(os.path.join(SUPPLEMENTARY_FIG_DIR, 'supp_figure_15h.pdf'), bbox_inches='tight', dpi=300)

# multimodal prediction plots
multi_dir = os.path.join(BASE_DIR, 'sequencing_data/multimodal_prediction')
baseline_results = pd.read_csv(os.path.join(multi_dir, 'baseline_results.csv'))
pre_nivo_results = pd.read_csv(os.path.join(multi_dir, 'pre_nivo_results.csv'))
on_nivo_results = pd.read_csv(os.path.join(multi_dir, 'on_nivo_results.csv'))

baseline_results['Timepoint'] = 'Baseline'
baseline_results = baseline_results.rename(columns={'auc_baseline_list': 'Combined', 'auc_rna_baseline_list': 'RNA', 'auc_protein_baseline_list': 'MIBI'})
pre_nivo_results['Timepoint'] = 'Pre nivo'
pre_nivo_results = pre_nivo_results.rename(columns={'auc_induction_list': 'Combined', 'auc_rna_induction_list': 'RNA', 'auc_protein_induction_list': 'MIBI'})
on_nivo_results['Timepoint'] = 'On nivo'
on_nivo_results = on_nivo_results.rename(columns={'auc_on_nivo_list': 'Combined', 'auc_rna_on_nivo_list': 'RNA', 'auc_protein_on_nivo_list': 'MIBI'})
all_results = pd.concat([baseline_results, pre_nivo_results, on_nivo_results])
all_results = all_results[['Timepoint', 'Combined']]

fig, ax = plt.subplots()
sns.boxplot(data=all_results, x='Timepoint', y='Combined', ax=ax, width=0.6,
            palette=sns.color_palette(["chocolate"]), showfliers=False)
sns.stripplot(data=all_results, x='Timepoint', y='Combined', ax=ax,
              palette=sns.color_palette(["chocolate"]), jitter=0.1)
fig.set_figheight(5)
fig.set_figwidth(6)
plt.xticks(rotation=45)
plt.title('Combined MIBI & RNA multivariate model accuracy')
plt.ylabel('AUC')
plt.xlabel('')
plt.ylim((0, 1))
sns.despine()
plt.savefig(os.path.join(SUPPLEMENTARY_FIG_DIR, 'supp_figure_15i.pdf'), bbox_inches='tight', dpi=300)
