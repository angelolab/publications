import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from matplotlib_venn import venn2


from python_files.supplementary_plot_helpers import random_feature_generation, calculate_feature_corr


BASE_DIR = "/Volumes/Shared/Noah Greenwald/TONIC_Cohort/"
SUPPLEMENTARY_FIG_DIR = os.path.join(BASE_DIR, "supplementary_figs")


# Analyse the significance scores of top features after randomization compared to the TONIC data.
fp_dir = os.path.join(SUPPLEMENTARY_FIG_DIR, 'supp_figure_9a_files')
if not os.path.exists(fp_dir):
    os.makedirs(fp_dir)

'''
# compute random feature sets
combined_df = pd.read_csv(os.path.join(BASE_DIR, 'analysis_files/timepoint_combined_features.csv'))
feature_df = pd.read_csv(os.path.join(BASE_DIR, 'analysis_files/feature_ranking.csv'))
feature_metadata = pd.read_csv(os.path.join(BASE_DIR, 'analysis_files/feature_metadata.csv'))

repeated_features, repeated_features_num, scores = [], [], []
overlapping_features, random_top_features = [], []

sample_num = 100
np.random.seed(13)

for i, seed in enumerate(random.sample(range(1, 2000), sample_num)):
    print(f'{i+1}/100')
    intersection_of_features, jaccard_score, top_random_features = random_feature_generation(combined_df, seed, feature_df[:100], feature_metadata)

    shared_df = pd.DataFrame({
        'random_seed': [seed] * len(intersection_of_features),
        'repeated_features' : list(intersection_of_features),
        'jaccard_score': [jaccard_score] * len(intersection_of_features)
    })
    overlapping_features.append(shared_df)

    top_random_features['seed'] = seed
    random_top_features.append(top_random_features)

results = pd.concat(overlapping_features)
top_features = pd.concat(random_top_features)
# add TONIC features to data with seed 0
top_features = pd.concat([top_features, feature_df[:100]])
top_features['seed'] = top_features['seed'].fillna(0)

results.to_csv(os.path.join(fp_dir, 'overlapping_features.csv'), index=False)
top_features.to_csv(os.path.join(fp_dir, 'top_features.csv'), index=False)
'''

top_features = pd.read_csv(os.path.join(fp_dir, 'top_features.csv'))
results = pd.read_csv(os.path.join(fp_dir, 'overlapping_features.csv'))

avg_scores = top_features[['seed', 'pval', 'log_pval', 'fdr_pval', 'med_diff']].groupby(by='seed').mean()
avg_scores['abs_med_diff'] = abs(avg_scores['med_diff'])
top_features['abs_med_diff'] = abs(top_features['med_diff'])

# log p-value & effect size plots
for name, metric in zip(['supp_figure_9a', 'supp_figure_9c'], ['log_pval', 'abs_med_diff']):
    # plot metric dist in top features for TONIC data and one random set
    TONIC = top_features[top_features.seed == 0]
    random = top_features[top_features.seed == 228]
    g = sns.distplot(TONIC[metric], kde=True, color='#1f77b4')
    g = sns.distplot(random[metric], kde=True, color='#ff7f0e')
    g.set(xlim=(0, None))
    plt.xlabel(metric)
    plt.title(f"{metric} Distribution in TONIC vs a Random")
    g.legend(labels=["TONIC", "Randomized"])
    sns.move_legend(g, "upper left", bbox_to_anchor=(0.9, 1))
    sns.despine()
    plt.tight_layout()
    plt.savefig(os.path.join(SUPPLEMENTARY_FIG_DIR, f"{name}.pdf"), dpi=300)
    plt.show()
    plt.close()

for name, metric in zip(['supp_figure_9b', 'supp_figure_9d'], ['log_pval', 'abs_med_diff']):
    # plot average metric across top features for each set
    g = sns.distplot(avg_scores[metric][1:], kde=True,  color='#ff7f0e')
    g.axvline(x=avg_scores[metric][0], color='#1f77b4')
    g.set(xlim=(0, avg_scores[metric][0]*1.2))
    plt.xlabel(f'Average {metric} of Top 100 Features')
    plt.title(f"Average {metric} in TONIC vs Random Sets")
    g.legend(labels=["Randomized", "TONIC"])
    sns.move_legend(g, "upper left", bbox_to_anchor=(0.9, 1))
    sns.despine()
    plt.tight_layout()
    plt.savefig(os.path.join(SUPPLEMENTARY_FIG_DIR, f"{name}.pdf"), dpi=300)
    plt.show()
    plt.close()

# # general feature overlap plots
# high_features = results.groupby(by='repeated_features').count().sort_values(by='random_seed', ascending=False).reset_index()
# high_features = high_features[high_features.random_seed>3].sort_values(by='random_seed')
# plt.barh(high_features.repeated_features, high_features.random_seed)
# plt.xlabel('How Many Random Sets Contain the Feature')
# plt.title('Repeated Top Features')
# sns.despine()
# plt.savefig(os.path.join(fp_dir, "Repeated_Top_Features.pdf"), dpi=300, bbox_inches='tight')
# plt.show()

repeated_features_num = results.groupby(by='random_seed').count().sort_values(by='repeated_features', ascending=False)
plt.hist(repeated_features_num.repeated_features)
plt.xlabel('Number of TONIC Top Features in each Random Set')
plt.title('Histogram of Overlapping Features')
sns.despine()
plt.savefig(os.path.join(SUPPLEMENTARY_FIG_DIR, f"supp_figure_9e.pdf"), dpi=300)
plt.show()
plt.close()

# plt.hist(results.jaccard_score, bins=10)
# plt.xlim((0, 0.10))
# plt.title('Histogram of Jaccard Scores')
# sns.despine()
# plt.xlabel('Jaccard Score')
# plt.savefig(os.path.join(fp_dir, "Histogram_of_Jaccard_Scores.pdf"), dpi=300)
# plt.show()

# compute the correlation between response-associated features
timepoint_features = pd.read_csv(os.path.join(BASE_DIR, 'analysis_files/timepoint_combined_features.csv'))
feature_ranking_df = pd.read_csv(os.path.join(BASE_DIR, 'analysis_files/feature_ranking.csv'))
feature_ranking_df = feature_ranking_df[np.isin(feature_ranking_df['comparison'], ['primary', 'baseline', 'pre_nivo' , 'on_nivo'])]
feature_ranking_df = feature_ranking_df.sort_values(by = 'feature_rank_global', ascending=True)

top_features = feature_ranking_df.iloc[:100, :].loc[:, ['feature_name_unique', 'comparison', 'feature_type']]
top_features.columns = ['feature_name_unique', 'Timepoint', 'feature_type']

remaining_features = feature_ranking_df.iloc[100:, :].loc[:, ['feature_name_unique', 'comparison', 'feature_type']]
remaining_features.columns = ['feature_name_unique', 'Timepoint', 'feature_type']


#C(100, 2) = 100! / [(100-2)! * 2!] = 4950 unique pairwise combinations top 100 features and each other
corr_within = calculate_feature_corr(timepoint_features, top_features, remaining_features, top=True)
corr_across = calculate_feature_corr(timepoint_features, top_features, remaining_features, top=False)

corr_across = corr_across[~np.isnan(corr_across)]
corr_within = corr_within[~np.isnan(corr_within)]

_, axes = plt.subplots(1, 1, figsize = (5, 4), gridspec_kw={'hspace': 0.45, 'wspace': 0.4, 'bottom':0.15})
g = sns.histplot(corr_within, color='#2089D5', ax = axes, kde = True, bins = 50, label = 'top-top', alpha = 0.5)
g = sns.histplot(corr_across, color='lightgrey', ax = axes, kde = True, bins = 50, label = 'top-remaining', alpha = 0.5)
g.tick_params(labelsize=12)
g.set_ylabel('number of comparisons', fontsize = 12)
g.set_xlabel('abs(correlation)', fontsize = 12)
plt.legend(prop={'size':9})
g.set_xlim(0, 1)
plt.savefig(os.path.join(SUPPLEMENTARY_FIG_DIR, 'supp_figure_9f.pdf'), bbox_inches = 'tight', dpi = 300)
plt.close()

# get overlap between static and evolution top features
ranked_features = pd.read_csv(os.path.join(BASE_DIR, 'analysis_files/feature_ranking.csv'))

overlap_type_dict = {'global': [['primary_', 'baseline', 'pre_nivo', 'on_nivo'],
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
    # keep_features = overlap_top_features.feature_name_unique.unique()[:100]
    # overlap_top_features = overlap_top_features.loc[overlap_top_features.feature_name_unique.isin(keep_features), :]
    # len(overlap_top_features.feature_name_unique.unique())
    static_ids = overlap_top_features.loc[
        overlap_top_features.comparison == 'static', 'feature_name_unique'].unique()
    evolution_ids = overlap_top_features.loc[
        overlap_top_features.comparison == 'evolution', 'feature_name_unique'].unique()

    overlap_results[overlap_type] = {'static_ids': static_ids, 'evolution_ids': evolution_ids}


# get counts of features in each category
static_ids = overlap_results['global']['static_ids']
evolution_ids = overlap_results['global']['evolution_ids']
venn2([set(static_ids), set(evolution_ids)], set_labels=('Static', 'Evolution'))
plt.savefig(os.path.join(SUPPLEMENTARY_FIG_DIR, 'supp_figure_9g.pdf'))
plt.close()
