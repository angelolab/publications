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

# read files
ranked_features_univariate = pd.read_csv(os.path.join(base_dir, 'analysis_files/feature_ranking.csv'))
ranked_features_univariate_genomic = pd.read_csv(os.path.join(base_dir, 'sequencing_data/genomics_outcome_ranking.csv'))
cv_scores = pd.read_csv(os.path.join(base_dir, 'multivariate_lasso', 'formatted_cv_scores.csv'))


# compare AUCs
ttest_ind(cv_scores.loc[np.logical_and(cv_scores.assay == 'RNA', cv_scores.variable == 'baseline'), 'value'],
          cv_scores.loc[np.logical_and(cv_scores.assay == 'RNA', cv_scores.variable == 'pre_nivo'), 'value'])


fig, ax = plt.subplots(1, 1, figsize=(6, 4))
order = ['primary', 'baseline', 'pre_nivo','on_nivo']
sns.stripplot(data=cv_scores, x='variable', y='value', hue='assay',
              order=order, ax=ax, dodge=True)
sns.boxplot(data=cv_scores, x='variable', y='value', hue='assay',
            order=order, ax=ax, showfliers=False)

ax.set_title('AUC')
ax.set_ylim([0, 1])
ax.set_xticklabels(ax.get_xticklabels(), rotation=90)

sns.despine()
plt.tight_layout()
plt.savefig(os.path.join(plot_dir, 'Figure5a.pdf'))
plt.close()


# compute CDF for model weights based on number of required channels
all_model_rankings = pd.read_csv(os.path.join(base_dir, 'multivariate_lasso/intermediate_results', 'all_model_rankings.csv'))

mibi_rankings = all_model_rankings.loc[np.logical_and(all_model_rankings.modality == 'MIBI',
                                       all_model_rankings.timepoint == 'on_nivo'), :]

# annotate required channels per feature
channels_reqs = {
    'Cancer_3__proportion_of__Cancer': ['ECAD', 'CK17', 'H3K27'],
    'T_Other__cell_cluster_density__cancer_border': ['CD3', 'CD4', 'CD8', 'CD45'],
    'Fe+__all': ['Fe'],
    'B__Structural__ratio__cancer_border': ['CD20', 'Collagen', 'Fibronectin', 'FAP', 'SMA'],
    'Smooth_Muscle__cell_cluster_density': ['SMA'],
    'NK__Structural__ratio__cancer_core': ['CD56'],
    'TCF1+__Cancer_1': ['TCF1'],
    'Cancer_diversity': ['ECAD',],
    'B__Cancer__ratio': ['CD20', 'ECAD', 'CK17'],
    'APC__cell_cluster_density__cancer_border': ['CD11c', 'HLA-DR'],
    'Endothelium__cell_cluster_density__cancer_border': ['CD31'],
    'Ki67+__Cancer_2__stroma_core': ['Ki67'],
    'cluster_2__proportion': ['FOXP3', 'CD14', 'CD68', 'CD163', 'ChyTr', 'Calprotectin']
}

channel_counts = mibi_rankings[['feature_name_unique', 'coef_norm', 'top_ranked']].copy()
channel_counts = channel_counts.sort_values('coef_norm', ascending=False)
channel_counts['feature_num'] = np.arange(channel_counts.shape[0]) + 1

total_channels = []
for feature in channel_counts.feature_name_unique.values:
    if feature in channels_reqs:
        total_channels += channels_reqs[feature]

    channel_counts.loc[channel_counts.feature_name_unique == feature, 'total_channels'] = len(set(total_channels))

channel_counts['coef_cdf'] = channel_counts['coef_norm'].cumsum() / channel_counts['coef_norm'].sum()

# plot
fig, ax = plt.subplots(1, 1, figsize=(3, 4))
sns.lineplot(data=channel_counts, x='total_channels', y='coef_cdf', ax=ax, estimator=None, errorbar=None)
ax.set_title('CDF of model weights based on required channels')
ax.set_ylim([0, 1])
ax.set_xlim([0, 30])

plt.tight_layout()
sns.despine()
plt.savefig(os.path.join(plot_dir, 'Figure5d.pdf'))
plt.close()


# barchart with number of required channels per feature
feature_names = ['total cd38', 'structural in stroma border', 'PDL1 APC', 'Cancer 3', 'B/Structural', 'CD8T border']
channel_counts = [3, 5, 5, 3, 6, 5]

fig, ax = plt.subplots(1, 1, figsize=(3, 4))
sns.barplot(y=feature_names, x=channel_counts, ax=ax, color='grey')

ax.set_title('Number of channels per feature')
plt.tight_layout()
sns.despine()
plt.savefig(os.path.join(plot_dir, 'Figure5b.pdf'))
plt.close()

# barchart with number of transcripts per signature
file_path = os.path.join(base_dir, 'sequencing_data/preprocessing/tme_gene_signatures.gmt')

# Initialize an empty list to store rows
rows = []

# Open the file and read it line by line
with open(file_path, 'r') as file:
    for line in file:
        # Split each line by tabs
        data = line.strip().split('\t')
        # Extract the header and genes
        header = data[0]
        genes = data[2:]
        # Append the header and genes as a list to rows
        rows.append([header] + genes)

# Create a DataFrame from the list of rows
gene_counts = pd.DataFrame(rows)

gene_counts = gene_counts.T
gene_counts.columns = gene_counts.iloc[0, :]
gene_counts = gene_counts.iloc[1:, :]
gene_counts['cytolytic activity'] = ['GZMA', 'PRF1'] + [None] * (len(gene_counts) - 2)

signatures = ['Coactivation_molecules', 'Th1_signature', 'Matrix', 'cytolytic activity',
              'Antitumor_cytokines', 'M1_signatures']

transcript_counts = [np.sum(gene_counts[x].values != None) for x in signatures]

fig, ax = plt.subplots(1, 1, figsize=(3, 4))
sns.barplot(y=signatures, x=transcript_counts, ax=ax, color='grey')

# set xlim
ax.set_xlim([0, 20])

ax.set_title('Number of transcripts per signature')
plt.tight_layout()
sns.despine()
plt.savefig(os.path.join(plot_dir, 'Figure5c.pdf'))
plt.close()
