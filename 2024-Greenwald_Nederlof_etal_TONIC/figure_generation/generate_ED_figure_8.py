import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


BASE_DIR = "/Volumes/Shared/Noah Greenwald/TONIC_Cohort/"
SUPPLEMENTARY_FIG_DIR = os.path.join(BASE_DIR, "supplementary_figs")

# plot key features from primary tumors
combined_df = pd.read_csv(os.path.join(BASE_DIR, 'analysis_files/timepoint_combined_features_with_outcomes.csv'))

for timepoint in ['primary', 'baseline', 'pre_nivo', 'on_nivo']:
    for feature, plot_name, lims in zip(['B__NK__ratio__cancer_core', 'CD68_Mac__cluster_density', 'Vim+__CD4T'],
                                        ['8a', '8b', '8c'], [[-5, 7], [0, 0.1], [0, 1.2]]):

        plot_df = combined_df.loc[(combined_df.feature_name_unique == feature) &
                                  (combined_df.Timepoint == timepoint), :]

        fig, ax = plt.subplots(1, 1, figsize=(2, 4))
        sns.stripplot(data=plot_df, x='Clinical_benefit', y='raw_mean', order=['Yes', 'No'],
                        color='black', ax=ax)
        sns.boxplot(data=plot_df, x='Clinical_benefit', y='raw_mean', order=['Yes', 'No'],
                        color='grey', ax=ax, showfliers=False, width=0.3)
        ax.set_title(feature + ' ' + timepoint)
        ax.set_ylim(lims)
        sns.despine()
        plt.tight_layout()
        plt.savefig(os.path.join(SUPPLEMENTARY_FIG_DIR, 'ED_figure_{}_{}.pdf'.format(plot_name, timepoint)))
        plt.close()

# top features from baseline
for timepoint in ['primary', 'baseline', 'pre_nivo', 'on_nivo']:
    for feature, plot_name, lims in zip(['Structural__cluster_broad_density__stroma_border', 'Other__Cancer__ratio__stroma_core', 'NK__Other__ratio__stroma_border'],
                                        ['8d', '8e', '8f'], [[0, 1], [-11, 4], [-10, 0]]):

        plot_df = combined_df.loc[(combined_df.feature_name_unique == feature) &
                                  (combined_df.Timepoint == timepoint), :]

        fig, ax = plt.subplots(1, 1, figsize=(2, 4))
        sns.stripplot(data=plot_df, x='Clinical_benefit', y='raw_mean', order=['Yes', 'No'],
                        color='black', ax=ax)
        sns.boxplot(data=plot_df, x='Clinical_benefit', y='raw_mean', order=['Yes', 'No'],
                        color='grey', ax=ax, showfliers=False, width=0.3)
        ax.set_title(feature + ' ' + timepoint)
        ax.set_ylim(lims)
        sns.despine()
        plt.tight_layout()
        plt.savefig(os.path.join(SUPPLEMENTARY_FIG_DIR, 'ED_figure_{}_{}.pdf'.format(plot_name, timepoint)))
        plt.close()
