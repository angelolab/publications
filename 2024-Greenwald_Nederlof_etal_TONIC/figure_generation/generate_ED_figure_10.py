import os
import numpy as np
import pandas as pd
import matplotlib

from python_files.utils import compare_timepoints

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

import matplotlib.pyplot as plt
import seaborn as sns


BASE_DIR = "/Volumes/Shared/Noah Greenwald/TONIC_Cohort/"
SUPPLEMENTARY_FIG_DIR = os.path.join(BASE_DIR, "supplementary_figs")
ANALYSIS_DIR = os.path.join(BASE_DIR, "analysis_files")

## 3.5 Baseline to On-nivo feature evolution ##

baseline_on_nivo_viz_dir = os.path.join(SUPPLEMENTARY_FIG_DIR, 'supp_figure_10')
harmonized_metadata = pd.read_csv(os.path.join(ANALYSIS_DIR, 'harmonized_metadata.csv'))
timepoint_features = pd.read_csv(os.path.join(ANALYSIS_DIR, 'timepoint_features_filtered.csv'))
timepoint_features = timepoint_features.merge(harmonized_metadata[['Patient_ID', 'Tissue_ID', 'Timepoint', 'baseline__on_nivo']].drop_duplicates(), on='Tissue_ID')
feature_subset = timepoint_features.loc[(timepoint_features.baseline__on_nivo) & (timepoint_features.Timepoint.isin(['baseline', 'on_nivo'])), :]

def summarize_timepoint_enrichment(input_df, feature_df, timepoints, output_dir, pval_thresh=2,
                                   diff_thresh=0.3, plot_type='strip', sort_by='mean_diff'):
    """Generate a summary of the timepoint enrichment results

    Args:
        input_df (pd.DataFrame): dataframe containing timepoint enrichment results
        feature_df (pd.DataFrame): dataframe containing feature information
        timepoints (list): list of timepoints to include
        output_dir (str): path to output directory
        pval_thresh (float): threshold for p-value
        diff_thresh (float): threshold for difference between timepoints
    """

    input_df_filtered = input_df.loc[(input_df.log_pval > pval_thresh) & (np.abs(input_df[sort_by]) > diff_thresh), :]

    input_df_filtered = input_df_filtered.sort_values(sort_by, ascending=False)

    # plot the results
    for idx, feature in enumerate(input_df_filtered.feature_name_unique):
        feature_subset = feature_df.loc[(feature_df.feature_name_unique == feature), :]
        feature_subset = feature_subset.loc[(feature_subset.Timepoint.isin(timepoints)), :]

        g = sns.catplot(data=feature_subset, x='Timepoint', y='raw_mean', kind=plot_type, color='grey')
        g.fig.suptitle(feature)
        g.savefig(os.path.join(output_dir, 'Evolution_{}_{}.pdf'.format(idx, feature)))
        plt.close()

    sns.catplot(data=input_df_filtered, x=sort_by, y='feature_name_unique', kind='bar', color='grey')
    plt.savefig(os.path.join(output_dir, 'Timepoint_summary.pdf'))
    plt.close()


for paired_status in ['baseline__on_nivo', None]:
    subdir_name = 'paired' if paired_status else 'unpaired'
    subdir_path = os.path.join(baseline_on_nivo_viz_dir, subdir_name)
    os.makedirs(subdir_path)

    primary_met_means = compare_timepoints(
        feature_df=timepoint_features, timepoint_1_name='baseline', timepoint_1_list=['baseline'],
        timepoint_2_name='on_nivo', timepoint_2_list=['on_nivo'], paired=paired_status, feature_suff='mean')

    summarize_timepoint_enrichment(input_df=primary_met_means, feature_df=timepoint_features,
                                   timepoints=['baseline', 'on_nivo'],
                                   pval_thresh=2, diff_thresh=0.3, output_dir=os.path.join(baseline_on_nivo_viz_dir, subdir_name))
