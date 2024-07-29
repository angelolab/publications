# File with code for generating supplementary plots
import os
import random
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
from ark.utils.plot_utils import cohort_cluster_plot
# from toffy import qc_comp, qc_metrics_plots
from alpineer import io_utils

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

import matplotlib.pyplot as plt
import seaborn as sns
from statsmodels.stats.multitest import multipletests


import python_files.supplementary_plot_helpers as supplementary_plot_helpers

BASE_DIR = "/Volumes/Shared/Noah Greenwald/TONIC_Cohort/"
ANALYSIS_DIR = "/Volumes/Shared/Noah Greenwald/TONIC_Cohort/analysis_files"
CHANNEL_DIR = '/Volumes/Shared/Noah Greenwald/TONIC_Cohort/image_data/samples/'
INTERMEDIATE_DIR = "/Volumes/Shared/Noah Greenwald/TONIC_Cohort/intermediate_files"
OUTPUT_DIR = "/Volumes/Shared/Noah Greenwald/TONIC_Cohort/output_files"
METADATA_DIR = "/Volumes/Shared/Noah Greenwald/TONIC_Cohort/intermediate_files/metadata"
SEG_DIR = '/Volumes/Shared/Noah Greenwald/TONIC_Cohort/segmentation_data/deepcell_output'
SUPPLEMENTARY_FIG_DIR = "/Volumes/Shared/Noah Greenwald/TONIC_Cohort/supplementary_figs"
SEQUENCE_DIR = "/Volumes/Shared/Noah Greenwald/TONIC_Cohort/sequencing_data"

# outlier patient analysis
timepoints = ['primary', 'baseline', 'pre_nivo' , 'on_nivo']

timepoint_features = pd.read_csv(os.path.join(ANALYSIS_DIR, 'timepoint_combined_features.csv'))
feature_ranking_df = pd.read_csv(os.path.join(ANALYSIS_DIR, 'feature_ranking.csv'))
feature_ranking_df = feature_ranking_df[np.isin(feature_ranking_df['comparison'], timepoints)]
feature_ranking_df = feature_ranking_df.sort_values(by = 'feature_rank_global', ascending=True)

#subset by the top response-associated features
feature_ranking_df = feature_ranking_df.iloc[:100, :]

#merge dataframes for patient-level analysis (feature, raw_mean, Patient_ID, Timepoint, Clinical_benefit) 
feature_ranking_df.rename(columns = {'comparison':'Timepoint'}, inplace = True)
merged_df = pd.merge(timepoint_features, feature_ranking_df, on = ['feature_name_unique', 'Timepoint'])

#create a dictionary mapping patients to their clinical benefit status
status_df = merged_df.loc[:, ['Patient_ID', 'Clinical_benefit']].drop_duplicates()
status_dict = dict(zip(status_df['Patient_ID'], status_df['Clinical_benefit']))

#for each feature, identify outlier patients 
outliers = dict()
feat_list = list(zip(feature_ranking_df['feature_name_unique'], feature_ranking_df['Timepoint']))
for i in range(0, len(merged_df['Clinical_benefit'].unique())): 
    for feat in feat_list:
        try:
            if i == 0:
                df_subset = merged_df.iloc[np.where((merged_df['Clinical_benefit'] == 'Yes') & (merged_df['feature_name_unique'] == feat[0]) & (merged_df['Timepoint'] == feat[1]))[0]].copy()
                df_subset_op = merged_df.iloc[np.where((merged_df['Clinical_benefit'] == 'No') & (merged_df['feature_name_unique'] == feat[0]) & (merged_df['Timepoint'] == feat[1]))[0]].copy()
            else:
                df_subset = merged_df.iloc[np.where((merged_df['Clinical_benefit'] == 'No') & (merged_df['feature_name_unique'] == feat[0]) & (merged_df['Timepoint'] == feat[1]))[0]].copy()
                df_subset_op = merged_df.iloc[np.where((merged_df['Clinical_benefit'] == 'Yes') & (merged_df['feature_name_unique'] == feat[0]) & (merged_df['Timepoint'] == feat[1]))[0]].copy()

            two_std = df_subset['raw_mean'].std() * 2

            #patient considered to be an outlier for this feature if 2 std from the mean in the direction of the opposite clinical benefit group
            outliers_indices = df_subset['raw_mean'] > df_subset['raw_mean'].mean() + two_std if df_subset_op['raw_mean'].mean() > df_subset['raw_mean'].mean() else df_subset['raw_mean'] < df_subset['raw_mean'].mean() - two_std
            outlier_patients = list(df_subset[outliers_indices]['Patient_ID'].values)

            for patient in outlier_patients:
                if patient not in outliers:
                    outliers[patient] = [(feat[0], feat[1])]
                else:
                    outliers[patient].append((feat[0], feat[1]))
        except:
            continue

#count the number of times a patient is an outlier for the top response-associated features
outlier_counts = pd.DataFrame()
counts = []
patients = []
for patient, features in outliers.items():
    counts.append(len(outliers[patient]))
    patients.append(patient.astype('int64'))

outlier_counts = pd.DataFrame(counts, index = patients, columns = ['outlier_counts'])
outlier_counts['Clinical_benefit'] = pd.Series(outlier_counts.index).map(status_dict).values

#plot the distribution indicating the number of times a patient has discordant feature values for the top response-associated features
sns.set_style('ticks')
_, axes = plt.subplots(1, 1, figsize = (5, 4), gridspec_kw={'hspace': 0.45, 'wspace': 0.4, 'bottom':0.15})
g = sns.histplot(data = outlier_counts, legend=False, ax = axes, bins = 30, color = '#AFA0BA', alpha=1)
g.tick_params(labelsize=10)
g.set_xlabel('number of outliers', fontsize = 12)
g.set_ylabel('number of patients', fontsize = 12)
plt.savefig(os.path.join(SUPPLEMENTARY_FIG_DIR, 'outlier_analysis', 'number_outliers_per_patient.pdf'), bbox_inches = 'tight', dpi =300)

#convert dictionary into a dataframe consisting of (Patient_ID, feature_name_unique, Timepoint, feature_type)
reshaped_data = []
for patient, records in outliers.items():
    for feature, timepoint in records:
        reshaped_data.append((patient, feature, timepoint))

outlier_df = pd.DataFrame(reshaped_data, columns = ['Patient_ID', 'feature_name_unique', 'Timepoint'])
outlier_df = pd.merge(outlier_df, feature_ranking_df.loc[:, ['feature_name_unique', 'feature_type', 'Timepoint']], on = ['feature_name_unique', 'Timepoint'])
outlier_df.to_csv(os.path.join(SUPPLEMENTARY_FIG_DIR, 'outlier_analysis', 'TONIC_outlier_counts.csv'), index=False)

#subset dataframe by patients that have discordant feature values in more than 4 response-associated features 
sig_outlier_patients = outlier_df.groupby('Patient_ID').size()[outlier_df.groupby('Patient_ID').size() > 4].index
sig_outlier_df = outlier_df[np.isin(outlier_df['Patient_ID'], sig_outlier_patients)].copy()

#plot the distribution of the feature classes for patients that have discordant feature values in more than 4 response-associated features 
df_pivot = sig_outlier_df.groupby(['Patient_ID', 'feature_type']).size().unstack().reset_index().melt(id_vars = 'Patient_ID').pivot(index='Patient_ID', columns='feature_type', values='value')

df_pivot.plot(kind='bar', stacked=True, figsize=(6,6))
plt.ylabel('Count', fontsize = 12)
plt.xlabel('Patient ID', fontsize = 12)
g.tick_params(labelsize=10)
plt.legend(bbox_to_anchor=(1, 1))
plt.savefig(os.path.join(SUPPLEMENTARY_FIG_DIR, 'outlier_analysis', 'feature_classes_outlier_patients.pdf'), bbox_inches = 'tight', dpi =300)

#are any features are consistently discordant for patients considered to be outliers (i.e. have discordant feature values in more than 4 response-associated features)?
cross_feat_counts = pd.DataFrame(sig_outlier_df.groupby('feature_name_unique').size().sort_values(ascending = False), columns = ['count'])
_, axes = plt.subplots(1, 1, figsize = (6, 16), gridspec_kw={'hspace': 0.45, 'wspace': 0.4, 'bottom':0.15})
g = sns.barplot(x = 'count', y =cross_feat_counts.index, data = cross_feat_counts, ax = axes, color = 'lightgrey')
g.tick_params(labelsize=10)
g.set_ylabel('Feature', fontsize = 12)
g.set_xlabel('Number of outlier patients', fontsize = 12)
plt.savefig(os.path.join(SUPPLEMENTARY_FIG_DIR, 'outlier_analysis', 'features_outlier_patients.pdf'), bbox_inches = 'tight', dpi =300)
