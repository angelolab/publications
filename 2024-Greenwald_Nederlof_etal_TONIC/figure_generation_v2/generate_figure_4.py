import os

import numpy as np
import pandas as pd
import matplotlib

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

import matplotlib.pyplot as plt
import seaborn as sns
from ark.utils.plot_utils import cohort_cluster_plot
import ark.settings as settings
import skimage.io as io


base_dir = '/Volumes/Shared/Noah Greenwald/TONIC_Cohort'
metadata_dir = os.path.join(base_dir, 'intermediate_files/metadata')
plot_dir = os.path.join(base_dir, 'figures')
seg_dir = os.path.join(base_dir, 'segmentation_data/deepcell_output')
image_dir = os.path.join(base_dir, 'image_data/samples/')

# load files
harmonized_metadata = pd.read_csv(os.path.join(metadata_dir, 'harmonized_metadata.csv'))
ranked_features_all = pd.read_csv(os.path.join(base_dir, 'analysis_files/feature_ranking.csv'))
ranked_features = ranked_features_all.loc[ranked_features_all.comparison.isin(['primary', 'baseline', 'pre_nivo', 'on_nivo'])]
top_features = ranked_features.loc[ranked_features.comparison.isin(['primary', 'baseline', 'pre_nivo', 'on_nivo']), :]
top_features = top_features.iloc[:100, :]


# summarize distribution of top features
top_features_by_comparison = top_features[['feature_name_unique', 'comparison']].groupby('comparison').count().reset_index()
top_features_by_comparison.columns = ['comparison', 'num_features']
top_features_by_comparison = top_features_by_comparison.sort_values('num_features', ascending=False)

fig, ax = plt.subplots(figsize=(4, 4))
sns.barplot(data=top_features_by_comparison, x='comparison', y='num_features', color='grey', ax=ax)
plt.xticks(rotation=90)
plt.tight_layout()
sns.despine()
plt.savefig(os.path.join(plot_dir, 'Figure4a_num_features.pdf'))
plt.close()

# heatmap of top features over time
timepoints = ['primary', 'baseline', 'pre_nivo' , 'on_nivo']

timepoint_features = pd.read_csv(os.path.join(base_dir, 'analysis_files/timepoint_combined_features.csv'))
feature_ranking_df = pd.read_csv(os.path.join(base_dir, 'analysis_files/feature_ranking.csv'))
feature_ranking_df = feature_ranking_df[np.isin(feature_ranking_df['comparison'], timepoints)]
feature_ranking_df = feature_ranking_df.sort_values(by = 'feature_rank_global', ascending=True)

#access the top response-associated features (unique because a feature could be in the top in multiple timepoints)
top_features = np.unique(feature_ranking_df.loc[:, 'feature_name_unique'][:100])

#compute the 90th percentile of importance scores and plot the distribution
perc = np.percentile(feature_ranking_df.importance_score, 90)
# _, axes = plt.subplots(1, 1, figsize = (4.5, 3.5), gridspec_kw={'hspace': 0.45, 'wspace': 0.4, 'bottom':0.15})
# g = sns.histplot(np.abs(feature_ranking_df.importance_score), ax = axes, color = '#1885F2')
# g.tick_params(labelsize=12)
# g.set_xlabel('importance score', fontsize = 12)
# g.set_ylabel('count', fontsize = 12)
# axes.axvline(perc, color = 'k', ls = '--', lw = 1, label = '90th percentile')
# g.legend(bbox_to_anchor=(0.98, 0.95), loc='upper right', borderaxespad=0, prop={'size':10})
# plt.show()

#subset data based on the 90th percentile
feature_ranking_df = feature_ranking_df[feature_ranking_df['importance_score'] > perc]

#min max scale the importance scores (scales features from 0 to 1)
from sklearn.preprocessing import MinMaxScaler
scaled_perc_scores = MinMaxScaler().fit_transform(feature_ranking_df['importance_score'].values.reshape(-1,1))
feature_ranking_df.loc[:, 'scaled_percentile_importance'] = scaled_perc_scores

#pivot the dataframe for plotting (feature x timepoint)
pivot_df = feature_ranking_df.loc[:, ['scaled_percentile_importance', 'feature_name_unique', 'comparison']].pivot(index = 'feature_name_unique', columns = 'comparison')
pivot_df.columns = pivot_df.columns.droplevel(0)
pivot_df = pivot_df.loc[:, timepoints] #reorder
pivot_df.fillna(0, inplace = True) #set features with nan importance scores (i.e. not in the top 90th percentile) to 0

#subset according to top response-associated features
pivot_df = pivot_df.loc[top_features, :]

#access the top 100 feature-timepoint pairs
pivot_df_top = feature_ranking_df[:100].loc[:, ['scaled_percentile_importance', 'feature_name_unique', 'comparison']].pivot(index = 'feature_name_unique', columns = 'comparison')
pivot_df_top.columns = pivot_df_top.columns.droplevel(0)
pivot_df_top = pivot_df_top.loc[:, timepoints] #reorder
pivot_df_top.fillna(0, inplace = True) #set features with nan importance scores (i.e. not in the top 90th percentile) to 0

#assign feature with delta label for plot order
feat_timepoint_dict = {'CD38+__Immune_Other': 'baseline',
 'diversity_cell_cluster__Cancer_2__stroma_border': 'on_nivo',
 'NK__Other__ratio__cancer_border':'baseline__pre_nivo__on_nivo',
 'CD45RO+__Fibroblast':'pre_nivo',
 'PDL1+__CD4T':'pre_nivo',
 'CD69+__CD4T':'pre_nivo',
 'TBET+__Treg':'pre_nivo',
 'diversity_cell_cluster__Monocyte__cancer_border':'pre_nivo',
 'PDL1+__Monocyte':'pre_nivo',
 'TBET+__CD4T':'pre_nivo__on_nivo',
 'CD45RO+__Monocyte':'pre_nivo',
 'PDL1+__CD163_Mac':'pre_nivo',
 'CD45RO+__CAF':'pre_nivo',
 'CD45RO+__Immune_Other':'pre_nivo',
 'CD45RO__CD69+__NK':'pre_nivo',
 'Ki67+__T_Other':'pre_nivo',
 'PDL1+__APC':'pre_nivo',
 'CD45RO+__CAF__cancer_border':'pre_nivo',
 'PDL1+__Mac_Other':'pre_nivo',
 'B__NK__ratio__cancer_core':'primary',
 'T__Cancer__ratio__cancer_border':'baseline__pre_nivo__on_nivo',
 'T__Cancer__ratio__cancer_core':'pre_nivo__on_nivo',
 'CD45RO+__CD68_Mac':'baseline__pre_nivo__on_nivo',
 'diversity_cell_cluster__Cancer_1__stroma_core':'baseline__pre_nivo__on_nivo',
 'Other__Structural__ratio':'baseline__pre_nivo__on_nivo',
 'Other__Cancer__ratio__stroma_border':'baseline__pre_nivo__on_nivo',
 'diversity_cell_cluster__Fibroblast__cancer_border':'pre_nivo__on_nivo',
 'Other__Cancer__ratio__cancer_border':'pre_nivo__on_nivo',
 'cluster_broad_diversity_cancer_border':'pre_nivo__on_nivo',
 'diversity_cell_cluster__APC':'baseline__pre_nivo__on_nivo',
 'T__Cancer__ratio__stroma_border':'baseline__on_nivo',
 'diversity_cell_cluster__Cancer_2__stroma_core':'baseline__on_nivo',
 'NK__T__ratio__cancer_border':'baseline__on_nivo',
 'Other__Cancer__ratio__stroma_core':'baseline__pre_nivo__on_nivo',
 'diversity_cell_cluster__CAF':'baseline__on_nivo',
 'cluster_3__proportion__stroma_border':'on_nivo',
 'B__Structural__ratio__cancer_border':'on_nivo',
 'Cancer_3__proportion_of__Cancer':'on_nivo',
 'T__cluster_broad_density__cancer_border':'on_nivo',
 'diversity_cell_cluster__Cancer_2':'on_nivo',
 'B__Cancer__ratio__stroma_border':'on_nivo',
 'CD8T__cluster_density__cancer_border':'on_nivo',
 'Cancer_Immune_mixing_score':'on_nivo',
 'NK__Structural__ratio__cancer_border':'on_nivo',
 'B__Cancer__ratio__cancer_border':'on_nivo',
 'diversity_cell_cluster__Cancer_1':'on_nivo',
 'NK__Structural__ratio__cancer_core':'on_nivo',
 'Other__Structural__ratio__cancer_core':'on_nivo',
 'Cancer_1__proportion_of__Cancer':'on_nivo',
 'cancer_diversity':'on_nivo',
 'CD8T__cluster_density':'on_nivo',
 'Monocyte__proportion_of__Mono_Mac__cancer_core':'on_nivo',
 'cluster_2__proportion':'on_nivo',
 'B__Cancer__ratio':'on_nivo',
 'B__Granulocyte__ratio':'on_nivo',
 'cluster_broad_diversity_cancer_core':'on_nivo',
 'CD69+__all':'on_nivo',
 'Mono_Mac__T__ratio':'on_nivo',
 'Cancer_3__proportion_of__Cancer__stroma_core':'on_nivo',
 'cancer_diversity_stroma_core':'on_nivo',
 'T__Cancer__ratio':'on_nivo',
 'Cancer_1__proportion_of__Cancer__stroma_core':'on_nivo',
 'Mono_Mac__Cancer__ratio__stroma_border':'on_nivo',
 'NK__Cancer__ratio__stroma_border':'on_nivo',
 'B__Structural__ratio':'on_nivo',
 'B__T__ratio__cancer_core':'on_nivo',
 'T_Other__cluster_density__cancer_border':'on_nivo',
 'PDL1+__Cancer_3__stroma_border':'on_nivo',
 'NK__Cancer__ratio':'on_nivo',
 'diversity_cell_cluster__Neutrophil':'on_nivo',
 'diversity_cell_cluster__Fibroblast':'on_nivo',
 'diversity_cell_cluster__Cancer_3':'on_nivo',
 'diversity_cell_cluster__Smooth_Muscle':'on_nivo',
 'Other__Cancer__ratio':'baseline__pre_nivo__on_nivo',
 'T__Structural__ratio__cancer_core':'on_nivo',
 'diversity_cell_cluster__Endothelium':'on_nivo',
 'T__Structural__ratio':'on_nivo',
 'Other__Cancer__ratio__cancer_core':'pre_nivo__on_nivo',
 'PDL1+__CAF__cancer_border':'pre_nivo__on_nivo',
 'PDL1+__CD68_Mac':'pre_nivo__on_nivo',
 'CD8T__Treg__ratio__cancer_core':'pre_nivo__on_nivo',
 'CD45RO+__all':'pre_nivo__on_nivo',
 'diversity_cell_cluster__CAF__cancer_border':'pre_nivo__on_nivo',
 'CD45RO+__Monocyte__cancer_border':'pre_nivo',
 'CD4T__cluster_density__cancer_border':'on_nivo',
 'Ki67+__Cancer_2__stroma_core':'on_nivo',
 'Mono_Mac__Structural__ratio__cancer_core': 'on_nivo',
 'Structural__cluster_broad_density__stroma_border': 'baseline__on_nivo',
 'T_Other__cluster_density':'on_nivo',
 'T__cluster_broad_density':'on_nivo',
 'stroma_core__stroma_border__log2_ratio': 'on_nivo',
 'diversity_cell_cluster__CD163_Mac__cancer_border': 'pre_nivo__on_nivo',
 'cluster_3__proportion__cancer_core': 'pre_nivo',
 'cancer_border__stroma_core__log2_ratio': 'on_nivo',
 'cancer_border__proportion': 'on_nivo',
 'TCF1+__all': 'on_nivo',
 'PDL1+__Neutrophil': 'pre_nivo',
 'PDL1+__CD8T': 'pre_nivo',
 'PD1+__all': 'on_nivo',
 'Other__cluster_broad_density': 'baseline__on_nivo',
 'NK__T__ratio__cancer_core': 'on_nivo',
 'NK__Cancer__ratio__cancer_border': 'on_nivo',
 'Mono_Mac__Cancer__ratio': 'on_nivo',
 'Ki67+__Fibroblast__cancer_border': 'on_nivo',
 'Ki67+__Cancer_1__stroma_border': 'on_nivo',
 'Immune_Other__cluster_density': 'on_nivo',
 'Cancer__Structural__ratio': 'on_nivo',
 'CD45RO+__Immune_Other__cancer_border': 'pre_nivo',
 'CD45RB+__all': 'on_nivo',
 'B__T__ratio': 'on_nivo',
 'B__Structural__ratio__cancer_core': 'on_nivo',
 'B__Other__ratio': 'on_nivo',
 'B__NK__ratio': 'on_nivo',
 'B__Mono_Mac__ratio': 'on_nivo'
}

#sort dataframe by delta group and get the order of the ticks
pivot_df["group"] = pivot_df.index.map(feat_timepoint_dict) 
pivot_df["group"] = pd.Categorical(pivot_df["group"], categories=['primary__baseline__pre_nivo__on_nivo', 'baseline__pre_nivo__on_nivo', 'pre_nivo__on_nivo', 'baseline__pre_nivo', 'baseline__on_nivo', 'primary__on_nivo', 'on_nivo', 'pre_nivo', 'baseline', 'primary'], ordered=True)
pivot_df.sort_values("group", inplace=True) 
xlabs = list(pivot_df.index)

#plot clustermap
from matplotlib.patches import Rectangle
cmap = ['#D8C198', '#D88484', '#5AA571', '#4F8CBE']
sns.set_style('ticks')

pivot_df_run = pivot_df.loc[xlabs, :].copy()
pivot_df_run.drop(columns  = ['group'], inplace=True)
pivot_df_top_run = pivot_df_top.loc[xlabs, :].copy()

g = sns.clustermap(data = pivot_df_run, yticklabels=True, cmap = 'Blues', vmin = 0, vmax = 1, row_cluster = False,
                   col_cluster = False, figsize = (5, 15), cbar_pos=(1, .03, .02, .1), dendrogram_ratio=0.1, colors_ratio=0.01,
                   col_colors=cmap)

g.tick_params(labelsize=12)

ax = g.ax_heatmap
ax.set_ylabel('Response-associated Features', fontsize = 12)
ax.set_xlabel('Timepoint', fontsize = 12)

ax.axvline(x=0, color='k',linewidth=2.5)
ax.axvline(x=1, color='k',linewidth=1.5)
ax.axvline(x=2, color='k',linewidth=1.5)
ax.axvline(x=3, color='k',linewidth=1.5)
ax.axvline(x=4, color='k',linewidth=2.5)
ax.axhline(y=0, color='k',linewidth=2.5)
ax.axhline(y=len(pivot_df), color='k',linewidth=2.5)

x0, _y0, _w, _h = g.cbar_pos
for spine in g.ax_cbar.spines:
    g.ax_cbar.spines[spine].set_color('k')
    g.ax_cbar.spines[spine].set_linewidth(1)

for i in range(0, pivot_df_top_run.shape[0]):
    row = pivot_df_top_run.astype('bool').iloc[i, :]
    ids = np.where(row == True)[0]
    for id in ids:
        #creates rectangle at given indices of top 100 feature timepoint pairs (x = timepoint_index, y = feature_index)
        rect = Rectangle((id, i), 1, 1, fill=False, edgecolor='red', lw=0.5, zorder = 10)

        # Add it to the plot
        g.ax_heatmap.add_patch(rect)

        # Redraw the figure
        plt.draw()

plt.savefig(os.path.join(plot_dir, 'Figure4b.pdf'), bbox_inches = 'tight', dpi =300)


# longitudinal T / Cancer ratios
combined_df = pd.read_csv(os.path.join(base_dir, 'analysis_files/timepoint_combined_features_outcome_labels.csv'))


# generate summary plots
for timepoint in ['primary', 'baseline', 'pre_nivo', 'on_nivo']:

    plot_df = combined_df.loc[(combined_df.feature_name_unique == 'T__Cancer__ratio__cancer_border') &
                              (combined_df.Timepoint == timepoint), :]

    fig, ax = plt.subplots(1, 1, figsize=(2, 4))
    sns.stripplot(data=plot_df, x='Clinical_benefit', y='raw_mean', order=['Yes', 'No'],
                    color='black', ax=ax)
    sns.boxplot(data=plot_df, x='Clinical_benefit', y='raw_mean', order=['Yes', 'No'],
                    color='grey', ax=ax, showfliers=False, width=0.3)
    ax.set_title('T/C ratio ' + ' ' + timepoint)
    ax.set_ylim([-15, 0])
    sns.despine()
    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir, 'Figure4c_T_C_ratio_{}.pdf'.format(timepoint)))
    plt.close()


# identify patients to show for visualization
# # check for longitudinal patients
# longitudinal_patients = combined_df.loc[combined_df.Timepoint.isin(['primary', 'baseline', 'pre_nivo', 'on_nivo',]), :]
# longitudinal_patients = longitudinal_patients.loc[longitudinal_patients.Clinical_benefit == 'Yes', :]
# longitudinal_patients = longitudinal_patients.loc[longitudinal_patients.feature_name_unique == 'T__Cancer__ratio__cancer_border', :]
#
# longitudinal_wide = longitudinal_patients.pivot(index=['Patient_ID'], columns='Timepoint', values='raw_mean')
#
#
# # corresponding overlays
cell_table_clusters = pd.read_csv(os.path.join(base_dir, 'analysis_files/cell_table_clusters.csv'))
annotations_by_mask = pd.read_csv(os.path.join(base_dir, 'intermediate_files/mask_dir', 'cell_annotation_mask.csv'))
annotations_by_mask = annotations_by_mask.rename(columns={'mask_name': 'tumor_region'})
cell_table_clusters = cell_table_clusters.merge(annotations_by_mask, on=['fov', 'label'], how='left')

#
tc_colormap = pd.DataFrame({'T_C_ratio': ['T', 'Cancer', 'Other_region', 'Other_cells'],
                         'color': ['yellow','white', 'grey', 'grey']})
#
# # pick patients for visualization
# subset = plot_df.loc[plot_df.raw_mean > -4, :]
#
# pats = [26, 33, 59, 62, 64, 65, 115, 118]
# fovs = harmonized_metadata.loc[(harmonized_metadata.Patient_ID.isin(pats) & harmonized_metadata.MIBI_data_generated.values), 'fov'].unique()
#
# # add column for T in cancer border, T elsewhere, and others
# cell_table_subset = cell_table_clusters.loc[(cell_table_clusters.fov.isin(fovs)), :]
# cell_table_subset['T_C_ratio'] = cell_table_subset.cell_cluster_broad
# cell_table_subset.loc[~cell_table_subset.T_C_ratio.isin(['T', 'Cancer']), 'T_C_ratio'] = 'Other_cells'
# cell_table_subset.loc[cell_table_subset.tumor_region != 'cancer_border', 'T_C_ratio'] = 'Other_region'
#
# #outside_t = (cell_table_subset.cell_cluster_broad == 'T') & (cell_table_subset.tumor_region != 'cancer_border')
# #outside_cancer = (cell_table_subset.cell_cluster_broad == 'Cancer') & (cell_table_subset.tumor_region != 'cancer_border')
# #cell_table_subset.loc[outside_t, 'T_C_ratio'] = 'T_outside'
# #cell_table_subset.loc[outside_cancer, 'T_C_ratio'] = 'Cancer_outside'
#
# tc_nivo_plot_dir = os.path.join(plot_dir, 'Figure5_tc_overlays_nivo')
# if not os.path.exists(tc_nivo_plot_dir):
#     os.mkdir(tc_nivo_plot_dir)
#
#
# for pat in pats:
#     pat_fovs = harmonized_metadata.loc[(harmonized_metadata.Patient_ID == pat) & (harmonized_metadata.MIBI_data_generated.values) & (harmonized_metadata.Timepoint == timepoint), 'fov'].unique()
#     pat_df = cell_table_subset.loc[cell_table_subset.fov.isin(pat_fovs), :]
#
#     pat_dir = os.path.join(tc_nivo_plot_dir, 'patient_{}'.format(pat))
#     if not os.path.exists(pat_dir):
#         os.mkdir(pat_dir)
#
#     cohort_cluster_plot(
#         fovs=pat_fovs,
#         seg_dir=seg_dir,
#         save_dir=pat_dir,
#         cell_data=pat_df,
#         erode=True,
#         fov_col=settings.FOV_ID,
#         label_col=settings.CELL_LABEL,
#         cluster_col='T_C_ratio',
#         seg_suffix="_whole_cell.tiff",
#         cmap=tc_colormap,
#         display_fig=False,
#     )
#
#
# ## check for induction FOVs
timepoint = 'pre_nivo'

plot_df = combined_df.loc[(combined_df.feature_name_unique == 'T__Cancer__ratio__cancer_border') &
                                    (combined_df.Timepoint == timepoint), :]
#
# fig, ax = plt.subplots(1, 1, figsize=(2, 4))
# sns.stripplot(data=plot_df, x='Clinical_benefit', y='raw_mean', order=['Yes', 'No'],
#                 color='black', ax=ax)
# sns.boxplot(data=plot_df, x='Clinical_benefit', y='raw_mean', order=['Yes', 'No'],
#                 color='grey', ax=ax, showfliers=False, width=0.3)
# ax.set_title(feature_name + ' ' + timepoint)
# ax.set_ylim([-15, 0])
# sns.despine()
# plt.tight_layout()
# plt.savefig(os.path.join(plot_dir, 'Figure5_feature_{}_{}.pdf'.format(feature_name, timepoint)))
# plt.close()
#
#
# # pick patients for visualization
subset = plot_df.loc[plot_df.raw_mean > -6, :]
pats = subset.loc[subset.Clinical_benefit == 'Yes', 'Patient_ID'].unique()

pats = [26, 4, 5, 11, 40, 37, 46, 56, 62, 64, 65, 102]
fovs = harmonized_metadata.loc[(harmonized_metadata.Patient_ID.isin(pats) & harmonized_metadata.MIBI_data_generated.values), 'fov'].unique()

# add column for T in cancer border, T elsewhere, and others
cell_table_subset = cell_table_clusters.loc[(cell_table_clusters.fov.isin(fovs)), :]
cell_table_subset['T_C_ratio'] = cell_table_subset.cell_cluster_broad
cell_table_subset.loc[~cell_table_subset.T_C_ratio.isin(['T', 'Cancer']), 'T_C_ratio'] = 'Other_cells'
cell_table_subset.loc[cell_table_subset.tumor_region != 'cancer_border', 'T_C_ratio'] = 'Other_region'

#outside_t = (cell_table_subset.cell_cluster_broad == 'T') & (cell_table_subset.tumor_region != 'cancer_border')
#outside_cancer = (cell_table_subset.cell_cluster_broad == 'Cancer') & (cell_table_subset.tumor_region != 'cancer_border')
#cell_table_subset.loc[outside_t, 'T_C_ratio'] = 'T_outside'
#cell_table_subset.loc[outside_cancer, 'T_C_ratio'] = 'Cancer_outside'

tc_induction_plot_dir = os.path.join(plot_dir, 'Figure5_tc_overlays_induction')
if not os.path.exists(tc_induction_plot_dir):
    os.mkdir(tc_induction_plot_dir)


for pat in pats:
    pat_fovs = harmonized_metadata.loc[(harmonized_metadata.Patient_ID == pat) & (harmonized_metadata.MIBI_data_generated.values) & (harmonized_metadata.Timepoint == timepoint), 'fov'].unique()
    pat_df = cell_table_subset.loc[cell_table_subset.fov.isin(pat_fovs), :]

    pat_dir = os.path.join(tc_induction_plot_dir, 'patient_{}'.format(pat))
    if not os.path.exists(pat_dir):
        os.mkdir(pat_dir)

    cohort_cluster_plot(
        fovs=pat_fovs,
        seg_dir=seg_dir,
        save_dir=pat_dir,
        cell_data=pat_df,
        erode=True,
        fov_col=settings.FOV_ID,
        label_col=settings.CELL_LABEL,
        cluster_col='T_C_ratio',
        seg_suffix="_whole_cell.tiff",
        cmap=tc_colormap,
        display_fig=False,
    )
#
# ## check for baseline FOVs
# timepoint = 'baseline'
#
# plot_df = combined_df.loc[(combined_df.feature_name_unique == feature_name) &
#                                     (combined_df.Timepoint == timepoint), :]
#
# fig, ax = plt.subplots(1, 1, figsize=(2, 4))
# sns.stripplot(data=plot_df, x='Clinical_benefit', y='raw_mean', order=['Yes', 'No'],
#                 color='black', ax=ax)
# sns.boxplot(data=plot_df, x='Clinical_benefit', y='raw_mean', order=['Yes', 'No'],
#                 color='grey', ax=ax, showfliers=False, width=0.3)
# ax.set_title(feature_name + ' ' + timepoint)
# ax.set_ylim([-15, 0])
# sns.despine()
# plt.tight_layout()
# plt.savefig(os.path.join(plot_dir, 'Figure5_feature_{}_{}.pdf'.format(feature_name, timepoint)))
# plt.close()
#
#
# # pick patients for visualization
# subset = plot_df.loc[plot_df.raw_mean > -6, :]
#
# pats = [26, 5, 11, 56, 64, 65, 84, 100, 102, 115, 118]
# fovs = harmonized_metadata.loc[(harmonized_metadata.Patient_ID.isin(pats) & harmonized_metadata.MIBI_data_generated.values), 'fov'].unique()
#
# # add column for T in cancer border, T elsewhere, and others
# cell_table_subset = cell_table_clusters.loc[(cell_table_clusters.fov.isin(fovs)), :]
# cell_table_subset['T_C_ratio'] = cell_table_subset.cell_cluster_broad
# cell_table_subset.loc[~cell_table_subset.T_C_ratio.isin(['T', 'Cancer']), 'T_C_ratio'] = 'Other_cells'
# cell_table_subset.loc[cell_table_subset.tumor_region != 'cancer_border', 'T_C_ratio'] = 'Other_region'
#
# #outside_t = (cell_table_subset.cell_cluster_broad == 'T') & (cell_table_subset.tumor_region != 'cancer_border')
# #outside_cancer = (cell_table_subset.cell_cluster_broad == 'Cancer') & (cell_table_subset.tumor_region != 'cancer_border')
# #cell_table_subset.loc[outside_t, 'T_C_ratio'] = 'T_outside'
# #cell_table_subset.loc[outside_cancer, 'T_C_ratio'] = 'Cancer_outside'
#
# tc_baseline_plot_dir = os.path.join(plot_dir, 'Figure5_tc_overlays_baseline')
# if not os.path.exists(tc_baseline_plot_dir):
#     os.mkdir(tc_baseline_plot_dir)
#
#
# for pat in pats:
#     pat_fovs = harmonized_metadata.loc[(harmonized_metadata.Patient_ID == pat) & (harmonized_metadata.MIBI_data_generated.values) & (harmonized_metadata.Timepoint == timepoint), 'fov'].unique()
#     pat_df = cell_table_subset.loc[cell_table_subset.fov.isin(pat_fovs), :]
#
#     pat_dir = os.path.join(tc_baseline_plot_dir, 'patient_{}'.format(pat))
#     if not os.path.exists(pat_dir):
#         os.mkdir(pat_dir)
#
#     cohort_cluster_plot(
#         fovs=pat_fovs,
#         seg_dir=seg_dir,
#         save_dir=pat_dir,
#         cell_data=pat_df,
#         erode=True,
#         fov_col=settings.FOV_ID,
#         label_col=settings.CELL_LABEL,
#         cluster_col='T_C_ratio',
#         seg_suffix="_whole_cell.tiff",
#         cmap=tc_colormap,
#         display_fig=False,
#     )
#
# ## check for primary FOVs
# timepoint = 'primary_untreated'
#
# plot_df = combined_df.loc[(combined_df.feature_name_unique == feature_name) &
#                                     (combined_df.Timepoint == timepoint), :]
#
# fig, ax = plt.subplots(1, 1, figsize=(2, 4))
# sns.stripplot(data=plot_df, x='Clinical_benefit', y='raw_mean', order=['Yes', 'No'],
#                 color='black', ax=ax)
# sns.boxplot(data=plot_df, x='Clinical_benefit', y='raw_mean', order=['Yes', 'No'],
#                 color='grey', ax=ax, showfliers=False, width=0.3)
# ax.set_title(feature_name + ' ' + timepoint)
# ax.set_ylim([-15, 0])
# sns.despine()
# plt.tight_layout()
# plt.savefig(os.path.join(plot_dir, 'Figure5_feature_{}_{}.pdf'.format(feature_name, timepoint)))
# plt.close()
#
#
# # pick patients for visualization
# # subset = plot_df.loc[plot_df.Clinical_benefit == "Yes", :]
#
# pats = [26, 59, 105, 4, 26, 11, 37, 14, 46, 62, 121, 85]
# fovs = harmonized_metadata.loc[(harmonized_metadata.Patient_ID.isin(pats) & harmonized_metadata.MIBI_data_generated.values), 'fov'].unique()
#
# # add column for T in cancer border, T elsewhere, and others
# cell_table_subset = cell_table_clusters.loc[(cell_table_clusters.fov.isin(fovs)), :]
# cell_table_subset['T_C_ratio'] = cell_table_subset.cell_cluster_broad
# cell_table_subset.loc[~cell_table_subset.T_C_ratio.isin(['T', 'Cancer']), 'T_C_ratio'] = 'Other_cells'
# cell_table_subset.loc[cell_table_subset.tumor_region != 'cancer_border', 'T_C_ratio'] = 'Other_region'
#
# #outside_t = (cell_table_subset.cell_cluster_broad == 'T') & (cell_table_subset.tumor_region != 'cancer_border')
# #outside_cancer = (cell_table_subset.cell_cluster_broad == 'Cancer') & (cell_table_subset.tumor_region != 'cancer_border')
# #cell_table_subset.loc[outside_t, 'T_C_ratio'] = 'T_outside'
# #cell_table_subset.loc[outside_cancer, 'T_C_ratio'] = 'Cancer_outside'
#
# tc_primary_plot_dir = os.path.join(plot_dir, 'Figure5_tc_overlays_primary')
# if not os.path.exists(tc_primary_plot_dir):
#     os.mkdir(tc_primary_plot_dir)
#
#
# for pat in pats:
#     pat_fovs = harmonized_metadata.loc[(harmonized_metadata.Patient_ID == pat) & (harmonized_metadata.MIBI_data_generated.values) & (harmonized_metadata.Timepoint == timepoint), 'fov'].unique()
#     pat_df = cell_table_subset.loc[cell_table_subset.fov.isin(pat_fovs), :]
#
#     pat_dir = os.path.join(tc_primary_plot_dir, 'patient_{}'.format(pat))
#     if not os.path.exists(pat_dir):
#         os.mkdir(pat_dir)
#
#     cohort_cluster_plot(
#         fovs=pat_fovs,
#         seg_dir=seg_dir,
#         save_dir=pat_dir,
#         cell_data=pat_df,
#         erode=True,
#         fov_col=settings.FOV_ID,
#         label_col=settings.CELL_LABEL,
#         cluster_col='T_C_ratio',
#         seg_suffix="_whole_cell.tiff",
#         cmap=tc_colormap,
#         display_fig=False,
#     )


# generate crops for selected FOVs

# nivo: 33, 65, 115
# pre nivo: 37
# baseline: 26
# primary: 4, 11, 37
fovs = ['TONIC_TMA12_R5C6', 'TONIC_TMA2_R11C6', 'TONIC_TMA5_R5C2', 'TONIC_TMA2_R8C4'] # 65 (nivo), 5 (pre nivo), 26 (baseline), 4 (primary)
#fovs = ['TONIC_TMA12_R6C2', 'TONIC_TMA11_R7C6', 'TONIC_TMA11_R8C1', 'TONIC_TMA7_R10C4', 'TONIC_TMA2_R11C6']
cell_table_clusters = pd.read_csv(os.path.join(base_dir, 'analysis_files/cell_table_clusters.csv'))
annotations_by_mask = pd.read_csv(os.path.join(base_dir, 'intermediate_files/mask_dir', 'cell_annotation_mask.csv'))
annotations_by_mask = annotations_by_mask.rename(columns={'mask_name': 'tumor_region'})
cell_table_clusters = cell_table_clusters.merge(annotations_by_mask, on=['fov', 'label'], how='left')

# add column for border location
cell_table_subset = cell_table_clusters.loc[(cell_table_clusters.fov.isin(fovs)), :]
cell_table_subset['T_C_ratio'] = cell_table_subset.cell_cluster_broad
cell_table_subset.loc[~cell_table_subset.T_C_ratio.isin(['T', 'Cancer']), 'T_C_ratio'] = 'Other_cells'
cell_table_subset.loc[cell_table_subset.tumor_region != 'cancer_border', 'T_C_ratio'] = 'Other_region'

tc_colormap = pd.DataFrame({'T_C_ratio': ['T', 'Cancer', 'Other_region', 'Other_cells'],
                         'color': ['yellow','white', 'grey', 'grey']})


subset_dir = os.path.join(plot_dir, 'Figure4c_tc_ratio_overlays')
if not os.path.exists(subset_dir):
    os.mkdir(subset_dir)

cohort_cluster_plot(
    fovs=fovs,
    seg_dir=seg_dir,
    save_dir=subset_dir,
    cell_data=cell_table_subset,
    erode=True,
    fov_col=settings.FOV_ID,
    label_col=settings.CELL_LABEL,
    cluster_col='T_C_ratio',
    seg_suffix="_whole_cell.tiff",
    cmap=tc_colormap,
    display_fig=False,
)


# same thing for compartment masks
compartment_colormap = pd.DataFrame({'tumor_region': ['cancer_core', 'cancer_border', 'stroma_border', 'stroma_core', 'immune_agg'],
                         'color': ['blue', 'deepskyblue', 'lightcoral', 'firebrick', 'firebrick']})
subset_mask_dir = os.path.join(plot_dir, 'Figure4c_tc_overlays_masks')
if not os.path.exists(subset_mask_dir):
    os.mkdir(subset_mask_dir)

cohort_cluster_plot(
    fovs=fovs,
    seg_dir=seg_dir,
    save_dir=subset_mask_dir,
    cell_data=cell_table_subset,
    erode=True,
    fov_col=settings.FOV_ID,
    label_col=settings.CELL_LABEL,
    cluster_col='tumor_region',
    seg_suffix="_whole_cell.tiff",
    cmap=compartment_colormap,
    display_fig=False,
)
# scale bars: 100 um = 2048 pixels / 8 = 256 pixels, 100um = 256 / 600 = 0.426666666666666 of image

# crop overlays
fov1 = fovs[0]
row_start, col_start = 1350, 200
row_len, col_len = 600, 900

for dir in [subset_dir, subset_mask_dir]:
    fov1_image = io.imread(os.path.join(dir, 'cluster_masks_colored', fov1 + '.tiff'))
    fov1_image = fov1_image[row_start:row_start + row_len, col_start:col_start + col_len, :]
    io.imsave(os.path.join(dir, 'cluster_masks_colored', fov1 + '_crop.tiff'), fov1_image)

fov2 = fovs[1]
row_start, col_start = 1200, 250
row_len, col_len = 600, 900

for dir in [subset_dir, subset_mask_dir]:
    fov2_image = io.imread(os.path.join(dir, 'cluster_masks_colored', fov2 + '.tiff'))
    fov2_image = fov2_image[row_start:row_start + row_len, col_start:col_start + col_len, :]
    io.imsave(os.path.join(dir, 'cluster_masks_colored', fov2 + '_crop.tiff'), fov2_image)

fov3 = fovs[2]
row_start, col_start = 700, 450
row_len, col_len = 900, 600

for dir in [subset_dir, subset_mask_dir]:
    fov3_image = io.imread(os.path.join(dir, 'cluster_masks_colored', fov3 + '.tiff'))
    fov3_image = fov3_image[row_start:row_start + row_len, col_start:col_start + col_len, :]
    io.imsave(os.path.join(dir, 'cluster_masks_colored', fov3 + '_crop.tiff'), fov3_image)

fov4 = fovs[3]
row_start, col_start = 1200, 750
row_len, col_len = 600, 900

for dir in [subset_dir, subset_mask_dir]:
    fov4_image = io.imread(os.path.join(dir, 'cluster_masks_colored', fov4 + '.tiff'))
    fov4_image = fov4_image[row_start:row_start + row_len, col_start:col_start + col_len, :]
    io.imsave(os.path.join(dir, 'cluster_masks_colored', fov4 + '_crop.tiff'), fov4_image)