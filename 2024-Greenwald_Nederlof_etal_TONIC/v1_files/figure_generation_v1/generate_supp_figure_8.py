import os

import numpy as np
import pandas as pd
import matplotlib


matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

import matplotlib.pyplot as plt
import seaborn as sns


BASE_DIR = "/Volumes/Shared/Noah Greenwald/TONIC_Cohort/"
SUPPLEMENTARY_FIG_DIR = os.path.join(BASE_DIR, "supplementary_figs")

# all feature plot by compartment
timepoint_df = pd.read_csv(os.path.join(BASE_DIR, 'analysis_files/timepoint_features.csv'))
timepoint_df['long_name'] = timepoint_df['Tissue_ID'] + '//' + timepoint_df['feature_name']

# subset df
t = timepoint_df.pivot(index='long_name', columns='compartment')['raw_mean']
t = t[t.isnull().sum(axis=1) < 4]
t = t[~t['all'].isna()]

# 2^x for previous log2 scores so that there aren't negative values
t[np.logical_or(t.index.str.contains('__ratio'), t.index.str.contains('H3K9ac_H3K27me3_ratio+'),
                t.index.str.contains('CD45RO_CD45RB_ratio+'))] =\
    2 ** t[np.logical_or(t.index.str.contains('__ratio'), t.index.str.contains('H3K9ac_H3K27me3_ratio+'),
                         t.index.str.contains('CD45RO_CD45RB_ratio+'))]

# normalize
comp_t = t.divide(t['all'], axis=0)
comp_t.index = [idx.split('//')[1] for idx in comp_t.index]
comp_t['feature_name'] = comp_t.index

df = comp_t.groupby(by=['feature_name']).mean()
df = np.log2(df)
df.replace([np.inf, -np.inf], np.nan, inplace=True)
df = df.dropna()
df = df[['all', 'cancer_core', 'cancer_border', 'stroma_border', 'stroma_core']]

# sns.set(font_scale=1)
# plt.figure(figsize=(8, 30))
# heatmap = sns.clustermap(
#     df, cmap="vlag", vmin=-2, vmax=2, col_cluster=False, cbar_pos=(1.03, 0.07, 0.015, 0.2),
# )
# heatmap.tick_params(labelsize=8)
# plt.setp(heatmap.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
# x0, _y0, _w, _h = heatmap.cbar_pos
# for spine in heatmap.ax_cbar.spines:
#     heatmap.ax_cbar.spines[spine].set_color('k')
#     heatmap.ax_cbar.spines[spine].set_linewidth(1)
#
# ax = heatmap.ax_heatmap
# ax.axvline(x=0, color='k', linewidth=0.8)
# ax.axvline(x=1, color='k', linewidth=0.8)
# ax.axvline(x=2, color='k', linewidth=0.8)
# ax.axvline(x=3, color='k', linewidth=0.8)
# ax.axvline(x=4, color='k', linewidth=0.8)
# ax.axvline(x=5, color='k', linewidth=0.8)
# ax.axhline(y=0, color='k', linewidth=1)
# ax.axhline(y=len(df), color='k', linewidth=1.5)
# ax.set_ylabel("Feature")
# ax.set_xlabel("Compartment")
#
# features_of_interest = [361, 107, 92, 110, 90, 258, 373, 311, 236, 266, 385, 83, 327, 61, 132, 150]
# feature_names = [df.index[i] for i in features_of_interest]
# reorder = heatmap.dendrogram_row.reordered_ind
# new_positions = [reorder.index(i) for i in features_of_interest]
# plt.setp(heatmap.ax_heatmap.yaxis.set_ticks(new_positions))
# plt.setp(heatmap.ax_heatmap.yaxis.set_ticklabels(feature_names))
# plt.tight_layout()
#
# plt.savefig(os.path.join(SUPPLEMENTARY_FIG_DIR, 'supp_figure_8a.pdf'), dpi=300, bbox_inches="tight")
#
# high/low standard deviation feature plot
df_copy = df.copy()
df_copy['row_std'] = df_copy.std(axis=1)
df_copy = df_copy.sort_values(by='row_std')
df_copy.to_csv(os.path.join(SUPPLEMENTARY_FIG_DIR, 'supp_figure_8a_values.csv'))

low_std = df_copy[:90]
high_std = df_copy[-90:]
all_std_data = pd.concat([high_std, low_std]).sort_values(by='row_std', ascending=False)
all_std_data = all_std_data[df.columns]

# sns.set(font_scale=1)
# plt.figure(figsize=(4, 17))
# heatmap = sns.heatmap(
#     all_std_data, cmap="vlag", vmin=-2, vmax=2, yticklabels=True, cbar_kws={'shrink': 0.1}
# )
# heatmap.tick_params(labelsize=6)
# heatmap.hlines([len(all_std_data)/2], *ax.get_xlim(), ls='--', color='black', linewidth=0.5,)
# ax.set_ylabel("Feature")
# ax.set_xlabel("Compartment")
# plt.tight_layout()
#
# plt.savefig(os.path.join(SUPPLEMENTARY_FIG_DIR, 'compartments-high_low_std.pdf'), dpi=300, bbox_inches="tight")

row_names = ['Smooth_Muscle__cluster_density', 'T__Cancer__ratio', 'Granulocyte__T__ratio', 'CD68_Mac__CD163_Mac__ratio',
             'CD8T__CD4T__ratio', 'HLA1+__CD4T', 'Ki67+__Cancer_2', 'Ki67+__Cancer_1', 'PDL1+__CD68_Mac', 'TBET+__Treg',
             'all__total_density', 'CD38+__Endothelium']
# only show row names for selected features
selected_std_data = all_std_data.copy()
selected_std_data = selected_std_data.iloc[:, 1:]

# set index to 0 for rows not in selected features
selected_std_data.index = np.where(selected_std_data.index.isin(row_names), selected_std_data.index, 0)

selected_std_data_1 = selected_std_data.iloc[:120]
selected_std_data_2 = selected_std_data.iloc[120:]

for df, name in zip([selected_std_data_1, selected_std_data_2], ['high', 'low']):
    sns.set(font_scale=1)
    plt.figure(figsize=(4, 17))
    heatmap = sns.heatmap(
        df, cmap="vlag", vmin=-2, vmax=2, yticklabels=True, cbar_kws={'shrink': 0.1}
    )
    heatmap.tick_params(labelsize=6)
    #heatmap.hlines([len(all_std_data)/2], *ax.get_xlim(), ls='--', color='black', linewidth=0.5,)
    #ax.set_ylabel("Feature")
    #ax.set_xlabel("Compartment")
    plt.tight_layout()

    plt.savefig(os.path.join(SUPPLEMENTARY_FIG_DIR, 'supp_figure_8a_{}.pdf'.format(name)), dpi=300, bbox_inches="tight")
    plt.close()

# histogram of standard deviations
plt.style.use("default")
g = sns.histplot(df_copy.row_std)
g.set(xlabel='Standard Deviation', ylabel='Feature Counts')
sns.despine()
plt.savefig(os.path.join(SUPPLEMENTARY_FIG_DIR, 'supp_figure_8c.pdf'), dpi=300, bbox_inches="tight")
plt.close()

# plot relative scores across specific features
feature_name_plot = 'CD8T__CD4T__ratio'
feature_plot = comp_t.loc[comp_t.feature_name == feature_name_plot, :]
feature_plot = feature_plot.drop(columns=['feature_name', 'all'])
feature_plot_long = feature_plot.melt(value_name='relative_score', var_name='compartment')

fig, ax = plt.subplots(1, 1, figsize=(4, 4))
sns.boxplot(data=feature_plot_long, x='compartment', y='relative_score',
            order=['cancer_core', 'cancer_border', 'stroma_border', 'stroma_core'],
            color='grey', ax=ax, showfliers=False, width=0.3)

ax.set_title('Relative scores for {}'.format(feature_name_plot))
sns.despine()
ax.set_ylim([0, 20])
plt.tight_layout()
plt.savefig(os.path.join(SUPPLEMENTARY_FIG_DIR, 'supp_figure_8d.pdf'), dpi=300, bbox_inches="tight")
plt.close()

# plot relative scores across specific features
feature_name_plot = 'Smooth_Muscle__cluster_density'
feature_plot = comp_t.loc[comp_t.feature_name == feature_name_plot, :]
feature_plot = feature_plot.drop(columns=['feature_name', 'all'])
feature_plot_long = feature_plot.melt(value_name='relative_score', var_name='compartment')

fig, ax = plt.subplots(1, 1, figsize=(4, 4))
sns.boxplot(data=feature_plot_long, x='compartment', y='relative_score',
            order=['cancer_core', 'cancer_border', 'stroma_border', 'stroma_core'],
            color='grey', ax=ax, showfliers=False, width=0.3)

ax.set_title('Relative scores for {}'.format(feature_name_plot))
sns.despine()
ax.set_ylim([0, 6])
plt.tight_layout()
plt.savefig(os.path.join(SUPPLEMENTARY_FIG_DIR, 'supp_figure_8e.pdf'), dpi=300, bbox_inches="tight")
plt.close()

feature_name_plot = 'Ki67+__Cancer_1'
feature_plot = comp_t.loc[comp_t.feature_name == feature_name_plot, :]
feature_plot = feature_plot.drop(columns=['feature_name', 'all'])
feature_plot_long = feature_plot.melt(value_name='relative_score', var_name='compartment')

fig, ax = plt.subplots(1, 1, figsize=(4, 4))
sns.boxplot(data=feature_plot_long, x='compartment', y='relative_score',
            order=['cancer_core', 'cancer_border', 'stroma_border', 'stroma_core'],
            color='grey', ax=ax, showfliers=False, width=0.3)

ax.set_title('Relative scores for {}'.format(feature_name_plot))
sns.despine()
ax.set_ylim([0, 2])
plt.tight_layout()
plt.savefig(os.path.join(SUPPLEMENTARY_FIG_DIR, 'supp_figure_8f.pdf'), dpi=300, bbox_inches="tight")
plt.close()

# CD57+__Cancer_1 (0.86), T__distance_to__Cancer (0.76), CD38+__Endothelium (0.62), TBET+__CD4T (0.5), all__total_density (0.4), TBET+__Treg (0.33),
# Ki67+__Cancer_1 (0.27), HLA1+__Cancer_1 (0.18), PDL1+__CD68_Mac (0.12)

feature_name_plot = 'PDL1+__CD68_Mac'
feature_plot = comp_t.loc[comp_t.feature_name == feature_name_plot, :]
feature_plot = feature_plot.drop(columns=['feature_name', 'all'])
feature_plot_long = feature_plot.melt(value_name='relative_score', var_name='compartment')

fig, ax = plt.subplots(1, 1, figsize=(4, 4))
sns.boxplot(data=feature_plot_long, x='compartment', y='relative_score',
            order=['cancer_core', 'cancer_border', 'stroma_border', 'stroma_core'],
            color='grey', ax=ax, showfliers=False, width=0.3)

ax.set_title('Relative scores for {}'.format(feature_name_plot))
sns.despine()
ax.set_ylim([0, 2])
plt.tight_layout()
plt.savefig(os.path.join(SUPPLEMENTARY_FIG_DIR, 'supp_figure_8g.pdf'), dpi=300, bbox_inches="tight")
plt.close()
