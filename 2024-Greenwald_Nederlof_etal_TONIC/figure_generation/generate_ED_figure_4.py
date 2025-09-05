import os
import numpy as np
import pandas as pd
import matplotlib

from sklearn.pipeline import make_pipeline
from sklearn.cluster import KMeans
from sklearn import preprocessing

from python_files.utils import generate_crop_sum_dfs, normalize_by_ecm_area

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

import matplotlib.pyplot as plt
import seaborn as sns


BASE_DIR = "/Volumes/Shared/Noah Greenwald/TONIC_Cohort/"
SUPPLEMENTARY_FIG_DIR = os.path.join(BASE_DIR, "supplementary_figs")
ANALYSIS_DIR = os.path.join(BASE_DIR, "analysis_files")
INTERMEDIATE_DIR = os.path.join(BASE_DIR, "intermediate_files")

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
# plt.savefig(os.path.join(SUPPLEMENTARY_FIG_DIR, 'ED_figure_4a.pdf'), dpi=300, bbox_inches="tight")
#
# high/low standard deviation feature plot
df_copy = df.copy()
df_copy['row_std'] = df_copy.std(axis=1)
df_copy = df_copy.sort_values(by='row_std')
df_copy.to_csv(os.path.join(SUPPLEMENTARY_FIG_DIR, 'ED_figure_4a_values.csv'))

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

    plt.savefig(os.path.join(SUPPLEMENTARY_FIG_DIR, 'ED_figure_4a_{}.pdf'.format(name)), dpi=300, bbox_inches="tight")
    plt.close()

# histogram of standard deviations
plt.style.use("default")
g = sns.histplot(df_copy.row_std)
g.set(xlabel='Standard Deviation', ylabel='Feature Counts')
sns.despine()
plt.savefig(os.path.join(SUPPLEMENTARY_FIG_DIR, 'ED_figure_4c.pdf'), dpi=300, bbox_inches="tight")
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
plt.savefig(os.path.join(SUPPLEMENTARY_FIG_DIR, 'ED_figure_4d.pdf'), dpi=300, bbox_inches="tight")
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
plt.savefig(os.path.join(SUPPLEMENTARY_FIG_DIR, 'ED_figure_4e.pdf'), dpi=300, bbox_inches="tight")
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
plt.savefig(os.path.join(SUPPLEMENTARY_FIG_DIR, 'ED_figure_4f.pdf'), dpi=300, bbox_inches="tight")
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
plt.savefig(os.path.join(SUPPLEMENTARY_FIG_DIR, 'ED_figure_4g.pdf'), dpi=300, bbox_inches="tight")
plt.close()


# ECM stats
ecm_dir = os.path.join(INTERMEDIATE_DIR, 'ecm')
tiled_crops = pd.read_csv(os.path.join(ecm_dir, 'tiled_crops.csv'))
harmonized_metadata = pd.read_csv(os.path.join(ANALYSIS_DIR, 'harmonized_metadata.csv'))

# plot distribution of clusters in each fov
cluster_counts = tiled_crops.groupby('fov').value_counts(['tile_cluster'])
cluster_counts = cluster_counts.reset_index()
cluster_counts.columns = ['fov', 'tile_cluster', 'count']
cluster_counts = cluster_counts.pivot(index='fov', columns='tile_cluster', values='count')
cluster_counts = cluster_counts.fillna(0)
cluster_counts = cluster_counts.apply(lambda x: x / x.sum(), axis=1)

# plot the cluster counts
cluster_counts_clustermap = sns.clustermap(cluster_counts, cmap='Reds', figsize=(10, 10))
plt.close()

# use kmeans to cluster fovs
fov_kmeans_pipe = make_pipeline(KMeans(n_clusters=5, random_state=0))
fov_kmeans_pipe.fit(cluster_counts.values)
kmeans_preds = fov_kmeans_pipe.predict(cluster_counts.values)

# get average cluster count in each kmeans cluster
cluster_counts['fov_cluster'] = kmeans_preds
rename_dict = {0: 'No_ECM', 1: 'Hot_Coll', 2: 'Cold_Coll', 3: 'Hot_Cold_Coll', 4: 'Hot_No_ECM'}
cluster_counts['fov_cluster'] = cluster_counts['fov_cluster'].replace(rename_dict)
cluster_means = cluster_counts.groupby('fov_cluster').mean()
cluster_means_clustermap = sns.clustermap(cluster_means, cmap='Reds', figsize=(10, 10))
plt.close()

# correlate ECM subtypes with patient data using hierarchically clustered heatmap as index
test_fov = 'TONIC_TMA10_R10C6'
test_fov_idx = np.where(cluster_counts.index == test_fov)[0][0]
stop_idx = np.where(cluster_counts_clustermap.dendrogram_row.reordered_ind == test_fov_idx)[0][0]
cluster_counts_subset = cluster_counts.iloc[cluster_counts_clustermap.dendrogram_row.reordered_ind[:stop_idx + 1], :]

test_fov_2 = 'TONIC_TMA11_R10C6'
test_fov_idx_2 = np.where(cluster_counts.index == test_fov_2)[0][0]
stop_idx_2 = np.where(cluster_counts_clustermap.dendrogram_row.reordered_ind == test_fov_idx_2)[0][0]

cluster_counts_subset_2 = cluster_counts.iloc[cluster_counts_clustermap.dendrogram_row.reordered_ind[stop_idx:stop_idx_2 + 1], :]
cluster_counts_subset_3 = cluster_counts.iloc[cluster_counts_clustermap.dendrogram_row.reordered_ind[stop_idx_2:], :]

harmonized_metadata['ecm_cluster'] = 'inflamed'
harmonized_metadata.loc[harmonized_metadata.fov.isin(cluster_counts_subset.index), 'ecm_cluster'] = 'no_ecm'
harmonized_metadata.loc[harmonized_metadata.fov.isin(cluster_counts_subset_2.index), 'ecm_cluster'] = 'cold_collagen'
df = harmonized_metadata[['Tissue_ID', 'Localization', 'ecm_cluster']].groupby(by=['Localization', 'ecm_cluster']).count().reset_index()
df['Tissue_ID'] = df['Tissue_ID'] / df.groupby('Localization')['Tissue_ID'].transform('sum')

# # plot the distribution of ECM subtypes in each patient by localization
# g = sns.barplot(x='Localization', y='Tissue_ID', hue='ecm_cluster', data=df)
# plt.xticks(rotation=45)
# plt.title('ECM cluster by localization')
# plt.ylim(0, 1)
# plt.ylabel('Proportion')
# sns.despine()
# plt.savefig(os.path.join(ecm_cluster_viz_dir, 'ECM_cluster_localization.pdf'), bbox_inches='tight', dpi=300)

# plot marker expression in each ECM subtype
core_df_func = pd.read_csv(os.path.join(os.path.join(BASE_DIR, 'output_files'), 'functional_df_per_core_filtered_deduped.csv'))
plot_df = core_df_func.merge(harmonized_metadata[['fov', 'ecm_cluster']], on='fov', how='left')
plot_df = plot_df[plot_df.Timepoint == 'primary']


# look at all fibroblasts
temp_df = plot_df[plot_df.subset == 'all']
temp_df = temp_df[temp_df.metric == 'cluster_broad_freq']
temp_df = temp_df[temp_df.cell_type == 'Structural']
temp_df = temp_df[temp_df.functional_marker.isin(['TIM3'])]

sns.set_style('white')
sns.set_style('ticks')
g = sns.catplot(x='ecm_cluster', y='value', data=temp_df,
                kind='box', col='functional_marker', sharey=False, color='grey', showfliers=False)
g.map_dataframe(sns.stripplot, x="ecm_cluster", y="value", data=temp_df, color='black', jitter=0.3)
plt.suptitle('Functional expression in Structural Cells')
plt.subplots_adjust(top=0.85)
plt.ylim(-0.05, 1)
plt.savefig(os.path.join(SUPPLEMENTARY_FIG_DIR, 'ED_figure_4i.pdf'), bbox_inches='tight', dpi=300)

# look at M2 macrophages
temp_df = plot_df[plot_df.subset == 'all']
temp_df = temp_df[temp_df.metric == 'cluster_freq']
temp_df = temp_df[temp_df.cell_type == 'CD68_Mac']
temp_df = temp_df[temp_df.functional_marker.isin(['TIM3'])]

g = sns.catplot(x='ecm_cluster', y='value', data=temp_df,
                kind='box', col='functional_marker', sharey=False, color='grey', showfliers=False)
g.map_dataframe(sns.stripplot, x="ecm_cluster", y="value", data=temp_df, color='black', jitter=0.3)
plt.suptitle('Functional expression in CD68 Macrophages')
plt.subplots_adjust(top=0.85)
plt.ylim(-0.05, 1)
plt.savefig(os.path.join(SUPPLEMENTARY_FIG_DIR, 'ED_figure_4j.pdf'), bbox_inches='tight', dpi=300)

temp_df = plot_df[plot_df.subset == 'all']
temp_df = temp_df[temp_df.metric == 'cluster_freq']
temp_df = temp_df[temp_df.cell_type == 'CD68_Mac']
temp_df = temp_df[temp_df.functional_marker.isin(['PDL1'])]

g = sns.catplot(x='ecm_cluster', y='value', data=temp_df,
                kind='box', col='functional_marker', sharey=False, color='grey', showfliers=False)
g.map_dataframe(sns.stripplot, x="ecm_cluster", y="value", data=temp_df, color='black', jitter=0.3)
plt.suptitle('Functional expression in CD68 Macrophages')
plt.subplots_adjust(top=0.85)
plt.ylim(-0.05, 1)
plt.savefig(os.path.join(SUPPLEMENTARY_FIG_DIR, 'ED_figure_4k.pdf'), bbox_inches='tight', dpi=300)

# ECM cluster expression
channel_dir = os.path.join(BASE_DIR, 'image_data/samples')
analysis_dir = os.path.join(BASE_DIR, 'analysis_files')
ecm_mask_dir = os.path.join(BASE_DIR, 'intermediate_files/ecm/masks')

# generate crop sums
channels = ['Collagen1', 'Fibronectin', 'FAP']
crop_size = 256

fov_subset = ['TONIC_TMA3_R2C5', 'TONIC_TMA20_R5C3', 'TONIC_TMA9_R7C6', 'TONIC_TMA9_R3C1',
              'TONIC_TMA6_R11C6', 'TONIC_TMA20_R5C4', 'TONIC_TMA13_R1C5', 'TONIC_TMA10_R5C4',
              'TONIC_TMA24_R5C1', 'TONIC_TMA9_R11C1', 'TONIC_TMA23_R12C2', 'TONIC_TMA22_R9C1',
              'TONIC_TMA13_R11C6', 'TONIC_TMA17_R8C5', 'TONIC_TMA12_R4C3', 'TONIC_TMA13_R10C6',
              'TONIC_TMA19_R3C6', 'TONIC_TMA24_R4C3', 'TONIC_TMA21_R9C2', 'TONIC_TMA11_R6C4',
              'TONIC_TMA13_R5C4', 'TONIC_TMA7_R4C5', 'TONIC_TMA21_R1C4', 'TONIC_TMA20_R8C2',
              'TONIC_TMA2_R10C6', 'TONIC_TMA8_R7C6', 'TONIC_TMA20_R10C5', 'TONIC_TMA16_R10C6',
              'TONIC_TMA14_R8C2', 'TONIC_TMA23_R9C4', 'TONIC_TMA12_R10C5', 'TONIC_TMA4_R2C3',
              'TONIC_TMA11_R8C6', 'TONIC_TMA11_R2C1', 'TONIC_TMA15_R1C5', 'TONIC_TMA9_R9C6',
              'TONIC_TMA15_R2C5', 'TONIC_TMA14_R4C1', 'TONIC_TMA7_R8C5', 'TONIC_TMA9_R6C3',
              'TONIC_TMA14_R8C1', 'TONIC_TMA2_R12C4']

tiled_crops = generate_crop_sum_dfs(channel_dir=channel_dir,
                                          mask_dir=ecm_mask_dir,
                                          channels=channels,
                                          crop_size=crop_size, fovs=fov_subset, cell_table=None)


tiled_crops = normalize_by_ecm_area(crop_sums=tiled_crops, crop_size=crop_size,
                                          channels=channels)

# create a pipeline for normalization and clustering the data
kmeans_pipe = make_pipeline(preprocessing.PowerTransformer(method='yeo-johnson', standardize=True),
                            KMeans(n_clusters=2, random_state=0))

# select subset of data to train on
no_ecm_mask = tiled_crops.ecm_fraction < 0.1
train_data = tiled_crops[~no_ecm_mask]
train_data = train_data.loc[:, channels]

# fit the pipeline on the data
kmeans_pipe.fit(train_data.values)

kmeans_preds = kmeans_pipe.predict(tiled_crops[channels].values)

# get the transformed intermediate data
transformed_data = kmeans_pipe.named_steps['powertransformer'].transform(tiled_crops[channels].values)
transformed_df = pd.DataFrame(transformed_data, columns=channels)
transformed_df['tile_cluster'] = kmeans_preds
tiled_crops['tile_cluster'] = kmeans_preds
tiled_crops.loc[no_ecm_mask, 'tile_cluster'] = -1
transformed_df.loc[no_ecm_mask, 'tile_cluster'] = -1

# generate average image for each cluster
tile_replace_dict = {0: 'cold_collagen', 1: 'inflamed', -1: 'no_ecm'}
transformed_df['tile_cluster'] = transformed_df['tile_cluster'].replace(tile_replace_dict)
cluster_means = transformed_df[~no_ecm_mask].groupby('tile_cluster').mean()

# plot the average images
cluster_means_clustermap = sns.clustermap(cluster_means, cmap='Reds', figsize=(10, 10))
plt.savefig(os.path.join(SUPPLEMENTARY_FIG_DIR, 'ED_figure_4h.pdf'), dpi=300)
plt.close()
