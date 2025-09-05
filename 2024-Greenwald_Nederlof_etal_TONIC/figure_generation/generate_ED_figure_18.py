import os
import random
import anndata
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines

import matplotlib

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

BASE_DIR = "/Volumes/Shared/Noah Greenwald/TONIC_Cohort/"
SUPPLEMENTARY_FIG_DIR = os.path.join(BASE_DIR, "supplementary_figs")
ANALYSIS_DIR = os.path.join(BASE_DIR, "analysis_files")
INTERMEDIATE_DIR = os.path.join(BASE_DIR, "intermediate_files")

## 3.9 Location bias for features associated with immune cells ##
ranked_features_all = pd.read_csv(os.path.join(ANALYSIS_DIR, 'feature_ranking.csv'))
ranked_features = ranked_features_all.loc[
    ranked_features_all.comparison.isin(['primary', 'baseline', 'pre_nivo', 'on_nivo'])]
feature_metadata = pd.read_csv(os.path.join(ANALYSIS_DIR, 'feature_metadata.csv'))

immune_cells = ['Immune', 'Mono_Mac', 'B', 'T', 'Granulocyte', 'NK', 'CD68_Mac', 'CD163_Mac',
                'Mac_Other', 'Monocyte', 'APC', 'CD4T', 'CD8T', 'Treg', 'T_Other', 'Neutrophil', 'Mast']
immmune_feature_names = []
for cell in (immune_cells):
    for feat in ranked_features.feature_name.unique():
        if cell in feat:
            immmune_feature_names.append(feat)

immune_features = ranked_features[ranked_features.feature_name.isin(immmune_feature_names)]
immune_features_comp = immune_features[immune_features.compartment != 'all']

fig, axs = plt.subplots(2, 2, figsize=(8, 8), layout='constrained')
fig.suptitle('Changes in compartment features when dropping immune specific features')
fig2, axs2 = plt.subplots(2, 2, figsize=(8, 8), layout='constrained')
fig2.suptitle('Changes in compartment features when dropping immune specific features')

for immune_drop, coords in [[0, (0, 0)], [0.10, (0, 1)], [0.25, (1, 0)], [0.50, (1, 1)]]:
    i, j = coords
    feature_df = immune_features_comp

    # subset features
    idx_list = list(feature_df.index)
    sample_perc = int(len(idx_list) * immune_drop)
    sub_idx_list = random.sample(idx_list, sample_perc)
    sub_df = ranked_features[~ranked_features.index.isin(sub_idx_list)]
    sub_df_comp = ranked_features[~ranked_features.index.isin(sub_idx_list)]
    sub_df_comp = sub_df_comp[sub_df_comp.compartment != 'all']

    # calculate abundance of each compartment across all features
    drop_features = ranked_features[ranked_features.index.isin(sub_idx_list)].feature_name_unique
    feature_metadata_sub = feature_metadata[~feature_metadata.feature_name_unique.isin(drop_features)]
    feature_metadata_comp = feature_metadata_sub[feature_metadata_sub.compartment!='all']
    total_counts = feature_metadata_comp.groupby('compartment').count().iloc[:, 0]
    sub_prop = total_counts / np.sum(total_counts)

    # create df
    ratio_df = pd.DataFrame({'compartment': sub_prop.index, 'ratio': sub_prop.values})
    ratio_df = ratio_df.sort_values(by='ratio', ascending=False)
    ratio_df.loc[ratio_df.compartment == 'all', 'color_order'] = 5
    ratio_df.loc[ratio_df.compartment == 'cancer_core', 'color_order'] = 1
    ratio_df.loc[ratio_df.compartment == 'cancer_border', 'color_order'] = 2
    ratio_df.loc[ratio_df.compartment == 'stroma_border', 'color_order'] = 3
    ratio_df.loc[ratio_df.compartment == 'stroma_core', 'color_order'] = 4

    cmap = ['blue', 'deepskyblue', 'lightcoral', 'firebrick', 'grey']
    sns.barplot(data=ratio_df, x='compartment', y='ratio', hue='color_order', palette=cmap, ax=axs[i][j])
    sns.despine()
    axs[i][j].set_ylim(0, 0.5)
    axs[i][j].set_ylabel('')
    axs[i][j].tick_params(axis='x', labelrotation=60)
    axs[i][j].get_legend().remove()
    axs[i][j].set_title(f'Drop {int(immune_drop * 100)}%')
    if immune_drop == 0.0:
        axs[i][j].set_title(f'All features')

    # look at enrichment by compartment
    top_counts = sub_df.iloc[:100, :].groupby('compartment').count().iloc[:, 0]
    total_counts = sub_df.groupby('compartment').count().iloc[:, 0]

    # calculate abundance of each compartment in the top 100 and across all features
    top_prop = top_counts / np.sum(top_counts)
    total_prop = total_counts / np.sum(total_counts)
    top_ratio = top_prop / total_prop
    top_ratio = np.log2(top_ratio)

    # create df
    ratio_df = pd.DataFrame({'compartment': top_ratio.index, 'ratio': top_ratio.values})
    ratio_df = ratio_df.sort_values(by='ratio', ascending=False)

    ratio_df.loc[ratio_df.compartment == 'all', 'color_order'] = 1
    ratio_df.loc[ratio_df.compartment == 'cancer_core', 'color_order'] = 2
    ratio_df.loc[ratio_df.compartment == 'cancer_border', 'color_order'] = 3
    ratio_df.loc[ratio_df.compartment == 'stroma_border', 'color_order'] = 4
    ratio_df.loc[ratio_df.compartment == 'stroma_core', 'color_order'] = 5

    cmap = ['grey', 'blue', 'deepskyblue', 'lightcoral', 'firebrick']

    sns.barplot(data=ratio_df, x='compartment', y='ratio', hue='color_order', palette=cmap, ax=axs2[i][j])
    sns.despine()
    axs2[i][j].set_ylim(-0.7, 1.6)
    axs2[i][j].set_ylabel('')
    axs2[i][j].tick_params(axis='x', labelrotation=60)
    axs2[i][j].get_legend().remove()
    axs2[i][j].set_title(f'Drop {int(immune_drop * 100)}%')
    if immune_drop == 0.0:
        axs2[i][j].set_title(f'All features')

fig.savefig(os.path.join(SUPPLEMENTARY_FIG_DIR, 'supp_figure_18a.pdf'), bbox_inches='tight', dpi=300)
fig2.savefig(os.path.join(SUPPLEMENTARY_FIG_DIR, 'supp_figure_18b.pdf'), bbox_inches='tight', dpi=300)

## 4.5 Other/Stroma_Collagen/Stroma_Fibronectin/SMA/VIM to Cancer reassignment ##
cancer_recluster = pd.read_csv(os.path.join(BASE_DIR, 'supplementary_figs/review_figures/Cancer_reclustering', 'reassigned_cell_table.csv'))

fig, axes = plt.subplots(1, 3, figsize=(10, 4))
fig.suptitle('Cancer cell proportion of cell neighborhoods')
fig.supxlabel('Cancer proportion')
fig.supylabel('Density')
for c, ax in zip(['Other', 'Stroma_Collagen', 'Stroma_Fibronectin', 'SMA', 'VIM'], axes.flat):
    cluster_table_sub = cancer_recluster[cancer_recluster.cell_meta_cluster.isin([c, 'Cancer_new'])]
    ax.hist(cluster_table_sub.cancer_neighbors_prop, bins=20, density=True, rwidth=0.9, edgecolor='black')
    ax.set_xlim((0, 1))
    ax.set_ylim((0, 10))
    ax.set_title(f'{c} neighborhoods')
plt.savefig(os.path.join(SUPPLEMENTARY_FIG_DIR, 'supp_figure_18c.pdf'), bbox_inches='tight')


fig, ax = plt.subplots(1, 1, figsize=(8, 6))
cancer_table = cancer_recluster[cancer_recluster.cell_cluster_broad_new=='Cancer'][['fov', 'cell_cluster_new']]
cancer_table_counts = cancer_table.groupby(by=['cell_cluster_new'], observed=True).count().reset_index()
cancer_table_counts = cancer_table_counts.sort_values(by='fov', ascending=False)
total_count = cancer_table_counts[['fov']].sum()
cancer_table_counts['cancer_prop'] = cancer_table_counts['fov'].div(total_count.values[0])
plt.bar(cancer_table_counts.cell_cluster_new, cancer_table_counts.cancer_prop, color=['lightgrey']*(len(cancer_table_counts)-1)+['green'])
plt.xticks(rotation=45)
plt.title('Cancer cell cluster proportions')
green_patch = mpatches.Patch(color='green', label='New Cancer cells')
grey_patch = mpatches.Patch(color='lightgrey', label='Original Cancer cells')
plt.legend(handles=[green_patch, grey_patch])
plt.savefig(os.path.join(SUPPLEMENTARY_FIG_DIR, 'supp_figure_18e.pdf'), bbox_inches='tight')

counts_table = cancer_recluster[cancer_recluster.cell_meta_cluster.isin(['Other', 'Stroma_Collagen', 'Stroma_Fibronectin', 'SMA', 'VIM'])][['cell_meta_cluster', 'cell_meta_cluster_new']]
counts_table = counts_table.groupby(by=['cell_meta_cluster', 'cell_meta_cluster_new'], observed=True).value_counts().reset_index()
counts_table.loc[counts_table.cell_meta_cluster_new!='Cancer_new', 'cell_meta_cluster_new'] = 'Same'
counts_table = counts_table.pivot(index='cell_meta_cluster', columns='cell_meta_cluster_new', values='count').reset_index()
row_sums = counts_table.select_dtypes(include='number').sum(axis=1)
counts_table.iloc[:, 1:] = counts_table.iloc[:, 1:].div(row_sums, axis=0)

counts_table.plot(x='cell_meta_cluster', kind='barh', stacked=True, title='Proportion of reassigned Cancer cells',
                  mark_right=True, color=['green', 'lightgrey'])
green_patch = mpatches.Patch(color='green', label='New Cancer cells')
grey_patch = mpatches.Patch(color='lightgrey', label='Non-cancer')
plt.legend(handles=[green_patch, grey_patch])
plt.ylabel('Cell Type')
plt.savefig(os.path.join(SUPPLEMENTARY_FIG_DIR, 'supp_figure_18d.pdf'), bbox_inches='tight')

## 4.6.1 immune_agg features ##
immune_agg_viz_dir = os.path.join(BASE_DIR, 'supplementary_figs/review_figures', "immune_agg_features")

'''
immune_agg_analysis_dir = os.path.join(immune_agg_viz_dir, 'analysis_files')
os.makedirs(immune_agg_analysis_dir, exist_ok=True)

# group by timepoint
adata_processed = anndata.read_h5ad(os.path.join(ANALYSIS_DIR, 'adata_processed.h5ad'))
harmonized_metadata = pd.read_csv(os.path.join(ANALYSIS_DIR, 'harmonized_metadata.csv'))

fov_data_df = adata_processed.uns['combined_feature_data']
fov_data_df = pd.merge(fov_data_df, harmonized_metadata[['Tissue_ID', 'fov']], on='fov', how='left')
grouped = fov_data_df.groupby(['Tissue_ID', 'feature_name', 'feature_name_unique', 'compartment',
                               'cell_pop_level', 'feature_type']).agg({'raw_value': ['mean', 'std'],
                                                                       'normalized_value': ['mean', 'std']})
grouped.columns = ['raw_mean', 'raw_std', 'normalized_mean', 'normalized_std']
grouped = grouped.reset_index()
grouped.to_csv(os.path.join(immune_agg_analysis_dir, 'timepoint_features.csv'), index=False)

fov_data_df_filtered = adata_processed.uns['combined_feature_data_filtered']
fov_data_df_filtered = pd.merge(fov_data_df_filtered, harmonized_metadata[['Tissue_ID', 'fov']], on='fov', how='left')
grouped = fov_data_df_filtered.groupby(['Tissue_ID', 'feature_name', 'feature_name_unique', 'compartment',
                                 'cell_pop_level', 'feature_type']).agg({'raw_value': ['mean', 'std'],
                                                                            'normalized_value': ['mean', 'std']})

grouped.columns = ['raw_mean', 'raw_std', 'normalized_mean', 'normalized_std']
grouped = grouped.reset_index()
grouped.to_csv(os.path.join(immune_agg_analysis_dir, 'timepoint_features_filtered.csv'), index=False)


## 7_create_evolution_df.py converted
study_name = 'TONIC'

harmonized_metadata = pd.read_csv(os.path.join(ANALYSIS_DIR, 'harmonized_metadata.csv'))
timepoint_features = pd.read_csv(os.path.join(immune_agg_analysis_dir, 'timepoint_features_filtered.csv'))
timepoint_features_agg = timepoint_features.merge(
    harmonized_metadata[['Tissue_ID', 'Timepoint', 'Patient_ID'] + TIMEPOINT_COLUMNS].drop_duplicates(), on='Tissue_ID',
    how='left')
patient_metadata = pd.read_csv(os.path.join(INTERMEDIATE_DIR, f'metadata/{study_name}_data_per_patient.csv'))

# add evolution features to get finalized features specified by timepoint
combine_features(immune_agg_analysis_dir, harmonized_metadata, timepoint_features, timepoint_features_agg, patient_metadata, 
                timepoint_columns=TIMEPOINT_COLUMNS, drop_immune_agg=False)


## nivo_outcomes.py converted
harmonized_metadata = pd.read_csv(os.path.join(ANALYSIS_DIR, 'harmonized_metadata.csv'))
patient_metadata = pd.read_csv(os.path.join(INTERMEDIATE_DIR, 'metadata/TONIC_data_per_patient.csv'))
feature_metadata = adata_processed.uns['feature_metadata']

#
# To generate the feature rankings, you must have downloaded the patient outcome data.
#
outcome_data = pd.read_csv(os.path.join(INTERMEDIATE_DIR, 'metadata/patient_clinical_data.csv'))

# load previously computed results
combined_df = pd.read_csv(os.path.join(immune_agg_analysis_dir, 'timepoint_combined_features.csv'))
combined_df = combined_df.merge(outcome_data, on='Patient_ID')
combined_df = combined_df.loc[combined_df.Clinical_benefit.isin(['Yes', 'No']), :]
combined_df.to_csv(os.path.join(immune_agg_analysis_dir, 'timepoint_combined_features_outcome_labels.csv'), index=False)

# generate  pvalues and feature ranking
generate_feature_rankings(immune_agg_analysis_dir, combined_df, feature_metadata)

# preprocess feature sets for modeling
df_feature = pd.read_csv(os.path.join(immune_agg_analysis_dir, f'timepoint_combined_features_outcome_labels.csv'))
prediction_dir = os.path.join(immune_agg_viz_dir, 'prediction_model')
os.makedirs(prediction_dir, exist_ok=True)

df_feature.to_csv(os.path.join(prediction_dir, 'timepoint_combined_features_outcome_labels.csv'), index=False)

prediction_preprocessing(df_feature, prediction_dir)
os.makedirs(os.path.join(prediction_dir, 'patient_outcomes'), exist_ok=True)
'''

# generate violin plots for top immune_agg features
ranked_features_df = pd.read_csv(os.path.join(immune_agg_viz_dir, 'feature_ranking.csv'))
combined_df = pd.read_csv(os.path.join(immune_agg_viz_dir, 'timepoint_combined_features_outcome_labels.csv'))
feat_list = list(ranked_features_df[ranked_features_df.compartment == 'immune_agg'].feature_name_unique.drop_duplicates()[:11]) + ['immune_agg__proportion']
input_df_filtered = ranked_features_df[np.isin(ranked_features_df['feature_name_unique'], feat_list)].copy()
timepoints = ['primary', 'baseline', 'pre_nivo', 'on_nivo']

ranked_features_all = pd.read_csv(os.path.join(immune_agg_viz_dir, 'feature_ranking.csv'))
ranked_features = ranked_features_all.loc[ranked_features_all.comparison.isin(['primary', 'baseline', 'pre_nivo', 'on_nivo'])]

# make volcano plots
for subset in ['immune_agg']:
    plot_title = 'All Features'
    if subset == 'immune_agg':
        ranked_features = ranked_features[ranked_features.compartment == 'immune_agg']
        plot_title = 'Immune Aggregate Features'

    # plot total volcano
    fig, ax = plt.subplots(figsize=(3,3))
    sns.scatterplot(data=ranked_features, x='med_diff', y='log_pval', alpha=1, hue='importance_score',
                    palette=sns.color_palette("icefire", as_cmap=True), s=2.5, edgecolor='none', ax=ax)
    ax.set_xlim(-3, 3)
    ax.set_ylim(0, 8)
    sns.despine()

    # add gradient legend
    norm = plt.Normalize(ranked_features.importance_score.min(), ranked_features.importance_score.max())
    sm = plt.cm.ScalarMappable(cmap="icefire", norm=norm)
    ax.get_legend().remove()
    ax.figure.colorbar(sm, ax=ax)
    plt.title(plot_title)
    plt.tight_layout()
    plt.savefig(os.path.join(SUPPLEMENTARY_FIG_DIR, 'supp_figure_18h.pdf'))

# generate prediction comparison boxplot
preds = pd.read_csv(os.path.join(immune_agg_viz_dir, 'all_timepoints_results_MIBI-immune_agg.csv'))
preds = preds[['auc_primary_list', 'auc_baseline_list', 'auc_post_induction_list', 'auc_on_nivo_list']]
preds = preds.rename(columns={'auc_primary_list': 'Primary', 'auc_baseline_list': 'Baseline',
                              'auc_post_induction_list': 'Pre nivo', 'auc_on_nivo_list': 'On nivo'})
preds = preds.melt()
preds['immune_agg'] = 'include'

old_preds = pd.read_csv(os.path.join(immune_agg_viz_dir, 'all_timepoints_results_MIBI.csv'))
old_preds = old_preds[['auc_primary_list', 'auc_baseline_list', 'auc_post_induction_list', 'auc_on_nivo_list']]
old_preds = old_preds.rename(columns={'auc_primary_list': 'Primary', 'auc_baseline_list': 'Baseline',
                                      'auc_post_induction_list': 'Pre nivo', 'auc_on_nivo_list': 'On nivo'})
old_preds = old_preds.melt()
old_preds['immune_agg'] = 'exclude'
all_preds = pd.concat([preds, old_preds])

fig, ax = plt.subplots()
sns.boxplot(data=all_preds, x='variable', y='value', ax=ax, width=0.6, hue='immune_agg',
            palette=sns.color_palette(["#1f77b4", "darkseagreen"]))
sns.stripplot(data=all_preds, x='variable', y='value', ax=ax, hue='immune_agg',
              palette=sns.color_palette(["#1f77b4", "darkseagreen"]), dodge=True, jitter=0.2)
plt.xticks(rotation=45)
plt.title('TONIC prediction including and excluding immune aggregate features')
plt.ylabel('AUC')
plt.xlabel('')
# Add the custom legend
blue_line = mlines.Line2D([], [], color="#1f77b4", marker="o", label="including", linestyle='None')
green_line = mlines.Line2D([], [], color="darkseagreen", marker="o", label="excluding", linestyle='None')
plt.legend(handles=[blue_line, green_line], loc='lower right')
sns.despine()

plt.savefig(os.path.join(immune_agg_viz_dir, 'supp_figure_18i.pdf'), bbox_inches='tight', dpi=300)

# immune_agg properties
agg_count_df = pd.read_csv(os.path.join(immune_agg_viz_dir, 'immune_agg_counts_per_fov.csv'))
sns.histplot(agg_count_df.immune_agg_count, binwidth=1)
plt.title('Number of immune aggregates per image')
plt.xlabel('Immune Aggregate Objects')
plt.xlim((0, 23))
sns.despine()
plt.savefig(os.path.join(SUPPLEMENTARY_FIG_DIR, 'supp_figure_18f.pdf'))

mask_areas = pd.read_csv(os.path.join(INTERMEDIATE_DIR, 'mask_dir/fov_annotation_mask_area.csv'))
new_df = mask_areas[mask_areas.compartment == 'immune_agg'][['fov', 'area']]
new_df = new_df.rename(columns={'area': 'immune_agg_area'})
new_df = new_df.merge(mask_areas[mask_areas.compartment == 'all'][['fov', 'area']], on='fov')
new_df = new_df.rename(columns={'area': 'total_area'})
new_df['immune_agg_proportion'] = new_df['immune_agg_area'] / new_df['total_area']

sns.histplot(new_df.immune_agg_proportion)
plt.title('Proportion of immune aggregates per image')
plt.xlabel('Tissue Proportion')
sns.despine()
plt.savefig(os.path.join(SUPPLEMENTARY_FIG_DIR, 'supp_figure_18g.pdf'))

## 3.2 Low cellularity ##
low_cellularity_viz_dir = os.path.join(BASE_DIR, 'supplementary_figs/review figures', "low_cellularity")

# EXPLORATORY ANALYSIS
harmonized_metadata = pd.read_csv(ANALYSIS_DIR, 'harmonized_metadata.csv')
cellularity_df = pd.read_csv(os.path.join(low_cellularity_viz_dir, 'low_cellularity_images.csv'))

# low cellularity by patient
patient_counts = cellularity_df.groupby(by=['Patient_ID', 'low_cellularity']).count().reset_index()
patient_counts_total = cellularity_df[['fov', 'Patient_ID']].groupby(by='Patient_ID').count().reset_index()

bar1 = sns.barplot(x="Patient_ID",  y="fov", data=patient_counts_total, color='lightblue')
norm_cell = patient_counts[patient_counts.low_cellularity == 'No']
bar2 = sns.barplot(x="Patient_ID", y="fov", data=norm_cell, estimator=sum, ci=None,  color='darkblue')
top_bar = mpatches.Patch(color='lightblue', label='Low cellularity images')
bottom_bar = mpatches.Patch(color='darkblue', label='Regular images')
plt.legend(handles=[top_bar, bottom_bar])
plt.xticks([])
plt.title("Low Cellularity by Patient")
plt.ylabel("# of Images")
plt.savefig(os.path.join(low_cellularity_viz_dir, 'supp_figure_18k.pdf'), bbox_inches='tight', dpi=300)

# low cellularity by timepoint
cellularity_df = pd.read_csv(os.path.join(low_cellularity_viz_dir, 'low_cellularity_images.csv'))
cellularity_df = cellularity_df.merge(harmonized_metadata[['fov', 'Timepoint']], on='fov')
timepoint_counts = cellularity_df.groupby(by=['Timepoint', 'low_cellularity']).count().reset_index()
timepoint_counts_total = cellularity_df[['fov', 'Timepoint']].groupby(by='Timepoint').count().reset_index()
sort_dict = {'baseline': 0, 'primary': 1, 'pre_nivo': 2, 'on_nivo': 3}
timepoint_counts_total = timepoint_counts_total.iloc[timepoint_counts_total['Timepoint'].map(sort_dict).sort_values().index]

bar1 = sns.barplot(x="Timepoint",  y="fov", data=timepoint_counts_total, color='lightblue')
norm_cell = timepoint_counts[timepoint_counts.low_cellularity == 'No']
bar2 = sns.barplot(x="Timepoint", y="fov", data=norm_cell, estimator=sum, ci=None,  color='darkblue')
top_bar = mpatches.Patch(color='lightblue', label='Low cellularity images')
bottom_bar = mpatches.Patch(color='darkblue', label='Regular images')
plt.legend(handles=[top_bar, bottom_bar])
plt.xticks(rotation=45)
plt.title("Low Cellularity by Timepoint")
plt.ylabel("# of Images")
plt.savefig(os.path.join(SUPPLEMENTARY_FIG_DIR, 'supp_figure_18l.pdf'), bbox_inches='tight', dpi=300)

# low cellularity vs regular image features
ranked_features_df = pd.readcsv(os.path.join(low_cellularity_viz_dir, 'cellularity_feature_ranking.csv'))

# top features by type
top_fts = ranked_features_df[:100]
top_fts_type = top_fts[['feature_name_unique', 'feature_type']].groupby(by='feature_type').count().reset_index()
top_fts_type = top_fts_type.sort_values(by='feature_name_unique')
sns.barplot(top_fts_type, y='feature_type', x='feature_name_unique')
plt.ylabel("Feature Type")
plt.xlabel("Count")
plt.title("Top Features Differing Between Low Cellularity and Regular Images")
plt.savefig(os.path.join(low_cellularity_viz_dir, 'supp_figure_18m.pdf'), bbox_inches='tight', dpi=300)

# DROPPING LOW CELL IMAGES
adata = anndata.read_h5ad(os.path.join(ANALYSIS_DIR, 'adata_processed.h5ad'))
cell_table = adata.obs

cluster_counts = np.unique(cell_table.fov, return_counts=True)[1]
plt.figure(figsize=(8, 6))
g = sns.histplot(data=cluster_counts, bins=32)
y = cell_table.groupby(by='fov', observed=True).count().sort_values(by='label').label

plt.axvline(x=y.iloc[int(round(len(y)*.05, 0))], color='pink', linestyle='--')
plt.axvline(x=y.iloc[int(round(len(y)*.1, 0))], color='red', linestyle='--')
plt.axvline(x=y.iloc[int(round(len(y)*.15, 0))], color='purple', linestyle='--')
plt.axvline(x=y.iloc[int(round(len(y)*.2, 0))], color='black', linestyle='--')
sns.despine()
plt.title("Histogram of Cell Counts per Image")
plt.xlabel("Number of Cells in an Image")

# Create custom legend handles
line_5 = mlines.Line2D([], [], color="pink", marker='_', label="5% of images", linestyle='None')
line_10 = mlines.Line2D([], [], color="red", marker='_', label="10% of images", linestyle='None')
line_15 = mlines.Line2D([], [], color="purple", marker='_', label="15% of images", linestyle='None')
line_20 = mlines.Line2D([], [], color="black", marker='_', label="20% of images", linestyle='None')
plt.legend(handles=[line_5, line_10, line_15, line_20], loc='upper right')
plt.tight_layout()
plt.savefig(os.path.join(SUPPLEMENTARY_FIG_DIR, 'supp_figure_18j.pdf'), bbox_inches='tight', dpi=300)
