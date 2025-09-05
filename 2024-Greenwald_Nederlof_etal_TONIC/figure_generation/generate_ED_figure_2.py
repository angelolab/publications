import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

import matplotlib.pyplot as plt
import seaborn as sns

import os
import anndata
import pandas as pd
from matplotlib.ticker import FuncFormatter, MultipleLocator

import python_files.supplementary_plot_helpers as supplementary_plot_helpers

BASE_DIR = "/Volumes/Shared/Noah Greenwald/TONIC_Cohort/"
SUPPLEMENTARY_FIG_DIR = os.path.join(BASE_DIR, "supplementary_figs")
ANALYSIS_DIR = os.path.join(BASE_DIR, 'analysis_files')

def format_k(x, pos):
    return f'{int(x / 1000)}k'


# Plot cell counts per clustering level
cell_table = pd.read_csv('/Volumes/Shared/Noah Greenwald/TONIC_Cohort/analysis_files/cell_table_clusters.csv')
meta_data = pd.read_csv('/Volumes/Shared/Noah Greenwald/TONIC_Cohort/analysis_files/harmonized_metadata.csv')
meta_data = meta_data[['fov', 'Patient_ID', 'Timepoint', 'Tissue_ID']]

clusters = [
    ('cell_cluster_broad', '#F3B91A', 'Broad clusters', slice(4, 12), '6b.pdf'),
    ('cell_cluster', '#F6D7ED', 'Intermediate clusters', slice(4, 26), '6c.pdf'),
    ('cell_meta_cluster', '#94C47D', 'Detailed clusters', slice(4, 37), '6d.pdf')
]

for cluster_level, color, title, col_range, filename in clusters:
    # get cell population and whole image counts
    cell_table_sub = cell_table[['fov', cluster_level, 'label']]
    cell_counts = cell_table_sub.groupby(by=['fov', cluster_level]).count().reset_index()
    all_data = cell_counts

    # get mean for each timepoint
    combined_df = all_data.merge(meta_data, on=['fov'], how='left')
    combined_df = combined_df.rename(columns={'label': 'cell_counts'})

    combined_df = combined_df[[cluster_level, 'cell_counts', 'Tissue_ID']]
    tp_df = combined_df.groupby(by=['Tissue_ID', cluster_level]).mean().reset_index()
    tp_df = tp_df.drop(columns='Tissue_ID').groupby(by=cluster_level).sum().reset_index()

    t_sorted = tp_df.sort_values(by='cell_counts', ascending=True)
    plt.figure(figsize=(3, 5))
    plt.barh(t_sorted[cluster_level], t_sorted['cell_counts'], color=color)
    plt.xlabel('Cell count')
    plt.ylabel('')
    plt.title(title)
    plt.gca().xaxis.set_major_formatter(FuncFormatter(format_k))
    plt.gca().xaxis.set_major_locator(MultipleLocator(100000))
    plt.grid(visible=True, linestyle='--', alpha=0.7)
    plt.tight_layout()
    file_path = os.path.join(SUPPLEMENTARY_FIG_DIR, 'supp_figure_' + filename)
    plt.savefig(file_path, dpi=300, bbox_inches='tight')

# Average cell counts by timepoint
for cluster_level, figure_name in zip(['cell_cluster_broad', 'cell_cluster'], ['e', 'f']):
    cell_table = pd.read_csv('/Volumes/Shared/Noah Greenwald/TONIC_Cohort/analysis_files/cell_table_clusters.csv')

    # get cell population and whole image counts
    cell_table = cell_table[['fov', cluster_level, 'cell_meta_cluster']]
    cell_counts = cell_table.groupby(by=['fov', cluster_level]).count().reset_index()
    fov_cell_counts = cell_counts[['fov', 'cell_meta_cluster']].groupby(by=['fov']).sum().reset_index()
    fov_cell_counts[cluster_level] = 'Total cells'
    all_data = pd.concat([cell_counts, fov_cell_counts])

    # get mean for each timepoint
    combined_df = all_data.merge(meta_data, on=['fov'], how='left')
    combined_df = combined_df.rename(columns={'cell_meta_cluster': 'cell_counts'})
    combined_df = combined_df[[cluster_level, 'cell_counts', 'Timepoint']]
    avg_df = combined_df.groupby(by=['Timepoint', cluster_level]).mean().reset_index()

    # order and reformat timepoint names
    order = pd.DataFrame({'Timepoint': ['primary', 'baseline', 'pre_nivo', 'on_nivo'], 'Order': [1, 2, 3, 4]})
    avg_df = avg_df.merge(order, on=['Timepoint']).sort_values(by='Order')
    avg_df = avg_df.replace({'primary': 'Primary', 'baseline': 'Baseline', 'pre_nivo': 'Pre-nivo', 'on_nivo': 'On-nivo'})

    plt.figure().set_figwidth(20)
    plt.figure().set_figheight(6)

    if cluster_level == 'cell_cluster_broad':
        color_map = {'Total cells': 'black', 'Cancer': 'dimgrey', 'Structural': 'darksalmon',
                     'Mono_Mac': 'red', 'T': 'navajowhite', 'B': 'darkviolet', 'Other': 'yellowgreen',
                     'Granulocyte': 'aqua', 'NK': 'dodgerblue'}
        dims = (1.4, 0.8)
    else:
        color_map = {'Total cells': 'black', 'Cancer_1': 'grey',
                     'B': 'darkviolet', 'CD8T': 'tan', 'CD4T': 'navajowhite', 'CAF': 'chocolate',
                     'Cancer_2': 'darkgrey', 'Cancer_3': 'lightgrey', 'Monocyte': 'orangered',
                     'Fibroblast': 'peru', 'CD163_Mac': 'firebrick', 'Immune_Other': 'olivedrab',
                     'Endothelium': 'darksalmon', 'Treg': 'antiquewhite', 'APC': 'indianred',
                     'CD68_Mac': 'red', 'Smooth_Muscle': 'sienna',
                     'T_Other': 'wheat', 'Neutrophil': 'aqua', 'Mac_Other': 'darkred',
                     'Other': 'yellowgreen', 'NK': 'dodgerblue', 'Mast': 'turquoise'}
        dims = (1.02, 0.98)

    # plot each cell population
    for group in list(color_map.keys()):
        sub_df = avg_df[avg_df[cluster_level] == group]
        plt.plot('Timepoint', 'cell_counts', data=sub_df, marker='o', linewidth=2, label=group, color=color_map[group])

    plt.yscale('log')
    plt.legend(loc='center right', bbox_to_anchor=dims, fontsize=9)
    plt.ylabel('Log mean cell count', fontsize=12)
    plt.xlabel('Timepoint', fontsize=12)
    plt.title('Cell counts over time', fontsize=12)
    plt.tight_layout()
    sns.despine()
    plt.savefig(os.path.join(SUPPLEMENTARY_FIG_DIR, "ED_figure_2{}.pdf".format(figure_name)), dpi=300)

# Proportions of Cancer cells
cell_table = pd.read_csv('/Volumes/Shared/Noah Greenwald/TONIC_Cohort/analysis_files/cell_table_clusters.csv')
meta_data = pd.read_csv('/Volumes/Shared/Noah Greenwald/TONIC_Cohort/analysis_files/harmonized_metadata.csv')
meta_data = meta_data[['fov', 'Patient_ID', 'Timepoint', 'Tissue_ID']]

cell_table_sub = cell_table[['fov', 'cell_cluster', 'label']]
cell_counts = cell_table_sub.groupby(by=['fov', 'cell_cluster']).count().reset_index()
all_data = cell_counts[cell_counts.cell_cluster.isin(['Cancer_1', 'Cancer_2', 'Cancer_3'])]

# get patient counts
combined_df = all_data.merge(meta_data, on=['fov'], how='left')
combined_df = combined_df.rename(columns={'label': 'cell_counts'})
combined_df = combined_df[['cell_cluster', 'cell_counts', 'Patient_ID']]
patient_df = combined_df.groupby(by=['Patient_ID', 'cell_cluster']).sum().reset_index()
patient_df['cell_prop'] = patient_df['cell_counts'] / patient_df.groupby('Patient_ID')['cell_counts'].transform('sum')

# plot
colors = ['#F8766D', '#00BA38', '#619CFF']
plt.figure(figsize=(4, 5))
sns.violinplot(data=patient_df, x='cell_cluster', y='cell_prop',
               palette=colors, cut=0, scale='width', inner=None)
sns.stripplot(data=patient_df, x='cell_cluster', y='cell_prop',
              jitter=True, size=3, alpha=1, color='black', dodge=False)
plt.title('Baseline cancer cell types', size=12)
plt.xlabel('Cell Type')
plt.ylabel('% of total cancer cells')
plt.ylim((0,1))
plt.tight_layout()
plt.savefig(os.path.join(SUPPLEMENTARY_FIG_DIR, 'ED_figure_2g.pdf'), dpi=300, bbox_inches='tight')

# Proportions of Structural cells
all_data = cell_counts[cell_counts.cell_cluster.isin(['Endothelium', 'CAF', 'Fibroblast', 'Smooth_Muscle'])]
all_data['cell_cluster'] = all_data['cell_cluster'].replace({'Fibroblast': 'CAF-Other', 'CAF': 'CAF-S1'})

# get patient counts
combined_df = all_data.merge(meta_data, on=['fov'], how='left')
combined_df = combined_df.rename(columns={'label': 'cell_counts'})
combined_df = combined_df[['cell_cluster', 'cell_counts', 'Patient_ID']]
patient_df = combined_df.groupby(by=['Patient_ID', 'cell_cluster']).sum().reset_index()
patient_df['cell_prop'] = patient_df['cell_counts'] / patient_df.groupby('Patient_ID')['cell_counts'].transform('sum')
patient_df = patient_df.sort_values(by='cell_cluster')

# plot
colors = ['#F8766D', '#00BA38', '#619CFF', '#B39EB5']
plt.figure(figsize=(4, 5))
sns.violinplot(data=patient_df, x='cell_cluster', y='cell_prop',
               palette=colors, cut=0, scale='width', inner=None)
sns.stripplot(data=patient_df, x='cell_cluster', y='cell_prop',
              jitter=True, size=3, alpha=1, color='black', dodge=False)
plt.title('Baseline structural cell types', size=12)
plt.xlabel('Cell Type')
plt.ylabel('% of total structural cells')
plt.ylim((0,1))
plt.tight_layout()
plt.savefig(os.path.join(SUPPLEMENTARY_FIG_DIR, 'ED_figure_2h.pdf'), dpi=300, bbox_inches='tight')

# Proportions of Immune cells
all_data = cell_counts[cell_counts.cell_cluster.isin(['Treg', 'T_Other', 'NK', 'Neutrophil', 'Monocyte', 'Mast', 'Immune_Other', 'Mac_Other', 'CD163_Mac', 'CD68_Mac', 'CD8T', 'CD4T', 'B', 'APC'])]

# get patient counts
combined_df = all_data.merge(meta_data, on=['fov'], how='left')
combined_df = combined_df.rename(columns={'label': 'cell_counts'})
combined_df = combined_df[['cell_cluster', 'cell_counts', 'Patient_ID']]
patient_df = combined_df.groupby(by=['Patient_ID', 'cell_cluster']).sum().reset_index()
patient_df['cell_prop'] = patient_df['cell_counts'] / patient_df.groupby('Patient_ID')['cell_counts'].transform('sum')
lst = ['APC', 'B', 'CD4T', 'CD8T', 'CD68_Mac', 'CD163_Mac', 'Mac_Other', 'Immune_Other', 'Mast', 'Monocyte', 'Neutrophil', 'NK', 'T_Other', 'Treg']
patient_df = patient_df.set_index('cell_cluster')
patient_df = patient_df.loc[lst]
patient_df.reset_index(inplace=True)

# plot
colors = sns.color_palette("husl", 14)
plt.figure(figsize=(4, 5))
sns.violinplot(data=patient_df, x='cell_cluster', y='cell_prop',
               palette=colors, cut=0, scale='width', inner=None)
sns.stripplot(data=patient_df, x='cell_cluster', y='cell_prop',
              jitter=True, size=3, alpha=1, color='black', dodge=False)
plt.title('Baseline structural cell types', size=12)
plt.xlabel('Cell Type')
plt.ylabel('% of total structural cells')
plt.ylim((0, 1))
plt.tight_layout()
plt.savefig(os.path.join(SUPPLEMENTARY_FIG_DIR, 'ED_figure_2i.pdf'), dpi=300, bbox_inches='tight')

# PD1 expression
cell_table = pd.read_csv(os.path.join(BASE_DIR, 'analysis_files/combined_cell_table_normalized_cell_labels_updated.csv'))
mean_df = cell_table.groupby(['cell_cluster', 'fov'])['PD1'].mean().reset_index()
cd4 = mean_df[mean_df.loc[:, 'cell_cluster'] == 'CD4T'].melt(id_vars=['cell_cluster', 'fov'])
cd8 = mean_df[mean_df.loc[:, 'cell_cluster'] == 'CD8T'].melt(id_vars=['cell_cluster', 'fov'])
total = pd.concat([cd4, cd8], axis=0)
total['cell_cluster'] = list(total['cell_cluster'].values)

plt.figure(figsize=(4, 4))
sns.boxplot(x='cell_cluster', y='value', data=total.reset_index(), width=0.5, fliersize=0, color='lightgray')

# Plot stripplot
strip = sns.stripplot(x='cell_cluster', y='value', data=total.reset_index(), color='black',
                      size=4, alpha=0.7, jitter=True)

plt.xlabel('cell type', fontsize=12)
plt.ylabel('PD1 Expression', fontsize=12)
plt.tick_params(labelsize=10)
plt.ylim(0, 0.002)
plt.title('TONIC', fontsize=12)
plt.savefig(os.path.join(SUPPLEMENTARY_FIG_DIR, 'ED_figure_2j_tonic.pdf'), bbox_inches='tight')

# PD1 expression on CD4T as compared to CD8T in Wang et al.
adata = anndata.read_h5ad('/Volumes/Shared/Noah Greenwald/NTPublic/adata/adata_preprocessed.h5ad')
cell_table_nt = adata.to_df()
cell_table_nt['cell_cluster'] = adata.obs['cell_cluster']
cell_table_nt['fov'] = adata.obs['fov']

mean_df = cell_table_nt.groupby(['cell_cluster', 'fov'])['PD-1'].mean().reset_index()
cd4 = mean_df[mean_df.loc[:, 'cell_cluster'] == 'CD4T'].melt(id_vars = ['cell_cluster', 'fov'])
cd8 = mean_df[mean_df.loc[:, 'cell_cluster'] == 'CD8T'].melt(id_vars = ['cell_cluster', 'fov'])
total = pd.concat([cd4, cd8], axis=0)
total['cell_cluster'] = list(total['cell_cluster'].values)

plt.figure(figsize = (4,4))
sns.boxplot(x='cell_cluster', y='value', data=total.reset_index(), width=0.5, fliersize=0, color='lightgray')

# Plot stripplot
strip = sns.stripplot(x='cell_cluster', y='value', data=total.reset_index(),
                      color='black', size=4, alpha=0.7, jitter=True)

plt.xlabel('cell type', fontsize=12)
plt.ylabel('PD1 Expression', fontsize=12)
plt.tick_params(labelsize=10)
plt.title('NeoTRIP', fontsize=12)
plt.savefig(os.path.join(SUPPLEMENTARY_FIG_DIR, 'ED_figure_2j_nt.pdf'), bbox_inches='tight')


'''
# Cell identification and classification
# placeholder, right now these plots are in R

cell_table = pd.read_csv('/Volumes/Shared/Noah Greenwald/TONIC_Cohort/analysis_files/cell_table_clusters.csv')
cluster_order = {'Cancer': 0, 'Cancer_EMT': 1, 'Cancer_Other': 2, 'CD4T': 3, 'CD8T': 4, 'Treg': 5,
                 'T_Other': 6, 'B': 7, 'NK': 8, 'M1_Mac': 9, 'M2_Mac': 10, 'Mac_Other': 11,
                 'Monocyte': 12, 'APC': 13, 'Mast': 14, 'Neutrophil': 15, 'Fibroblast': 16,
                 'Stroma': 17, 'Endothelium': 18, 'Other': 19, 'Immune_Other': 20}
cell_table = cell_table.sort_values(by=['cell_cluster'], key=lambda x: x.map(cluster_order))

save_dir = '/Volumes/Shared/Noah Greenwald/TONIC_Cohort/supplementary_figs'

cluster_stats_dir = os.path.join(save_dir, "cluster_stats")
if not os.path.exists(cluster_stats_dir):
    os.makedirs(cluster_stats_dir)

## cell cluster counts
sns.histplot(data=cell_table, x="cell_cluster")
sns.despine()
plt.title("Cell Cluster Counts")
plt.xlabel("Cell Cluster")
plt.xticks(rotation=75)
plt.tight_layout()
plt.savefig(os.path.join(cluster_stats_dir, "cells_per_cluster.pdf"), dpi=300)


## cell type composition by tissue location of met and timepoint
meta_data = pd.read_csv('/Volumes/Shared/Noah Greenwald/TONIC_Cohort/analysis_files/harmonized_metadata.csv')
# meta_data = meta_data[['fov', 'Patient_ID', 'Timepoint', 'Localization']]
meta_data = meta_data[['fov', 'Patient_ID', 'Timepoint',]]

all_data = cell_table.merge(meta_data, on=['fov'], how='left')

# for metric in ['Localization', 'Timepoint']:
for metric in ['Timepoint']:
    data = all_data[all_data.Timepoint == 'baseline'] if metric == 'Localization' else all_data

    groups = np.unique(data.Localization) if metric == 'Localization' else \
        ['primary', 'baseline', 'pre_nivo', 'on_nivo']
    dfs = []
    for group in groups:
        sub_data = data[data[metric] == group]

        df = sub_data.groupby("cell_cluster_broad").count().reset_index()
        df = df.set_index('cell_cluster_broad').transpose()
        sub_df = df.iloc[:1].reset_index(drop=True)
        sub_df.insert(0, metric, [group])
        sub_df[metric] = sub_df[metric].map(str)
        sub_df = sub_df.set_index(metric)

        dfs.append(sub_df)
    prop_data = pd.concat(dfs).transform(func=lambda row: row / row.sum(), axis=1)

    color_map = {'Cancer': 'dimgrey', 'Stroma': 'darksalmon', 'T': 'navajowhite',
                 'Mono_Mac': 'red', 'B': 'darkviolet', 'Other': 'yellowgreen',
                 'Granulocyte': 'aqua', 'NK': 'dodgerblue'}

    means = prop_data.mean(axis=0).reset_index()
    means = means.sort_values(by=[0], ascending=False)
    prop_data = prop_data[means.cell_cluster_broad]

    colors = [color_map[cluster] for cluster in means.cell_cluster_broad]
    prop_data.plot(kind='bar', stacked=True, color=colors)
    sns.despine()
    plt.ticklabel_format(style='plain', useOffset=False, axis='y')
    plt.gca().set_ylabel("Cell Proportions")
    xlabel = "Tissue Location" if metric == 'Localization' else "Timepoint"
    plt.gca().set_xlabel(xlabel)
    plt.xticks(rotation=30)
    plt.title(f"Cell Type Composition by {xlabel}")
    handles, labels = plt.gca().get_legend_handles_labels()
    plt.legend(handles[::-1], labels[::-1],
               bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0)
    plt.tight_layout()
    plot_name = "cell_props_by_tissue_loc.pdf" if metric == 'Localization' else "cell_props_by_timepoint.pdf"
    plt.savefig(os.path.join(cluster_stats_dir, plot_name), dpi=300)


## colored cell cluster masks from random subset of 20 FOVs
random.seed(13)
seg_dir = '/Volumes/Shared/Noah Greenwald/TONIC_Cohort/segmentation_data/deepcell_output'

all_fovs = list(cell_table['fov'].unique())
fovs = random.sample(all_fovs, 20)
cell_table_subset = cell_table[cell_table.fov.isin(fovs)]

cohort_cluster_plot(
    fovs=fovs,
    seg_dir=seg_dir,
    save_dir=save_dir,
    cell_data=cell_table_subset,
    erode=True,
    fov_col='fov',
    label_col='label',
    cluster_col='cell_cluster_broad',
    seg_suffix="_whole_cell.tiff",
    cmap=color_map,
    display_fig=False,
)


# Functional marker thresholding
cell_table = pd.read_csv(
    os.path.join(ANALYSIS_DIR, "combined_cell_table_normalized_cell_labels_updated.csv")
)
functional_marker_viz_dir = os.path.join(SUPPLEMENTARY_FIG_DIR, "functional_marker_dist_thresholds_test")
if not os.path.exists(functional_marker_viz_dir):
    os.makedirs(functional_marker_viz_dir)


supplementary_plot_helpers.functional_marker_thresholding(
    cell_table, functional_marker_viz_dir, marker_info=marker_info,
    figsize=(20, 40)
)
'''
