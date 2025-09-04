import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

import matplotlib.pyplot as plt
import seaborn as sns

import numpy as np
import os
import pandas as pd

import python_files.supplementary_plot_helpers as supplementary_plot_helpers

BASE_DIR = "/Volumes/Shared/Noah Greenwald/TONIC_Cohort/"
SUPPLEMENTARY_FIG_DIR = os.path.join(BASE_DIR, "supplementary_figs")

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


## average cell counts by timepoint
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
    plt.savefig(os.path.join(cluster_stats_dir, "supp_figure_6{}.pdf".format(figure_name)), dpi=300)


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
