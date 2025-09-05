import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import anndata
import scanpy as sc
from sklearn.preprocessing import StandardScaler
import os
from matplotlib.gridspec import GridSpec

import matplotlib

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42


%reload_ext autoreload
%load_ext autoreload
%autoreload 2
%matplotlib inline

BASE_DIR = "/Volumes/Shared/Noah Greenwald/TONIC_Cohort/"
SUPPLEMENTARY_FIG_DIR = os.path.join(BASE_DIR, "supplementary_figs")
ANALYSIS_DIR = os.path.join(BASE_DIR, "analysis_files")

# TONIC DATA
markers = ["ECAD", "CK17", "CD45", "CD3", "CD4", "CD8", "FOXP3", "CD20", "CD56", "CD14", "CD68",
           "CD163", "CD11c", "HLADR", "ChyTr", "Calprotectin", "FAP",  "Fibronectin",
           "Collagen1", "Vim", "SMA", "CD31"]

cell_ordering = ['Cancer_1', 'Cancer_2', 'Cancer_3', 'CD4T', 'CD8T', 'Treg', 'T_Other', 'B',
                 'NK', 'CD68_Mac', 'CD163_Mac', 'Mac_Other', 'Monocyte', 'APC','Mast', 'Neutrophil',
                 'CAF', 'Fibroblast', 'Smooth_Muscle','Endothelium']
cols = ['Ki67','GLUT1','CD45RO','CD69', 'PD1','CD57','TBET', 'TCF1',
        'CD45RB', 'TIM3', 'PDL1','Fe','IDO','CD38']

cell_table = pd.read_csv(os.path.join(ANALYSIS_DIR, 'combined_cell_table_normalized_cell_labels_updated.csv'))

cell_table['cell_cluster_revised'] = cell_table['cell_cluster'].copy()
cell_table.loc[cell_table[cell_table['cell_meta_cluster'] == 'Cancer_Vim'].index, 'cell_cluster_revised'] = 'Cancer_Vim'
cell_table.loc[cell_table[cell_table['cell_meta_cluster'] == 'Cancer_SMA'].index, 'cell_cluster_revised'] = 'Cancer_SMA'
cell_table.loc[cell_table[cell_table['cell_meta_cluster'] == 'Cancer_Mono'].index, 'cell_cluster_revised'] = 'Cancer_HLADR'
cell_table.loc[cell_table[cell_table['cell_meta_cluster'] == 'Cancer_CK17'].index, 'cell_cluster_revised'] = 'Cancer_CK17'
cell_table.loc[cell_table[cell_table['cell_meta_cluster'] == 'Cancer_CD56'].index, 'cell_cluster_revised'] = 'Cancer_CD56'
cell_table.loc[cell_table[cell_table['cell_meta_cluster'] == 'Cancer_Ecad'].index, 'cell_cluster_revised'] = 'Cancer_ECAD'
cell_table.loc[cell_table[cell_table['cell_meta_cluster'] == 'Cancer_Other'].index, 'cell_cluster_revised'] = 'Cancer_Dim'

metadata_dir = os.path.join(BASE_DIR, 'intermediate_files/metadata')
harmonized_metadata = pd.read_csv(os.path.join(metadata_dir, 'harmonized_metadata.csv'))

marker_list = cell_table.columns[:84]
marker_list = [x for x in marker_list if '_nuclear' not in x]
cell_table_counts = cell_table.loc[:, ['fov', 'label', 'cell_meta_cluster', 'cell_cluster', 'cell_cluster_broad', 'cell_cluster_revised'] + marker_list]

markers = ["ECAD", "CK17", "CD45", "CD3", "CD4", "CD8", "FOXP3", "CD20", "CD56", "CD14", "CD68",
           "CD163", "CD11c", "HLADR", "ChyTr", "Calprotectin", "FAP",  "Fibronectin",
           "Collagen1", "Vim", "SMA", "CD31"]

cell_ordering = ['Cancer_ECAD', 'Cancer_CK17', 'Cancer_Vim', 'Cancer_SMA', 'Cancer_HLADR', 'Cancer_CD56', 'Cancer_Dim', 'CD4T', 'CD8T', 'Treg', 'T_Other', 'B',
                 'NK', 'CD68_Mac', 'CD163_Mac', 'Mac_Other', 'Monocyte', 'APC','Mast', 'Neutrophil',
                 'CAF', 'Fibroblast', 'Smooth_Muscle','Endothelium']

study_fovs = harmonized_metadata.loc[harmonized_metadata.Timepoint.isin(['primary', 'baseline', 'pre_nivo', 'on_nivo']), 'fov'].values
study_fovs = np.delete(study_fovs, np.where(study_fovs == 'TONIC_TMA14_R1C1'))
phenotype_col_name = "cell_cluster_revised"
cell_table_counts = cell_table_counts.loc[cell_table_counts.fov.isin(study_fovs), :]
mean_counts = cell_table_counts.groupby(phenotype_col_name)[markers].mean()
mean_counts = mean_counts.reindex(cell_ordering)

# set column order
mean_counts = mean_counts[markers]

plot_df = cell_table.loc[:, cols+['fov', 'cell_cluster_revised']].groupby(['fov', 'cell_cluster_revised']).mean().reset_index()
plot_df = plot_df.loc[:, cols + ['cell_cluster_revised']].groupby('cell_cluster_revised').mean()
plot_df = plot_df.reindex(cell_ordering)


cols = ['Ki67','GLUT1','CD45RO','CD69', 'PD1','CD57','TBET', 'TCF1',
        'CD45RB', 'TIM3', 'PDL1','Fe','IDO','CD38']
plot_df = plot_df[cols]

# combine together
combined_df = pd.concat([mean_counts, plot_df], axis=1)
X_ = combined_df.reset_index(drop = True)
X_ = pd.DataFrame(StandardScaler().fit_transform(X_), columns = combined_df.columns)
adata_plot = anndata.AnnData(X_)
adata_plot.obs['cell_cluster_revised'] = combined_df.index

phenotypic_markers = ['ECAD', 'CK17', 'CD45', 'CD3', 'CD4', 'CD8', 'FOXP3', 'CD20', 'CD56', 'CD14', 'CD68',
                    'CD163', 'CD11c', 'HLADR', 'ChyTr', 'Calprotectin', 'FAP', 'SMA', 'Vim', 'Fibronectin',
                    'Collagen1', 'CD31']


functional_markers = ['PDL1','Ki67','GLUT1','CD45RO','CD69', 'PD1','CD57','TBET', 'TCF1',
                        'CD45RB', 'TIM3','IDO', 'CD38']

fig = plt.figure(figsize=(6.75, 5.5), dpi = 400)
gs = GridSpec(1, 2, width_ratios=[1, 0.6], wspace=0, hspace=0.45, bottom=0.15)  #adjust this depending on how many phenotypic/functional markers you ahve
ax1 = fig.add_subplot(gs[0])
ax2 = fig.add_subplot(gs[1])
#phenotypic markers
mp1 = sc.pl.matrixplot(adata_plot,
                      var_names=phenotypic_markers,
                      groupby='cell_cluster_revised',
                      vmin=-3, vmax=3, cmap='vlag',
                      categories_order=cell_ordering,
                      ax=ax1,
                      colorbar_title='avg. expression \n (z-score)',
                      return_fig=True)

mp1.legend(show = False)
mp1.add_totals(size = 2, color = 'lightgrey', show = False).style(edge_color='black', cmap = 'vlag') #only have this here for edge colors, show = False so we don't see the totals on the first
ax1 = mp1.get_axes()['mainplot_ax']
y_ticks = ax1.get_yticks()

ax1.set_title('Phenotypic', fontsize = 14)

#functional markers
mp2 = sc.pl.matrixplot(adata_plot,
                      var_names=functional_markers,
                      groupby='cell_cluster_revised',
                      vmin=-3, vmax=3, cmap='vlag',
                      categories_order=cell_ordering,
                      ax=ax2,
                      colorbar_title='avg. expression \n (z-score)',
                      return_fig=True)
mp2.legend(show = False)
mp2.add_totals(size = 2, color = 'lightgrey', show = False).style(edge_color='black', cmap = 'vlag')
ax2 = mp2.get_axes()['mainplot_ax']
ax2.set_title('Functional', fontsize = 14)
ax2.set_yticklabels([])
ax2.set_yticks([])
plt.savefig(os.path.join(SUPPLEMENTARY_FIG_DIR, 'supp_figure_7a.pdf'), bbox_inches = 'tight')

feature_ranking = pd.read_csv(os.path.join('Cancer_reclustering', 'SpaceCat', 'feature_ranking.csv'))
feature_ranking_df = feature_ranking[np.isin(feature_ranking['comparison'], ['primary', 'baseline', 'pre_nivo', 'on_nivo'])]
feature_ranking_df['top_feature'] = False
feature_ranking_df['top_feature'][:100] = True

timepoint_features = pd.read_csv(os.path.join('Cancer_reclustering', 'SpaceCat', 'timepoint_combined_features_outcome_labels.csv'))
# feature_ranking_df['feature_name_unique'] = [i.replace('cell_cluster', 'cell_cluster_revised') for i in feature_ranking_df['feature_name_unique'] ]
feature_ranking_df['Timepoint'] = feature_ranking_df['comparison'].copy()
timepoint_features_old = pd.read_csv(os.path.join('SpaceCat_final', 'SpaceCat', 'timepoint_combined_features_outcome_labels.csv'))

merged_df = pd.merge(timepoint_features, feature_ranking_df, on = ['feature_name_unique', 'Timepoint'])
colors_dict = {'No': '#2089D5', 'Yes': 'lightgrey'}
feature_list = ['Ki67+__Cancer_HLADR', 'Ki67+__Cancer_Vim']
feature_list_nice = ['Cancer_HLADR', 'Cancer_Vim']
for i in range(0, len(feature_list)):
    feature = feature_list[i]
    f = feature_list_nice[i]
    fig, axes = plt.subplots(1, 1,
                             figsize=(4, 4),
                             gridspec_kw={'hspace': 0.3,
                                          'wspace': 0.3,
                                          'bottom': 0.15})

    df_subset = merged_df.iloc[
        np.where((merged_df['feature_name_unique'] == feature) & (merged_df['Timepoint'] == 'pre_nivo'))[0]].copy()
    g = sns.boxplot(y='raw_mean', x='Clinical_benefit', data=df_subset, palette=colors_dict, order=['No', 'Yes'],
                    ax=axes, linewidth=1, fliersize=0, width=0.6)
    g = sns.stripplot(y='raw_mean', x='Clinical_benefit', data=df_subset, linewidth=0.8, size=5, edgecolor="black",
                      jitter=True, ax=axes, palette=colors_dict, order=['No', 'Yes'])

    logq = df_subset['log10_qval'].values[0]
    pval = df_subset['pval'].values[0]
    med = df_subset['med_diff'].values[0]
    axes.margins(y=0.1)  # Adds 20% space on top (and bottom) of the tallest data point

    axes.text(
        0.75, 0.95,
        f"med_diff = {med:.3g}\np = {pval:.3g} \nlog10(FDR) = {logq:.3g}",
        transform=axes.transAxes,  # position is relative to the Axes (0 to 1)
        ha='center', va='top',
        fontsize=10,
        color='black'
    )

    g.tick_params(labelsize=10)
    g.set_xlabel('Clinical Benefit', fontsize=10)
    g.set_ylabel('Proportion', fontsize=10)
    axes.set_title(f'{f} Ki67+', fontsize=10)
    plt.ylim(-0.05, 1)
    plt.savefig(os.path.join(SUPPLEMENTARY_FIG_DIR, 'supp_figure_7c_Ki67+.pdf'), bbox_inches='tight')

colors_dict = {'No': '#2089D5', 'Yes': 'lightgrey'}
feature_list = ['Cancer_HLADR__cell_cluster_revised_density', 'Cancer_Vim__cell_cluster_revised_density']
feature_list_nice = ['Cancer_HLADR', 'Cancer_Vim']
for i in range(0, len(feature_list)):
    feature = feature_list[i]
    f = feature_list_nice[i]
    fig, axes = plt.subplots(1, 1,
                             figsize=(4, 4),
                             gridspec_kw={'hspace': 0.3,
                                          'wspace': 0.3,
                                          'bottom': 0.15})

    df_subset = merged_df.iloc[
        np.where((merged_df['feature_name_unique'] == feature) & (merged_df['Timepoint'] == 'pre_nivo'))[0]].copy()
    g = sns.boxplot(y='raw_mean', x='Clinical_benefit', data=df_subset, palette=colors_dict, order=['No', 'Yes'],
                    ax=axes, linewidth=1, fliersize=0, width=0.6)
    g = sns.stripplot(y='raw_mean', x='Clinical_benefit', data=df_subset, linewidth=0.8, size=5, edgecolor="black",
                      jitter=True, ax=axes, palette=colors_dict, order=['No', 'Yes'])

    logq = df_subset['log10_qval'].values[0]
    pval = df_subset['pval'].values[0]
    med = df_subset['med_diff'].values[0]
    axes.margins(y=0.1)  # Adds 20% space on top (and bottom) of the tallest data point

    axes.text(
        0.75, 0.95,
        f"med_diff = {med:.3g}\np = {pval:.3g} \nlog10(FDR) = {logq:.3g}",
        transform=axes.transAxes,  # position is relative to the Axes (0 to 1)
        ha='center', va='top',
        fontsize=10,
        color='black'
    )

    g.tick_params(labelsize=10)
    g.set_xlabel('Clinical Benefit', fontsize=10)
    g.set_ylabel('Mean', fontsize=10)
    axes.set_title(f'{f} density', fontsize=10)
    plt.savefig(os.path.join(SUPPLEMENTARY_FIG_DIR, 'supp_figure_7c_density.pdf'), bbox_inches='tight')

# NEOTRIP DATA
adata = anndata.read_h5ad('/Volumes/Shared/Noah Greenwald/NTPublic/adata/adata_preprocessed.h5ad')
cell_table = pd.read_csv('/Volumes/Shared/Noah Greenwald/NTPublic/data/derived_ark/final_cell_table.csv')
cell_table = cell_table[np.isin(cell_table['BiopsyPhase'], ['On-treatment', 'Baseline'])]
cell_table = cell_table[np.isin(cell_table['cellAnnotation'], ['TME', 'invasive'])]

markers = ['CK5/14', 'CK8/18', 'panCK', 'AR','CD45', 'CD3', 'CD4', 'CD8', 'FOXP3', 'CD20','CD79a', 'CD56', 'CD68', 'CD163', 'CD11c', 'HLA-DR',  'CD15', 'MPO', 'Calponin', 'SMA', 'Vimentin', 'PDGFRB','PDPN', 'CD31']
cols = ['PD-L1 (SP142)', 'PD-L1 (73-10)', 'IDO', 'PD-1', 'OX40', 'ICOS', 'CA9', 'c-PARP', 'Ki67', 'pH2AX', 'Helios', 'GATA3', 'T-bet', 'TCF1', 'TOX', 'GZMB', 'HLA-ABC']

adata_new = anndata.AnnData(pd.merge(cell_table.loc[:, ['fov', 'label'] + markers + cols], adata.obs, on = ['fov', 'label']).loc[:, markers+ cols])
adata_new.obs = adata.obs

cell_table_new = adata_new.to_df()
cell_table_new['cell_meta_cluster'] = adata_new.obs['cell_meta_cluster']
cell_table_new['fov'] = adata_new.obs['fov']
cell_table_new['label'] = adata_new.obs['label']

cell_ordering = ['CK^{hi}GATA3^{+}', 'MHC I&II^{hi}', 'PD-L1^{+}GZMB^{+}', 'PD-L1^{+}IDO^{+}', 'AR^{+}LAR',
                  'Apoptosis',  'CK8/18^{med}', 'panCK^{med}','CK^{lo}GATA3^{+}', 'TCF1^{+}', 'Helios^{+}', 'CA9^{+}Hypoxia',
                   'CD15^{+}', 'pH2AX^{+}DSB','Basal', 'CD56^{+}NE', 'Vimentin^{+}EMT',
                   'Treg','CD4^+PD1^+T','CD4^+TCF1^+T', 'CD8^+T','CD8^+PD1^+T_{Ex}', 'CD8^+TCF1^+T', 'CD8^+GZMB^+T',
                    'CD20^+B', 'CD79a^+Plasma', 'M2 Mac', 'DCs','CD56^+NK', 'CA9^+', 'PD-L1^+IDO^+APCs', 'PD-L1^+APCs','Neutrophils','Fibroblasts',
                    'Myofibroblasts', 'PDPN^+Stromal','Endothelial']

cell_ordering = [r'CK^{hi}GATA3^{+}', r'MHC I&II^{hi}', r'PD-L1^{+}GZMB^{+}', r'PD-L1^{+}IDO^{+}', r'AR^{+}LAR',
                  r'Apoptosis',  r'CK8/18^{med}', r'panCK^{med}',r'CK^{lo}GATA3^{+}', r'TCF1^{+}', r'Helios^{+}', r'CA9^{+}Hypoxia',
                   r'CD15^{+}', r'pH2AX^{+}DSB', r'Basal', r'CD56^{+}NE', r'Vimentin^{+}EMT',
                   r'Treg',r'CD4^+PD1^+T',r'CD4^+TCF1^+T', r'CD8^+T',r'CD8^+PD1^+T_{Ex}', r'CD8^+TCF1^+T', r'CD8^+GZMB^+T',
                    r'CD20^+B', r'CD79a^+Plasma', r'M2 Mac', r'DCs','CD56^+NK', r'CA9^+', r'PD-L1^+IDO^+APCs', r'PD-L1^+APCs',r'Neutrophils',r'Fibroblasts',
                    r'Myofibroblasts',r'PDPN^+Stromal',r'Endothelial']

original_labels = [
    'CK^{hi}GATA3^{+}', 'MHC I&II^{hi}', 'PD-L1^{+}GZMB^{+}', 'PD-L1^{+}IDO^{+}', 'AR^{+}LAR',
    'Apoptosis', 'CK8/18^{med}', 'panCK^{med}', 'CK^{lo}GATA3^{+}', 'TCF1^{+}', 'Helios^{+}',
    'CA9^{+}Hypoxia', 'CD15^{+}', 'pH2AX^{+}DSB', 'Basal', 'CD56^{+}NE', 'Vimentin^{+}EMT',
    'Treg', 'CD4^+PD1^+T', 'CD4^+TCF1^+T', 'CD8^+T', 'CD8^+PD1^+T_{Ex}', 'CD8^+TCF1^+T',
    'CD8^+GZMB^+T', 'CD20^+B', 'CD79a^+Plasma', 'M2 Mac', 'DCs', 'CD56^+NK', 'CA9^+',
    'PD-L1^+IDO^+APCs', 'PD-L1^+APCs', 'Neutrophils', 'Fibroblasts', 'Myofibroblasts',
    'PDPN^+Stromal', 'Endothelial'
]

unicode_mapping = {
    'CK^{hi}GATA3^{+}':   'CKhiGATA3⁺',
    'MHC I&II^{hi}':      'MHC I&IIhi',
    'PD-L1^{+}GZMB^{+}':   'PD-L1⁺GZMB⁺',
    'PD-L1^{+}IDO^{+}':    'PD-L1⁺IDO⁺',
    'AR^{+}LAR':          'AR⁺LAR',
    'Apoptosis':          'Apoptosis',
    'CK8/18^{med}':       'CK8/18med',
    'panCK^{med}':        'panCKmed',
    'CK^{lo}GATA3^{+}':   'CKloGATA3⁺',
    'TCF1^{+}':           'TCF1⁺',
    'Helios^{+}':         'Helios⁺',
    'CA9^{+}Hypoxia':     'CA9⁺Hypoxia',
    'CD15^{+}':           'CD15⁺',
    'pH2AX^{+}DSB':       'pH2AX⁺DSB',
    'Basal':              'Basal',
    'CD56^{+}NE':         'CD56⁺NE',
    'Vimentin^{+}EMT':     'Vimentin⁺EMT',
    'Treg':               'Treg',
    'CD4^+PD1^+T':        'CD4⁺PD1⁺T',
    'CD4^+TCF1^+T':       'CD4⁺TCF1⁺T',
    'CD8^+T':             'CD8⁺T',
    'CD8^+PD1^+T_{Ex}':    'CD8⁺PD1⁺Tₑₓ',
    'CD8^+TCF1^+T':       'CD8⁺TCF1⁺T',
    'CD8^+GZMB^+T':       'CD8⁺GZMB⁺T',
    'CD20^+B':           'CD20⁺B',
    'CD79a^+Plasma':      'CD79a⁺Plasma',
    'M2 Mac':             'M2 Mac',
    'DCs':                'DCs',
    'CD56^+NK':          'CD56⁺NK',
    'CA9^+':             'CA9⁺',
    'PD-L1^+IDO^+APCs':    'PD-L1⁺IDO⁺APCs',
    'PD-L1^+APCs':        'PD-L1⁺APCs',
    'Neutrophils':        'Neutrophils',
    'Fibroblasts':        'Fibroblasts',
    'Myofibroblasts':     'Myofibroblasts',
    'PDPN^+Stromal':      'PDPN⁺Stromal',
    'Endothelial':        'Endothelial'
}

cell_table_counts = cell_table_new.loc[:, ['fov', 'label', 'cell_meta_cluster'] + markers + cols]

phenotype_col_name = "cell_meta_cluster"
mean_counts = cell_table_counts.groupby(phenotype_col_name)[markers].mean()
mean_counts = mean_counts.reindex(cell_ordering)

# set column order
mean_counts = mean_counts[markers]

plot_df = cell_table_new.loc[:, cols+['fov', 'cell_meta_cluster']].groupby(['fov', 'cell_meta_cluster']).mean().reset_index()
plot_df = plot_df.loc[:, cols + ['cell_meta_cluster']].groupby('cell_meta_cluster').mean()
plot_df = plot_df.reindex(cell_ordering)

plot_df = plot_df[cols]

# combine together
combined_df = pd.concat([mean_counts, plot_df], axis=1)
X_ = combined_df.reset_index(drop = True)
X_ = pd.DataFrame(StandardScaler().fit_transform(X_), columns = combined_df.columns)
adata_plot = anndata.AnnData(X_)
adata_plot.obs['cell_meta_cluster'] = combined_df.index

fig = plt.figure(figsize=(7, 7.5), dpi = 400)
gs = GridSpec(1, 2, width_ratios=[1, 0.6], wspace=0, hspace=0.45, bottom=0.15)  #adjust this depending on how many phenotypic/functional markers you ahve
ax1 = fig.add_subplot(gs[0])
ax2 = fig.add_subplot(gs[1])
#phenotypic markers
mp1 = sc.pl.matrixplot(adata_plot,
                      var_names=markers,
                      groupby='cell_meta_cluster',
                      vmin=-3, vmax=3, cmap='vlag',
                      categories_order=cell_ordering,
                      ax=ax1,
                      colorbar_title='avg. expression \n (z-score)',
                      return_fig=True)

mp1.legend(show = False)
mp1.add_totals(size = 2, color = 'lightgrey', show = False).style(edge_color='black', cmap = 'vlag') #only have this here for edge colors, show = False so we don't see the totals on the first
ax1 = mp1.get_axes()['mainplot_ax']
y_ticks = ax1.get_yticks()
plt.draw()

ax1.set_title('Phenotypic', fontsize = 14)
current_labels = ax1.get_yticklabels()
new_labels = [ unicode_mapping.get(label.get_text(), label.get_text()) for label in current_labels]

ax1.set_yticklabels(new_labels)

#functional markers
mp2 = sc.pl.matrixplot(adata_plot,
                      var_names=cols,
                      groupby='cell_meta_cluster',
                      vmin=-3, vmax=3, cmap='vlag',
                      categories_order=cell_ordering,
                      ax=ax2,
                      colorbar_title='avg. expression \n (z-score)',
                      return_fig=True)
mp2.legend(show = False)
mp2.add_totals(size = 2, color = 'lightgrey', show = False).style(edge_color='black', cmap = 'vlag')
ax2 = mp2.get_axes()['mainplot_ax']
ax2.set_title('Functional', fontsize = 14)
ax2.set_yticklabels([])
ax2.set_yticks([])
plt.savefig(os.path.join(SUPPLEMENTARY_FIG_DIR, 'supp_figure_7b.pdf'), bbox_inches = 'tight')

cell_densities = pd.read_csv('cell_densities.csv')
ki67 = pd.read_csv('Ki67_positivity.csv')

cell_densities[cell_densities['BiopsyPhase'] == 'Baseline']
cell_densities = cell_densities[cell_densities['Arm'] == 'C&I']
cell_densities = cell_densities[np.isin(cell_densities['cellAnnotation'], ['TME', 'invasive'])]

ki67[ki67['BiopsyPhase'] == 'Baseline']
ki67 = ki67[ki67['Arm'] == 'C&I']

feature_list = ['MHC I&II^{hi}', 'Vimentin^{+}EMT']
feature_list_nice = ['MHC I&IIhi', 'Vimentin+EMT']
colors_dict = {'RD':'#2089D5','pCR':'lightgrey'}
for i in range(0, len(feature_list)):
    feature = feature_list[i]
    f = feature_list_nice[i]
    fig, axes = plt.subplots(1, 1,
                         figsize=(4, 4),
                         gridspec_kw={'hspace': 0.3,
                                      'wspace': 0.3,
                                      'bottom': 0.15})
    df_subset = ki67[(ki67['BiopsyPhase'] == 'Baseline') & (ki67['Label'] == feature)]
    # Boxplot + stripplot
    sns.boxplot(
        y='proportionKi67ByVar', x='pCR', data=df_subset, order = ['RD', 'pCR'],
        palette=colors_dict, ax=axes,
        linewidth=1, fliersize=0, width=0.6,
    )
    sns.stripplot(
        y='proportionKi67ByVar', x='pCR', data=df_subset, order = ['RD', 'pCR'],
        linewidth=0.8, size=5, edgecolor="black",
        jitter=True, ax=axes, palette=colors_dict
    )

    axes.set_title(f'{f} Ki67+', fontsize=10)
    axes.set_xlabel("Timepoint", fontsize = 10)  # optional: remove x-axis label
    axes.set_ylabel("Proportion", fontsize = 10)  # optional: label your y-axis
    plt.ylim(-0.05, 1)
    plt.savefig(os.path.join(SUPPLEMENTARY_FIG_DIR, 'supp_figure_7d_Ki67+.pdf'), bbox_inches = 'tight')

feature_list = ['MHC I&II^{hi}', 'Vimentin^{+}EMT']
feature_list_nice = ['MHC I&IIhi', 'Vimentin+EMT']
colors_dict = {'RD':'#2089D5','pCR':'lightgrey'}
for i in range(0, len(feature_list)):
    feature = feature_list[i]
    f = feature_list_nice[i]
    fig, axes = plt.subplots(1, 1,
                         figsize=(4, 4),
                         gridspec_kw={'hspace': 0.3,
                                      'wspace': 0.3,
                                      'bottom': 0.15})
    df_subset = cell_densities[(cell_densities['BiopsyPhase'] == 'Baseline') & (cell_densities['Label'] == feature)]
    # Boxplot + stripplot
    sns.boxplot(
        y='Density', x='pCR', data=df_subset, order = ['RD', 'pCR'],
        palette=colors_dict, ax=axes,
        linewidth=1, fliersize=0, width=0.6,
    )
    sns.stripplot(
        y='Density', x='pCR', data=df_subset, order = ['RD', 'pCR'],
        linewidth=0.8, size=5, edgecolor="black",
        jitter=True, ax=axes, palette=colors_dict
    )

    axes.set_title(f'{f} Density', fontsize=10)
    axes.set_xlabel("Timepoint", fontsize = 10)  # optional: remove x-axis label
    axes.set_ylabel("raw_mean", fontsize = 10)  # optional: label your y-axis
    plt.savefig(os.path.join(SUPPLEMENTARY_FIG_DIR, 'supp_figure_7d_density.pdf'), bbox_inches = 'tight')