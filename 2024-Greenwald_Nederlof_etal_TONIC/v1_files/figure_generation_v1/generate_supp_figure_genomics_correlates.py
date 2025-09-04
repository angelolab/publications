import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from statsmodels.stats.multitest import multipletests
from scipy.stats import spearmanr


BASE_DIR = "/Volumes/Shared/Noah Greenwald/TONIC_Cohort/"
raw_dir = "/Volumes/Shared/Noah Greenwald/TONIC_Acquisition/"
SUPPLEMENTARY_FIG_DIR = os.path.join(BASE_DIR, "supplementary_figs")
sequence_dir = os.path.join(BASE_DIR, 'sequencing_data')

# transcript-level correlation
corr_df = pd.read_csv(os.path.join(sequence_dir, 'analysis/MIBI_RNA_correlation.csv'))

# annotate correlations based on marker
annot_dict = {'T cell': ['CD3', 'CD4', 'CD8', 'FOXP3'],
              'B cell': ['CD20'],
              'NK cell': ['CD56'],
              'Granulocyte': ['Calprotectin', 'ChyTr'],
              'Monocyte': ['CD14', 'CD68', 'CD163', 'CD11c'],
              'Stroma': ['Vim', 'CD31', 'SMA'],
              'ECM': ['FAP', 'Collagen1', 'Fibronectin'],
              'Checkpoints': ['PDL1', 'PD1', 'IDO', 'TIM3'],
              'Tumor': ['ECAD', 'CK17'],
              'Functional': ['CD38', 'HLA1', 'HLADR', 'Ki67', 'GLUT1', 'CD45RO', 'CD45RB',
                             'CD57', 'TCF1', 'TBET', 'CD69'],
              'Other': ['H3K9ac', 'H3K27me3', 'CD45']}

for key, value in annot_dict.items():
    corr_df.loc[corr_df.marker.isin(value), 'group'] = key

corr_df['single_group'] = 'Single'
corr_df = corr_df.sort_values('correlation', ascending=False)

# annotated box plot
fig, ax = plt.subplots(1, 1, figsize=(2, 4))
sns.stripplot(data=corr_df, x='single_group', y='correlation', color='black', ax=ax)
sns.boxplot(data=corr_df, x='single_group', y='correlation', color='grey', ax=ax, showfliers=False,
            width=0.3)

ax.set_ylabel('RNA vs image marker correlation')
plt.xticks(rotation=45)
plt.tight_layout()
fig.savefig(os.path.join(SUPPLEMENTARY_FIG_DIR, 'supp_figure_9a.pdf'), dpi=300)
plt.close()


# plot lineage agreement from RNA data
rna_correlations = pd.read_csv(os.path.join(sequence_dir, 'genomics_image_correlation_RNA.csv'))
rna_correlations = rna_correlations.loc[~(rna_correlations.genomic_features.apply(lambda x: 'rna' in x)), :]


populations, rna_feature, image_feature, values = [], [], [], []

pairings = {'T cells': [['T_cells', 'T_cell_traffic'], ['T__cluster_broad_density']],
            'CD8T cells': [['Effector_cells'], ['CD8T__cluster_density']],
            'APC': [['MHCII', 'Macrophage_DC_traffic'], ['APC__cluster_density']],
            'B cells': [['B_cells'], ['B__cluster_broad_density']],
            'NK cells': [['NK_cells'], ['NK__cluster_broad_density']],
            'T regs': [['T_reg_traffic', 'Treg'], ['Treg__cluster_density']],
            'Endothelium': [['Endothelium'], ['Endothelium__cluster_density']],
            'Macrophages': [['Macrophages', 'Macrophage_DC_traffic', 'M1_signatures'],
                            ['M1_Mac__cluster_density', 'M2_Mac__cluster_density', 'Monocyte__cluster_density', 'Mac_Other__cluster_density']],
            'Granulocytes': [['Granulocyte_traffic', 'Neutrophil_signature'], ['Neutrophil__cluster_density', 'Mast__cluster_density']],
            'ECM': [['CAF', 'Matrix_remodeling', 'Matrix'], ['Fibroblast__cluster_density', 'Stroma__cluster_density']]}

# look through  pairings and pull out correlations
for pop, pairs in pairings.items():
    for pair_1 in pairs[0]:
        for pair_2 in pairs[1]:
            pop_df = rna_correlations.loc[(rna_correlations.image_feature == pair_2) &
                                          (rna_correlations.genomic_features == pair_1), :]
            if len(pop_df) == 0:
                print('No data for: ', pop, pair_1, pair_2)
                continue
            populations.append(pop)
            rna_feature.append(pair_1)
            image_feature.append(pair_2)
            values.append(pop_df.cor.values[0])

pop_correlations = pd.DataFrame({'populations': populations, 'correlation': values, 'rna_feature': rna_feature,
                                 'image_feature': image_feature})

# plot results
fig, ax = plt.subplots(1, 1, figsize=(7, 5))
sns.stripplot(data=pop_correlations, x='populations', y='correlation', dodge=True, ax=ax)
ax.set_ylabel('RNA vs image lineage correlation')
plt.xticks(rotation=45)
plt.tight_layout()
fig.savefig(os.path.join(SUPPLEMENTARY_FIG_DIR, 'supp_figure_9b.pdf'), dpi=300)
plt.close()


# Same thing for DNA data
dna_correlations = pd.read_csv(os.path.join(sequence_dir, 'genomics_image_correlation_DNA.csv'))
dna_correlations['log_qval'] = -np.log10(dna_correlations.fdr)

lineage_features = [pairings[x][1] for x in pairings.keys()]
lineage_features = [item for sublist in lineage_features for item in sublist]

SNVs = [x for x in dna_correlations.genomic_features.unique() if 'SNV' in x]
amps = [x for x in dna_correlations.genomic_features.unique() if '_cn' in x]
alterations = SNVs + amps

dna_correlations = dna_correlations.loc[dna_correlations.image_feature.isin(lineage_features), :]
dna_correlations = dna_correlations.loc[dna_correlations.genomic_features.isin(alterations), :]

# update fdr calculation
dna_correlations['qval_subset'] = multipletests(dna_correlations.pval.values, method='fdr_bh')[1]
dna_correlations['log_qval_subset'] = -np.log10(dna_correlations.qval_subset)

sns.scatterplot(data=dna_correlations, x='cor', y='log_qval_subset', edgecolor='none', s=7.5)
plt.xlabel('Spearman Correlation')
plt.ylabel('-log10(q-value)')
plt.title('Correlation between DNA and Image Features')
plt.ylim(0, 0.2)
sns.despine()
plt.savefig(os.path.join(SUPPLEMENTARY_FIG_DIR, 'supp_figure_9d.pdf'), dpi=300)
plt.close()

# compare sTILs with MIBI densities
patient_metadata = pd.read_csv(os.path.join(BASE_DIR, 'intermediate_files/metadata/TONIC_data_per_patient.csv'))
image_feature_df = pd.read_csv(os.path.join(BASE_DIR, 'analysis_files/timepoint_combined_features.csv'))

tils = patient_metadata[['Patient_ID', 'sTIL_(%)_revised']].drop_duplicates()
mibi_tils = image_feature_df.loc[image_feature_df.Timepoint == 'baseline', :]
mibi_tils = mibi_tils.loc[mibi_tils.feature_name_unique.isin(['T__cluster_broad_density', 'B__cluster_broad_density']), :]
mibi_tils = mibi_tils[['Patient_ID', 'raw_mean']]

mibi_tils = mibi_tils.groupby('Patient_ID').sum().reset_index()
mibi_tils = mibi_tils.rename(columns={'raw_mean': 'MIBI_density'})

combined_tils = pd.merge(tils, mibi_tils, on='Patient_ID', how='inner')
combined_tils = combined_tils.dropna(subset=['sTIL_(%)_revised', 'MIBI_density'])

# plot
sns.scatterplot(data=combined_tils, x='sTIL_(%)_revised', y='MIBI_density')

plt.xlabel('sTIL (%)')
plt.ylabel('MIBI density')
plt.title('sTIL vs MIBI density')
plt.savefig(os.path.join(SUPPLEMENTARY_FIG_DIR, 'supp_figure_9c.pdf'), dpi=300)
plt.close()

spearmanr(combined_tils['sTIL_(%)_revised'], combined_tils['MIBI_density'])
