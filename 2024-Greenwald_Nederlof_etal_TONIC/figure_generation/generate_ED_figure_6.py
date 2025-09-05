import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import gseapy as gp
from gseapy import Msigdb
import os
import numpy as np
import scipy.stats as sp

import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

BASE_DIR = "/Volumes/Shared/Noah Greenwald/TONIC_Cohort/"
SUPPLEMENTARY_FIG_DIR = os.path.join(BASE_DIR, "supplementary_figs")
ANALYSIS_DIR = os.path.join(BASE_DIR, "analysis_files")


gene_revisions_dir = os.path.join(BASE_DIR, "sequencing_data/preprocessing/")
'''
# Load gene expression data
gene_features = pd.read_csv(os.path.join(gene_revisions_dir, "TONIC_gene_expression_TPM_table.tsv"), sep='\t', index_col=0)

# Get genesets
msig = Msigdb()
categories_keep = ["h.all", "c2.all", "c5.go"]
all_gene_sets = {}
for one_category in categories_keep:
    one_set = msig.get_gmt(category=one_category, dbver="2023.2.Hs")
    all_gene_sets.update(one_set)

# Keep specific gene sets
gene_sets_keep = pd.read_csv("msigdb_genesets.csv")
gene_sets_names_keep = gene_sets_keep['name'].values
gene_sets = {x: all_gene_sets[x] for x in gene_sets_names_keep}

# Do ssgsea
ss = gp.ssgsea(data = gene_features,
               gene_sets = gene_sets,
               sample_norm_method = 'rank')
ss_out_tab = ss.res2d
ss_out_tab.head()
ss_out_tab.to_csv(os.path.join(gene_revisions_dir, "msigdb_ssgsea_es.csv"), index=False)

save_dat = ss_out_tab.pivot(
    index='Name',
    columns='Term',
    values='NES'
).reset_index()
save_dat = save_dat.rename(columns={'Name':'rna_seq_sample_id'})
save_dat.to_csv(os.path.join(gene_revisions_dir, "msigdb_ssgsea_nes_pivot.csv"), index=False)
'''

results = pd.read_csv(os.path.join(gene_revisions_dir, "msigdb_ssgsea_es.csv"))
results = results.sort_values('NES', ascending=False)
medians = results.groupby('Term')['NES'].median().sort_values(ascending=False)

plt.figure(figsize=(14,5))
sns.boxplot(data=results, x='Term', y='NES', order=medians.index, showfliers=False)
sns.stripplot(data=results, x='Term', y='NES', order=medians.index, jitter=True, color='gray', alpha=0.5, size=2)
plt.xticks(rotation=90, ha='center')
plt.savefig(os.path.join(SUPPLEMENTARY_FIG_DIR, 'ED_figure_6a.pdf'), dpi=300, bbox_inches='tight')

## Get correlations between MIBI and RNAseq data
# Load MIBI data
mibi_features = pd.read_csv("/mnt/nas_share/Noah Greenwald/TONIC_Cohort/analysis_files/timepoint_combined_features.csv")
mibi_features_names = mibi_features['feature_name_unique'].values

# Only keep baseline timepoints
mibi_features = mibi_features[mibi_features['Timepoint']=='baseline']
mibi_dat = mibi_features.pivot(index='Patient_ID',columns='feature_name_unique',values='normalized_mean')

# Load metadata
metadata = pd.read_csv("/mnt/nas_share/Noah Greenwald/TONIC_Cohort/analysis_files/harmonized_metadata.csv")
metadata_sub = metadata[["Patient_ID","rna_seq_sample_id","Timepoint"]]
metadata_sub = metadata_sub.dropna()
metadata_sub = metadata_sub.drop_duplicates()
metadata_sub = metadata_sub[metadata_sub['Timepoint']=='baseline']

# Load ssgsea scores
ssgsea = pd.read_csv("msigdb_ssgsea_nes_pivot.csv")

# Load literture gene set scores
lit_gene_sets = pd.read_csv("literature_gene_set_scores.csv")

# Load immune sigs/CIBERSORT data
cibersort = pd.read_csv("TONIC_immune_sig_score_table_addCBX_250121.tsv", sep='\t')
cibersort = cibersort.rename(columns={'sample_identifier':'rna_seq_sample_id'})

# Combine all RNAseq scores
all_rnaseq = pd.merge(ssgsea, lit_gene_sets, on='rna_seq_sample_id', how='inner')
all_rnaseq = pd.merge(all_rnaseq, cibersort, on='rna_seq_sample_id', how='inner')

# Only keep baseline timepoints
rnaseq_sub = all_rnaseq.merge(metadata_sub[['Patient_ID','rna_seq_sample_id']], on='rna_seq_sample_id', how='inner')
rnaseq_sub = rnaseq_sub.set_index('Patient_ID')
rnaseq_sub = rnaseq_sub.drop('rna_seq_sample_id', axis=1)

# Only keep samples that have both measurements
common_indices = mibi_dat.index.intersection(rnaseq_sub.index)
mibi_dat_filtered = mibi_dat.loc[common_indices]
rnaseq_dat_filtered = rnaseq_sub.loc[common_indices]

# Only keep samples that have both measurements
common_indices = mibi_dat.index.intersection(rnaseq_sub.index)
mibi_dat_filtered = mibi_dat.loc[common_indices]
rnaseq_dat_filtered = rnaseq_sub.loc[common_indices]

# Get correlations
corr_matrix = pd.DataFrame(index=rnaseq_dat_filtered.columns, columns=mibi_dat_filtered.columns)
pval_matrix = pd.DataFrame(index=rnaseq_dat_filtered.columns, columns=mibi_dat_filtered.columns)

# Calculate correlations for each pair
for img_feature in mibi_dat_filtered.columns:
    for rnaseq_feature in rnaseq_dat_filtered.columns:
        # Get complete cases for this pair
        mask = ~(rnaseq_dat_filtered[rnaseq_feature].isna() | mibi_dat_filtered[img_feature].isna())
        if mask.sum() > 9:  # Filter for at least 10 points
            corr, pval = sp.spearmanr(rnaseq_dat_filtered.loc[mask, rnaseq_feature],
                                      mibi_dat_filtered.loc[mask, img_feature],
                                      nan_policy='omit')
            corr_matrix.loc[rnaseq_feature, img_feature] = corr
            pval_matrix.loc[rnaseq_feature, img_feature] = pval
        else:
            corr_matrix.loc[rnaseq_feature, img_feature] = np.nan
            pval_matrix.loc[rnaseq_feature, img_feature] = np.nan

corr_df = corr_matrix.reset_index().melt(
    id_vars=['index'],
    var_name='mibi_feature',
    value_name='corr'
)
corr_df = corr_df.rename(columns={'index':'rnaseq_feature'})

pval_df = pval_matrix.reset_index().melt(
    id_vars=['index'],
    var_name='mibi_feature',
    value_name='pval'
)
pval_df = pval_df.rename(columns={'index':'rnaseq_feature'})

corr_pval_df = pd.merge(corr_df, pval_df, on=['rnaseq_feature', 'mibi_feature'])
corr_pval_df = corr_pval_df.dropna()
corr_pval_df['corr'] = pd.to_numeric(corr_pval_df['corr'])
corr_pval_df['pval'] = pd.to_numeric(corr_pval_df['pval'])
corr_pval_df['log10_pval'] = -np.log10(corr_pval_df['pval'])

corr_pval_df['corr_nosign'] = abs(corr_pval_df['corr'])
corr_pval_df['rank_val'] = corr_pval_df['corr_nosign'].rank()

fig, ax = plt.subplots(figsize=(5,4))
# color pallete options: Greys, magma, vlag, icefire
sns.scatterplot(data=corr_pval_df, x='corr', y='log10_pval', alpha=1, hue='rank_val', palette=sns.color_palette("icefire", as_cmap=True),
                s=2.5, edgecolor='none', ax=ax)
ax.set_xlim(-1, 1)
sns.despine()

# add gradient legend
norm = plt.Normalize(corr_pval_df.rank_val.min(), corr_pval_df.rank_val.max())
sm = plt.cm.ScalarMappable(cmap="icefire", norm=norm)
ax.get_legend().remove()
ax.figure.colorbar(sm, ax=ax)
plt.savefig(os.path.join(SUPPLEMENTARY_FIG_DIR, 'ED_figure_6b.pdf'))
plt.tight_layout()
plt.show()

# enrichment by compartment
mibi_feature_metadata = pd.read_csv("/mnt/nas_share/Noah Greenwald/TONIC_Cohort/analysis_files/feature_metadata.csv")
mibi_feature_metadata = mibi_feature_metadata.rename(columns={'feature_name_unique': 'mibi_feature'})

spatial_features = ['mixing_score', 'cell_diversity', 'compartment_area_ratio', 'pixie_ecm',
                    'compartment_area', 'fiber', 'linear_distance', 'ecm_fraction', 'ecm_cluster']
spatial_mask = np.logical_or(mibi_feature_metadata.feature_type.isin(spatial_features), mibi_feature_metadata.compartment != 'all')
mibi_feature_metadata['spatial_feature'] = spatial_mask

corr_cutoff = 0.6
pval_cutoff = 0.05
corr_pval_df_with_features = pd.merge(corr_pval_df, mibi_feature_metadata[['mibi_feature','compartment','feature_type','spatial_feature']], on='mibi_feature')
pass_threshold_mask = (corr_pval_df_with_features['corr_nosign'] > corr_cutoff) & (corr_pval_df_with_features['pval'] < pval_cutoff)
corr_pval_df_with_features['pass_threshold'] = pass_threshold_mask
good_corrs = corr_pval_df_with_features[pass_threshold_mask]

# look at enrichment by compartment
top_counts = good_corrs.groupby('compartment').size()
total_counts = corr_pval_df_with_features.groupby('compartment').size()

# calculate abundance of each compartment in the thresholded features and across all features
top_prop = top_counts / np.sum(top_counts)
total_prop = total_counts / np.sum(total_counts)
top_ratio = top_prop / total_prop
top_ratio = np.log2(top_ratio)

# create df
ratio_df = pd.DataFrame({'compartment': top_ratio.index, 'ratio': top_ratio.values})
ratio_df = ratio_df.sort_values(by='ratio', ascending=False)
ratio_df = ratio_df.dropna()

matplotlib.rcParams['pdf.fonttype'] = 42
fig, ax = plt.subplots(figsize=(3, 4))
sns.barplot(data=ratio_df, x='compartment', y='ratio', color='grey', ax=ax)
sns.despine()
ax.set_ylim(-2, 2)
plt.xticks(rotation=45, ha='right')
plt.tight_layout()
plt.savefig(os.path.join(SUPPLEMENTARY_FIG_DIR, 'supp_figure_13g.pdf'))
plt.show()

## feature correlation plots
rnaseq_feature = "tgfb_avg_exp"
mibi_feature = "HLADR+__B__cancer_border"

corr, pval = sp.spearmanr(rnaseq_dat_filtered[rnaseq_feature], mibi_dat_filtered[mibi_feature], nan_policy='omit')

plt.figure(figsize=(5,5))
sns.regplot(x=rnaseq_dat_filtered[rnaseq_feature], y=mibi_dat_filtered[mibi_feature])
# Add correlation annotation
plt.text(0.05, 0.95, f'Spearman r = {corr:.3f}',
         transform=plt.gca().transAxes,  # Position relative to plot
         verticalalignment='top',
         bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

plt.xlabel('RNA-seq feature: '+rnaseq_feature)
plt.ylabel('MIBI feature: '+mibi_feature)
plt.tight_layout()
os.makedirs(os.path.join(SUPPLEMENTARY_FIG_DIR, 'ED_figure_6'), exist_ok=True)
plt.savefig(os.path.join(SUPPLEMENTARY_FIG_DIR, 'ED_figure_6', 'corr_'+rnaseq_feature+'_'+mibi_feature+'.pdf'))
plt.show()
