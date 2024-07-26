import os
import re

import gseapy as gp
from gseapy import Msigdb


import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import numpy as np


base_dir = '/Volumes/Shared/Noah Greenwald/TONIC_Cohort'
sequence_dir = os.path.join(base_dir, 'sequencing_data')

harmonized_metadata = pd.read_csv(os.path.join(base_dir, 'analysis_files/all_fovs/harmonized_metadata.csv'))

# anonymize patient IDs
patient_ID_mapping = pd.read_csv(os.path.join(base_dir, 'final_patient_mapping.csv'))
harmonized_metadata = pd.merge(harmonized_metadata, patient_ID_mapping, on='Patient_ID', how='left')
harmonized_metadata = harmonized_metadata.drop(columns=['Patient_ID'])
harmonized_metadata = harmonized_metadata.rename(columns={'new_Patient_ID': 'Patient_ID'})

#
# calculate gene set enrichment scores
#

# Load gene expression data
gene_features = pd.read_csv(os.path.join(sequence_dir, 'preprocessing/TONIC_gene_expression_TPM_table.tsv'), sep='\t', index_col=0)

# Calculate msigdb genesets
msig = Msigdb()
categories_keep = ["h.all", "c2.all", "c5.go"]
all_gene_sets = {}
for one_category in categories_keep:
    one_set = msig.get_gmt(category=one_category, dbver="2023.2.Hs")
    all_gene_sets.update(one_set)

# Keep specific gene sets
gene_sets_keep = pd.read_csv(os.path.join(sequence_dir, 'preprocessing/msigdb_genesets.csv'))
gene_sets_names_keep = gene_sets_keep['name'].values
gene_sets = {x: all_gene_sets[x] for x in gene_sets_names_keep}

# Do ssgsea
ss = gp.ssgsea(data = gene_features,
               gene_sets = gene_sets,
               sample_norm_method = 'rank')

ss.res2d.to_csv(os.path.join(sequence_dir, 'preprocessing/msigdb_ssgsea_es.csv'), index=False)

# calculate other gene sets

# Load data
metadata_sub = harmonized_metadata[["rna_seq_sample_id","Timepoint"]]
metadata_sub = metadata_sub.dropna()
metadata_sub = metadata_sub.drop_duplicates()

all_genes = gene_features.index.to_list()

# Initialize data table to store all the gene signature scores
all_gene_scores_dt = pd.DataFrame({'rna_seq_sample_id': metadata_sub.rna_seq_sample_id.values})

# TGFb signature from Mariathasan et al.
# TGFb in fibroblasts
geneset_name = "pan_f_tbrs"
# Direct from paper
gene_set = ["ACTA2", "ACTG2", "ADAM12", "ADAM19", "CNN1", "COL4A1", "CTGF", "CTPS1", "FAM101B", "FSTL3", "HSPB1", "IGFBP3", "PXDC1", "SEMA7A", "SH3PXD2A", "TAGLN", "TGFBI", "TNS1", "TPM1"]
# Change a few gene names to match data, CTGF -> CCN2, FAM101B -> RFLNB
gene_set = ["ACTA2", "ACTG2", "ADAM12", "ADAM19", "CNN1", "COL4A1", "CCN2", "CTPS1", "RFLNB", "FSTL3", "HSPB1", "IGFBP3", "PXDC1", "SEMA7A", "SH3PXD2A", "TAGLN", "TGFBI", "TNS1", "TPM1"]
# Check to make sure they are all included in feature table
set(gene_set) <= set(all_genes)


# First PC as gene score (as described in paper)
# Subset data for gene set
gene_features_sub = gene_features.loc[gene_set]
# Take transpose
gene_features_sub = gene_features_sub.transpose()

# PCA
pca_in = StandardScaler().fit_transform(gene_features_sub)
pca = PCA(n_components=2)
pca_out = pca.fit_transform(pca_in)

# Save as dataframe
out_dat = pd.DataFrame({"rna_seq_sample_id": gene_features_sub.index.values,
                        "pc1": pca_out[:, 0]})

# Save to dataframe to store all outputs
all_gene_scores_dt = pd.merge(all_gene_scores_dt, out_dat, on='rna_seq_sample_id')
all_gene_scores_dt = all_gene_scores_dt.rename(columns={'pc1': 'tgfb_pc1'})

# Merge with metadata
out_dat = pd.merge(out_dat, metadata_sub, on='rna_seq_sample_id')

# Get average of genes
gene_features_sub = gene_features.loc[gene_set]
means = gene_features_sub.mean()

# Save as dataframe
means_df = pd.DataFrame({"rna_seq_sample_id":gene_features_sub.columns.values,
                         "avg_exp":means.values})

# Save to dataframe to store all outputs
all_gene_scores_dt = pd.merge(all_gene_scores_dt, means_df, on='rna_seq_sample_id')
all_gene_scores_dt = all_gene_scores_dt.rename(columns={'avg_exp':'tgfb_avg_exp'})

# Merge with metadata
out_dat = pd.merge(out_dat, means_df, on='rna_seq_sample_id')

# cDC1/cDC2 signatures
cDC1_genes = ["IRF8","BATF3","THBD"]
cDC2_genes = ["ZEB2","IRF4","KLF4","CD1C","ITGAX","ITGAM"]

# Get average of genes
#gene_set = cDC1_genes
gene_set = cDC2_genes
gene_features_sub = gene_features.loc[gene_set]
means = gene_features_sub.mean()

# Save as dataframe
means_df = pd.DataFrame({"rna_seq_sample_id":gene_features_sub.columns.values,
                         "avg_exp":means.values})

# Save to dataframe to store all outputs
all_gene_scores_dt = pd.merge(all_gene_scores_dt, means_df, on='rna_seq_sample_id')
#all_gene_scores_dt = all_gene_scores_dt.rename(columns={'avg_exp':'cDC1_avg_exp'})
all_gene_scores_dt = all_gene_scores_dt.rename(columns={'avg_exp':'cDC2_avg_exp'})

# Trm signature
# from Savas et al, https://doi.org/10.1038/s41591-018-0078-7
trm_sig = pd.read_csv(os.path.join(sequence_dir, 'preprocessing/trm_sig_savas.csv'))
gene_set = trm_sig['gene_names'].values

# Check to make sure they are all included in feature table
set(gene_set) <= set(all_genes)
# Only keep those that exist in our dataset
gene_set = [x for x in gene_set if x in all_genes]
set(gene_set) <= set(all_genes)

# Get average of genes
gene_features_sub = gene_features.loc[gene_set]
means = gene_features_sub.mean()

# Save as dataframe
means_df = pd.DataFrame({"rna_seq_sample_id":gene_features_sub.columns.values,
                         "avg_exp":means.values})

# Save to dataframe to store all outputs
all_gene_scores_dt = pd.merge(all_gene_scores_dt, means_df, on=['rna_seq_sample_id'])
all_gene_scores_dt = all_gene_scores_dt.rename(columns={'avg_exp':'Trm_avg_exp'})

# Merge with metadata
out_dat = pd.merge(metadata_sub, means_df, on='rna_seq_sample_id')

# from Tietscher et al, https://doi.org/10.1038/s41467-022-35238-w
t_cell_attractive = ["CCL21","CCL17","CCL2","CXCL9","CXCL10","CXCL11","CXCL12","CCL3","CCL4","CCL5","CXCL16"]
t_cell_suppressive = ["IDO1","CD274","PDCD1LG2","TNFSF10","HLA-E","HLA-G","VTCN1","IL10","TGFB1","TGFB2","PTGS2","PTGES","LGALS9","CCL22","CD80","CD86"]

# Get average of genes
#gene_set = t_cell_attractive
gene_set = t_cell_suppressive
gene_features_sub = gene_features.loc[gene_set]
means = gene_features_sub.mean()

# Save as dataframe and merge with metadata
means_df = pd.DataFrame({"rna_seq_sample_id":gene_features_sub.columns.values,
                         "avg_exp":means.values})

# Save to dataframe to store all outputs
all_gene_scores_dt = pd.merge(all_gene_scores_dt, means_df, on=['rna_seq_sample_id'])
#all_gene_scores_dt = all_gene_scores_dt.rename(columns={'avg_exp':'tcell_attractive_avg_exp'})
all_gene_scores_dt = all_gene_scores_dt.rename(columns={'avg_exp':'tcell_suppressive_avg_exp'})

## T cell exhaustion
# from Doering et al, https://doi.org/10.1016/j.immuni.2012.08.021
gene_set = ["RTP4", "FOXP1", "IKZF2", "ZEB2", "LASS6", "TOX", "EOMES"]
set(gene_set) <= set(all_genes)
# change one gene name
gene_set = ["RTP4", "FOXP1", "IKZF2", "ZEB2", "CERS6", "TOX", "EOMES"]
set(gene_set) <= set(all_genes)

# Get average of genes
gene_features_sub = gene_features.loc[gene_set]
means = gene_features_sub.mean()

# Save as dataframe
means_df = pd.DataFrame({"rna_seq_sample_id":gene_features_sub.columns.values,
                         "avg_exp":means.values})
# Save to dataframe to store all outputs
all_gene_scores_dt = pd.merge(all_gene_scores_dt, means_df, on=['rna_seq_sample_id'])
all_gene_scores_dt = all_gene_scores_dt.rename(columns={'avg_exp':'tcell_exhaustion_avg_exp'})

# Merge with metadata
out_dat = pd.merge(metadata_sub, means_df, on='rna_seq_sample_id')
all_gene_scores_dt.to_csv(os.path.join(sequence_dir, 'preprocessing/misc_gene_set_scores.csv'), index=False)


#
# clean up genomics features
#

# load data
genomics_feature_df = pd.read_csv(os.path.join(sequence_dir, 'preprocessing/TONIC_WES_meta_table.tsv'), sep='\t')
genomics_feature_df = genomics_feature_df.rename(columns={'Individual.ID': 'Patient_ID', 'timepoint': 'Timepoint'})

# anonymize patient IDs
patient_ID_mapping = pd.read_csv(os.path.join(base_dir, 'final_patient_mapping.csv'))
genomics_feature_df = pd.merge(genomics_feature_df, patient_ID_mapping, on='Patient_ID', how='left')
genomics_feature_df = genomics_feature_df.drop(columns=['Patient_ID'])
genomics_feature_df = genomics_feature_df.rename(columns={'new_Patient_ID': 'Patient_ID'})

# add in new IC10 predictions
new_ic10 = pd.read_csv(os.path.join(sequence_dir, 'preprocessing/tonic_predictions.txt'), sep='\t')
new_ic10['Experiment.System.ID'] = new_ic10['Sample'].apply(lambda x: x.split('_vs')[0])
new_ic10 = new_ic10.rename(columns={'voting': 'Eniclust_2'})
genomics_feature_df = pd.merge(genomics_feature_df, new_ic10[['Experiment.System.ID', 'Eniclust_2']], on='Experiment.System.ID', how='left')

# drop columns which aren't needed for analysis
drop_cols = ['Localization', 'induction_arm', 'Experiment.System.ID', 'B2', 'Clinical_benefit']
genomics_feature_df = genomics_feature_df.drop(columns=drop_cols)

# check which columns have non-numeric values
non_numeric_cols = []
for col in genomics_feature_df.columns:
    if genomics_feature_df[col].dtype == 'object':
        non_numeric_cols.append(col)

# recode mutation features to yes/no
binarize_muts = [col for col in genomics_feature_df.columns if '_SNV' in col]
for mut in binarize_muts:
    genomics_feature_df[mut] = (~genomics_feature_df[mut].isna()).astype(int)

# change binary features to 0/1
genomics_feature_df['tail'] = (genomics_feature_df['tail'] == True).astype(int)
genomics_feature_df['wgd'] = (genomics_feature_df['wgd'] == True).astype(int)
genomics_feature_df['hla_loh'] = (genomics_feature_df['hla_loh'] == 'yes').astype(int)

# recode copy number features
amp_cols = [col for col in genomics_feature_df.columns if '_AMP' in col]
del_cols = [col for col in genomics_feature_df.columns if '_DEL' in col]
loh_cols = [col for col in genomics_feature_df.columns if '_LOH' in col]
cn_cols = amp_cols + del_cols + loh_cols

for col in cn_cols:
    genomics_feature_df[col] = (genomics_feature_df[col].values == 'yes').astype(int)

# simplify categorical features
np.unique(genomics_feature_df['Eniclust_2'].values, return_counts=True)
keep_dict = {'IC10_rna': [4, 9, 10], 'IC10_rna_cna': [4., 9., 10.], 'PAM50': ['Basal', 'Her2'],
             'IntClust': ['ic10', 'ic4', 'ic9'], 'A1': ['A*01:01', 'A*02:01', 'A*03:01'],
             'A2': ['A*02:01', 'A*03:01', 'A*24:02'], 'B1': ['B*07:02', 'B*08:01', 'B*40:01'],
             'C1': ['C*03:04', 'C*04:01', 'C*07:01', 'C*07:02'], 'C2': ['C*07:01', 'C*07:02'],
             'Eniclust_2': ['ic4', 'ic9', 'ic10']}

for col in keep_dict.keys():
    genomics_feature_df[col] = genomics_feature_df[col].apply(lambda x: x if x in keep_dict[col] else 'Other')
    genomics_feature_df[col] = genomics_feature_df[col].astype('category')

# add computed features
genomics_feature_df['subclonal_clonal_ratio'] = genomics_feature_df['SUBCLONAL'] / genomics_feature_df['CLONAL']

# change categorical features to dummy variables
genomics_feature_df = pd.get_dummies(genomics_feature_df, columns=keep_dict.keys(), drop_first=True)

# save processed genomics features
genomics_feature_df.to_csv(os.path.join(sequence_dir, 'preprocessing/TONIC_WES_meta_table_processed.csv'), index=False)


#
# process RNA-seq data
#

RNA_feature_df = pd.read_csv(os.path.join(sequence_dir, 'preprocessing/TONIC_immune_sig_score_table.tsv'), sep='\t')
RNA_feature_df_genes = pd.read_csv(os.path.join(sequence_dir, 'preprocessing/TONIC_immune_sig_gene_TPM_table.tsv'), sep='\t')
RNA_feature_df = pd.merge(RNA_feature_df, RNA_feature_df_genes, on='sample_identifier', how='left')
RNA_feature_df = RNA_feature_df.rename(columns={'sample_identifier': 'rna_seq_sample_id'})
RNA_feature_df = RNA_feature_df.merge(harmonized_metadata[['rna_seq_sample_id', 'Timepoint', 'Patient_ID', 'Tissue_ID']].drop_duplicates(),
                                    on='rna_seq_sample_id', how='left')
RNA_feature_df = RNA_feature_df.drop(columns=['rna_seq_sample_id'])
RNA_feature_df['TME_subtype_immune'] = RNA_feature_df['TME_subtype'].apply(lambda x: 0 if x == 'F' or x == 'D' else 1)
RNA_feature_df['TME_subtype_fibrotic'] = RNA_feature_df['TME_subtype'].apply(lambda x: 1 if x == 'F' or x == 'IE/F' else 0)
RNA_feature_df = RNA_feature_df.drop(columns=['TME_subtype'])

RNA_feature_df.to_csv(os.path.join(sequence_dir, 'preprocessing/TONIC_immune_sig_score_and_genes_processed.csv'), index=False)

# proccess other gene scores files
msigdb_scores = pd.read_csv(os.path.join(sequence_dir, 'preprocessing/msigdb_ssgsea_es.csv'))
msigdb_scores = msigdb_scores.rename(columns={'Name': 'rna_seq_sample_id', 'Term': 'feature_name', 'NES': 'feature_value'})
msigdb_scores = msigdb_scores.merge(harmonized_metadata[['rna_seq_sample_id', 'Timepoint', 'Patient_ID', 'Tissue_ID']].drop_duplicates(),
                                    on='rna_seq_sample_id', how='left')
msigdb_scores = msigdb_scores.drop(columns=['ES', 'rna_seq_sample_id'])

other_scores = pd.read_csv(os.path.join(sequence_dir, 'preprocessing/misc_gene_set_scores.csv'))
other_scores_long = pd.melt(other_scores, id_vars=['rna_seq_sample_id'], var_name='feature_name', value_name='feature_value')
other_scores_long = other_scores_long.merge(harmonized_metadata[['rna_seq_sample_id', 'Timepoint', 'Patient_ID', 'Tissue_ID']].drop_duplicates(),
                                    on='rna_seq_sample_id', how='left')
other_scores_long = other_scores_long.drop(columns=['rna_seq_sample_id'])

# add in RS and cgas scores
rs_cgas = pd.read_csv(os.path.join(sequence_dir, 'preprocessing/tonic_RS_cGAS_STING_signatures.txt'), sep='\t')
rs_cgas = rs_cgas.rename(columns={'Sample': 'Experiment.System.ID', 'RS': 'replication_stress', 'STING': 'cgas_sting'})
rna_ids = pd.read_csv(os.path.join(sequence_dir, 'preprocessing/TONIC_tissue_rna_id.tsv'), sep='\t')
rna_ids = rna_ids.rename(columns={'system_id2': 'Experiment.System.ID'})
rs_cgas = rs_cgas.merge(rna_ids, on='Experiment.System.ID', how='left')
rs_cgas = rs_cgas.drop(columns=['Experiment.System.ID', 'sample_identifier'])
rs_cgas_long = pd.melt(rs_cgas, id_vars=['Tissue_ID'], var_name='feature_name', value_name='feature_value')

rs_cgas_long = rs_cgas_long.merge(harmonized_metadata[['Timepoint', 'Patient_ID', 'Tissue_ID']].drop_duplicates(),
                                    on='Tissue_ID', how='left')

# combine together
combined_rna_scores = pd.concat([msigdb_scores, other_scores_long, rs_cgas_long])
combined_rna_scores.to_csv(os.path.join(sequence_dir, 'preprocessing/msigdb_misc_processed.csv'), index=False)

#
# Aggregate sequencing data and label with appropriate metaedata
#

# transform to long format
genomics_feature_df = pd.read_csv(os.path.join(sequence_dir, 'preprocessing/TONIC_WES_meta_table_processed.csv'))

genomics_feature_df = pd.merge(genomics_feature_df, harmonized_metadata[['Patient_ID', 'Timepoint', 'Tissue_ID']].drop_duplicates(),
                               on=['Patient_ID', 'Timepoint'], how='left')

genomics_feature_df_long = pd.melt(genomics_feature_df, id_vars=['Patient_ID', 'Timepoint', 'Tissue_ID'],
                                   var_name='feature_name', value_name='feature_value')

# label summary metrics
summary_metrics = ['purity', 'ploidy', 'fga', 'wgd', 'frac_loh']
genomics_feature_df_long.loc[genomics_feature_df_long.feature_name.isin(summary_metrics), 'data_type'] = 'DNA'
genomics_feature_df_long.loc[genomics_feature_df_long.feature_name.isin(summary_metrics), 'feature_type'] = 'cn_summary'

# label gene cns
gene_cns = [col for col in genomics_feature_df.columns if col.endswith('_cn')]
genomics_feature_df_long.loc[genomics_feature_df_long.feature_name.isin(gene_cns), 'data_type'] = 'DNA'
genomics_feature_df_long.loc[genomics_feature_df_long.feature_name.isin(gene_cns), 'feature_type'] = 'gene_cn'

# label gene rna
gene_rna = [col for col in genomics_feature_df.columns if col.endswith('_rna')]
genomics_feature_df_long.loc[genomics_feature_df_long.feature_name.isin(gene_rna), 'data_type'] = 'RNA'
genomics_feature_df_long.loc[genomics_feature_df_long.feature_name.isin(gene_rna), 'feature_type'] = 'gene_rna'

# label mutation summary
mut_summary = ['snv_count', 'missense_count', 'synonymous_count', 'frameshift_count']
genomics_feature_df_long.loc[genomics_feature_df_long.feature_name.isin(mut_summary), 'data_type'] = 'DNA'
genomics_feature_df_long.loc[genomics_feature_df_long.feature_name.isin(mut_summary), 'feature_type'] = 'mut_summary'

# label mutations
mutations = [col for col in genomics_feature_df.columns if col.endswith('_SNV')]
genomics_feature_df_long.loc[genomics_feature_df_long.feature_name.isin(mutations), 'data_type'] = 'DNA'
genomics_feature_df_long.loc[genomics_feature_df_long.feature_name.isin(mutations), 'feature_type'] = 'mutation'

# label clonality
clones = ['CLONAL', 'SUBCLONAL', 'INDETERMINATE', 'clones', 'tail', 'subclonal_clonal_ratio']
genomics_feature_df_long.loc[genomics_feature_df_long.feature_name.isin(clones), 'data_type'] = 'DNA'
genomics_feature_df_long.loc[genomics_feature_df_long.feature_name.isin(clones), 'feature_type'] = 'clonality'

# label antigen information
antigens = ['antigens', 'antigens_clonal', 'antigens_subclonal', 'antigens_indeterminate',
             'antigens_variants', 'HBMR_norm', 'HBMR_norm_p', 'HBMR_norm_cilo', 'HBMR_norm_ciup',
             'hla_tcn', 'hla_lcn', 'hla_loh']
genomics_feature_df_long.loc[genomics_feature_df_long.feature_name.isin(antigens), 'data_type'] = 'DNA'
genomics_feature_df_long.loc[genomics_feature_df_long.feature_name.isin(antigens), 'feature_type'] = 'antigen'

# label hlas
hlas = [col for col in genomics_feature_df.columns if re.match('[ABC][12]_', col) is not None]
genomics_feature_df_long.loc[genomics_feature_df_long.feature_name.isin(hlas), 'data_type'] = 'DNA'
genomics_feature_df_long.loc[genomics_feature_df_long.feature_name.isin(hlas), 'feature_type'] = 'hla'

# label mutation signatures
mut_signatures = [col for col in genomics_feature_df.columns if col.startswith('Signature')]
genomics_feature_df_long.loc[genomics_feature_df_long.feature_name.isin(mut_signatures), 'data_type'] = 'DNA'
genomics_feature_df_long.loc[genomics_feature_df_long.feature_name.isin(mut_signatures), 'feature_type'] = 'mut_signature'

# label chromosome changes
chrom_cn = [col for col in genomics_feature_df.columns if col.endswith('AMP') or col.endswith('DEL') or col.endswith('LOH')]
genomics_feature_df_long.loc[genomics_feature_df_long.feature_name.isin(chrom_cn), 'data_type'] = 'DNA'
genomics_feature_df_long.loc[genomics_feature_df_long.feature_name.isin(chrom_cn), 'feature_type'] = 'chrom_cn'

# label IC subtypes
ic_subtypes = [col for col in genomics_feature_df.columns if col.startswith('IC') or col.startswith('IntClust') or col.startswith('Eniclust')]
genomics_feature_df_long.loc[genomics_feature_df_long.feature_name.isin(ic_subtypes), 'data_type'] = 'RNA'
genomics_feature_df_long.loc[genomics_feature_df_long.feature_name.isin(ic_subtypes), 'feature_type'] = 'ic_subtype'

# label PAM50 subtypes
pam = ['PAM50_Her2', 'PAM50_Other']
genomics_feature_df_long.loc[genomics_feature_df_long.feature_name.isin(pam), 'data_type'] = 'RNA'
genomics_feature_df_long.loc[genomics_feature_df_long.feature_name.isin(pam), 'feature_type'] = 'pam_subtype'


# transform RNA features to long format
RNA_feature_df = pd.read_csv(os.path.join(sequence_dir, 'preprocessing/TONIC_immune_sig_score_and_genes_processed.csv'))
RNA_feature_df_long = pd.melt(RNA_feature_df, id_vars=['Patient_ID', 'Timepoint', 'Tissue_ID'],
                                   var_name='feature_name', value_name='feature_value')

# label cell type signatures
cell_signatures = ['NK_cells', 'T_cells', 'B_cells', 'Treg', 'Neutrophil_signature', 'MDSC', 'Macrophages',
                   'CAF', 'Matrix', 'Endothelium', 'TME_subtype_immune', 'TME_subtype_fibrotic']
# recode fibrotic to immune yes/no, fibrotic yes/no
RNA_feature_df_long.loc[RNA_feature_df_long.feature_name.isin(cell_signatures), 'data_type'] = 'RNA'
RNA_feature_df_long.loc[RNA_feature_df_long.feature_name.isin(cell_signatures), 'feature_type'] = 'cell_signature'

# label functional signatures
functional_signatures = ['T_cell_traffic', 'MHCI', 'MHCII', 'Coactivation_molecules', 'Effector_cells',
                         'T_cell_traffic', 'M1_signatures', 'Th1_signature', 'Antitumor_cytokines',
                         'Checkpoint_inhibition', 'T_reg_traffic', 'Granulocyte_traffic', 'MDSC_traffic',
                         'Macrophage_DC_traffic', 'Th2_signature', 'Protumor_cytokines', 'Matrix_remodeling',
                         'Angiogenesis', 'Proliferation_rate', 'EMT_signature', 'CAscore', 'GEP_mean']
RNA_feature_df_long.loc[RNA_feature_df_long.feature_name.isin(functional_signatures), 'data_type'] = 'RNA'
RNA_feature_df_long.loc[RNA_feature_df_long.feature_name.isin(functional_signatures), 'feature_type'] = 'functional_signature'

# label individual genes
unlabeled = RNA_feature_df_long.loc[RNA_feature_df_long.data_type.isna(), 'feature_name'].unique()
RNA_feature_df_long.loc[RNA_feature_df_long.feature_name.isin(unlabeled), 'data_type'] = 'RNA'
RNA_feature_df_long.loc[RNA_feature_df_long.feature_name.isin(unlabeled), 'feature_type'] = 'gene_rna'
RNA_feature_df_long.loc[RNA_feature_df_long.feature_type == 'gene_rna', 'feature_name'] = RNA_feature_df_long.loc[RNA_feature_df_long.feature_type == 'gene_rna', 'feature_name'].apply(lambda x: x + '_rna')

# find RNAs already present in genomics features
RNA_names = genomics_feature_df_long.loc[genomics_feature_df_long.feature_type == 'gene_rna', 'feature_name'].unique()
RNA_feature_df_long = RNA_feature_df_long.loc[~RNA_feature_df_long.feature_name.isin(RNA_names), :]

# read in other RNA features
other_rna_features = pd.read_csv(os.path.join(sequence_dir, 'preprocessing/msigdb_misc_processed.csv'))
other_rna_features['data_type'] = 'RNA'
other_rna_features['feature_type'] = 'functional_signature_2'

# combine together
genomics_feature_df_long_combined = pd.concat([genomics_feature_df_long, RNA_feature_df_long, other_rna_features])

# normalize

# compute z-scores for each feature
genomics_feature_df_long_combined = genomics_feature_df_long_combined.rename(columns={'feature_value': 'raw_value'})
genomics_feature_df_wide = genomics_feature_df_long_combined.pivot(index=['Tissue_ID'], columns='feature_name', values='raw_value')

# remove columns that sum to 0
genomics_feature_df_wide = genomics_feature_df_wide.loc[:, genomics_feature_df_wide.sum() != 0]

zscore_df = (genomics_feature_df_wide - genomics_feature_df_wide.mean()) / genomics_feature_df_wide.std()

# add z-scores to original df
zscore_df = zscore_df.reset_index()
zscore_df_long = pd.melt(zscore_df, id_vars='Tissue_ID', var_name='feature_name', value_name='normalized_value')
genomics_feature_df_long_combined = pd.merge(genomics_feature_df_long_combined, zscore_df_long, on=['Tissue_ID', 'feature_name'], how='left')


genomics_feature_df_long_combined.to_csv(os.path.join(sequence_dir, 'processed_genomics_features.csv'), index=False)
