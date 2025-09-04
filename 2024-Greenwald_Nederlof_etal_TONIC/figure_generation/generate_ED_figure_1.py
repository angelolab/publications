import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

import os
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from upsetplot import from_contents, plot

from pathlib import Path
import python_files.supplementary_plot_helpers as supplementary_plot_helpers
from venny4py.venny4py import venny4py


BASE_DIR = "/Volumes/Shared/Noah Greenwald/TONIC_Cohort/"
SUPPLEMENTARY_FIG_DIR = os.path.join(BASE_DIR, "supplementary_figs")
seg_dir = os.path.join(BASE_DIR, "segmentation_data/deepcell_output")
metadata_dir = os.path.join(BASE_DIR, 'intermediate_files/metadata')
sequence_dir = os.path.join(BASE_DIR, 'sequencing_data')

# HnE Core, FOV and Segmentation Overlays
hne_fovs = [
    "TONIC_TMA2_R7C4",
    "TONIC_TMA4_R12C4"]

hne_path = Path(SUPPLEMENTARY_FIG_DIR) / "supp_figure_1b_masks"

save_dir = Path(SUPPLEMENTARY_FIG_DIR) / "supp_figure_1b"
save_dir.mkdir(exist_ok=True, parents=True)
for fov in hne_fovs:
    supplementary_plot_helpers.CorePlot(
        fov=fov, hne_path=hne_path, seg_dir=seg_dir
    ).make_plot(save_dir=save_dir)

# ROIs per timepoint
metadata = pd.read_csv(os.path.join(BASE_DIR, 'analysis_files/harmonized_metadata.csv'))
metadata = metadata.loc[metadata.MIBI_data_generated, :]
metadata = metadata.loc[metadata.Timepoint.isin(['primary', 'baseline', 'pre_nivo', 'on_nivo']), :]

# FOV counts per patient and timepoint
clinical_data = pd.read_csv(os.path.join(BASE_DIR, 'intermediate_files/metadata/patient_clinical_data.csv'))
metadata = metadata[['Patient_ID', 'Timepoint', 'MIBI_data_generated', 'Tissue_ID', 'fov', 'rna_seq_sample_id']].merge(clinical_data, on='Patient_ID')
metadata = metadata[metadata.Clinical_benefit.isin(['Yes', 'No'])]

fov_counts = metadata.groupby('Tissue_ID').size().values
fov_counts = pd.DataFrame(fov_counts, columns=['FOV Count'])
sns.histplot(data=fov_counts, x='FOV Count', bins=np.array(range(1,8))-0.5, shrink=0.7)
sns.despine()
plt.title("Number of FOVs per Timepoint")
plt.xlabel("Number of FOVs")
plt.tight_layout()

plt.savefig(os.path.join(SUPPLEMENTARY_FIG_DIR, "supp_figure_1c.pdf"), dpi=300)
plt.close()

# create upset plot of timepoint sample overlap
timepoint_metadata = pd.read_csv(os.path.join(metadata_dir, 'TONIC_data_per_timepoint.csv'))
clinical_data = pd.read_csv(os.path.join(BASE_DIR, 'intermediate_files/metadata/patient_clinical_data.csv'))

timepoint_metadata = timepoint_metadata.loc[timepoint_metadata.MIBI_data_generated, :]
timepoint_metadata = timepoint_metadata[timepoint_metadata.Timepoint.isin(['primary', 'baseline', 'pre_nivo', 'on_nivo'])]
timepoint_metadata = timepoint_metadata[['Patient_ID', 'Timepoint']].merge(clinical_data, on='Patient_ID')
timepoint_metadata = timepoint_metadata[timepoint_metadata.Clinical_benefit.isin(['Yes', 'No'])]
timepoint_metadata = timepoint_metadata.drop_duplicates()

# create default upset plot
contents_dict = {'Primary': timepoint_metadata.loc[timepoint_metadata.Timepoint == 'primary', 'Patient_ID'].values,
                 'Baseline': timepoint_metadata.loc[timepoint_metadata.Timepoint == 'baseline', 'Patient_ID'].values,
                 'Pre nivo': timepoint_metadata.loc[timepoint_metadata.Timepoint == 'pre_nivo', 'Patient_ID'].values,
                 'On nivo': timepoint_metadata.loc[timepoint_metadata.Timepoint == 'on_nivo', 'Patient_ID'].values}

upset_data = from_contents(contents_dict)
plot(upset_data, sort_categories_by='-input')
plt.savefig(os.path.join(SUPPLEMENTARY_FIG_DIR, 'supp_figure_1d.pdf'), dpi=300)


# venn diagram of modalities across timepoints
timepoint_metadata = pd.read_csv(os.path.join(metadata_dir, 'TONIC_data_per_timepoint.csv'))
harmonized_metadata = pd.read_csv(os.path.join(metadata_dir, 'harmonized_metadata.csv'))

# load data
mibi_metadata = timepoint_metadata.loc[timepoint_metadata.MIBI_data_generated, :]
wes_metadata = pd.read_csv(os.path.join(sequence_dir, 'preprocessing/TONIC_WES_meta_table.tsv'), sep='\t')
rna_metadata = pd.read_csv(os.path.join(sequence_dir, 'preprocessing/TONIC_tissue_rna_id.tsv'), sep='\t')
rna_metadata = rna_metadata.merge(harmonized_metadata[['Patient_ID', 'Tissue_ID', 'Timepoint']].drop_duplicates(), on='Tissue_ID', how='left')

clinical_data = pd.read_csv(os.path.join(BASE_DIR, 'intermediate_files/metadata/patient_clinical_data.csv'))

mibi_metadata = mibi_metadata[mibi_metadata.Timepoint.isin(['primary', 'baseline', 'pre_nivo', 'on_nivo'])]
mibi_metadata = mibi_metadata[['Patient_ID', 'Timepoint']].merge(clinical_data, on='Patient_ID')
mibi_metadata = mibi_metadata[mibi_metadata.Clinical_benefit.isin(['Yes', 'No'])]
wes_metadata = wes_metadata.rename(columns={'Individual.ID': 'Patient_ID', 'timepoint': 'Timepoint'})
wes_metadata = wes_metadata[wes_metadata.Clinical_benefit.isin(['Yes', 'No'])]
rna_metadata = rna_metadata[['Patient_ID', 'Timepoint']].merge(clinical_data, on='Patient_ID').drop_duplicates()
rna_metadata = rna_metadata[rna_metadata.Clinical_benefit.isin(['Yes', 'No'])]

# separate venn diagram per timepoint
for timepoint, plot_name in zip(['baseline', 'pre_nivo', 'on_nivo'], ["e", "f", "g"]):
    mibi_ids = set(mibi_metadata.loc[mibi_metadata.Timepoint == timepoint, 'Patient_ID'].values)
    wes_ids = set(wes_metadata.loc[wes_metadata.Timepoint == timepoint, 'Patient_ID'].values)
    rna_ids = set(rna_metadata.loc[rna_metadata.Timepoint == timepoint, 'Patient_ID'].values)

    if timepoint == 'baseline':
        sets = {
            'MIBI': mibi_ids,
            'WES': wes_ids,
            'RNA': rna_ids}
        venny4py(sets=sets, colors="bgr")
    else:
        sets = {
            'MIBI': mibi_ids,
            'RNA': rna_ids}
        venny4py(sets=sets, colors="br")

    plt.savefig(os.path.join(SUPPLEMENTARY_FIG_DIR, 'supp_figure_1{}.pdf'.format(plot_name)), dpi=300,
                bbox_inches='tight')
    plt.close()
