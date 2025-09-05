import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

import os
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

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

fov_counts = metadata.groupby('Tissue_ID').size().values
fov_counts = pd.DataFrame(fov_counts, columns=['FOV Count'])
sns.histplot(data=fov_counts, x='FOV Count')
sns.despine()
plt.title("Number of FOVs per Timepoint")
plt.xlabel("Number of FOVs")
plt.tight_layout()

plt.savefig(os.path.join(SUPPLEMENTARY_FIG_DIR, "supp_figure_1c.pdf"), dpi=300)
plt.close()


# venn diagram of modalities across timepoints
timepoint_metadata = pd.read_csv(os.path.join(metadata_dir, 'TONIC_data_per_timepoint.csv'))
harmonized_metadata = pd.read_csv(os.path.join(metadata_dir, 'harmonized_metadata.csv'))

# load data
mibi_metadata = timepoint_metadata.loc[timepoint_metadata.MIBI_data_generated, :]
wes_metadata = pd.read_csv(os.path.join(sequence_dir, 'preprocessing/TONIC_WES_meta_table.tsv'), sep='\t')
rna_metadata = pd.read_csv(os.path.join(sequence_dir, 'preprocessing/TONIC_tissue_rna_id.tsv'), sep='\t')
rna_metadata = rna_metadata.merge(harmonized_metadata[['Patient_ID', 'Tissue_ID', 'Timepoint']].drop_duplicates(), on='Tissue_ID', how='left')

# separate venn diagram per timepoint
for timepoint, plot_name in zip(['baseline', 'pre_nivo', 'on_nivo'], ["e", "f", "g"]):
    mibi_ids = set(mibi_metadata.loc[mibi_metadata.Timepoint == timepoint, 'Patient_ID'].values)
    wes_ids = set(wes_metadata.loc[wes_metadata.timepoint == timepoint, 'Individual.ID'].values)
    rna_ids = set(rna_metadata.loc[rna_metadata.Timepoint == timepoint, 'Patient_ID'].values)

    sets = {
        'MIBI': mibi_ids,
        'WES': wes_ids,
        'RNA': rna_ids}

    venny4py(sets=sets)
    plt.savefig(os.path.join(SUPPLEMENTARY_FIG_DIR, 'supp_figure_1{}.pdf'.format(plot_name)), dpi=300, bbox_inches='tight')
    plt.close()