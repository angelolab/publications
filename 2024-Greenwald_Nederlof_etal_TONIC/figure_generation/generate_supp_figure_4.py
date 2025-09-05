import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

import matplotlib.pyplot as plt
import seaborn as sns

from pathlib import Path
import numpy as np
import os
import pandas as pd

import python_files.supplementary_plot_helpers as supplementary_plot_helpers

BASE_DIR = "/Volumes/Shared/Noah Greenwald/TONIC_Cohort/"
raw_dir = "/Volumes/Shared/Noah Greenwald/TONIC_Acquisition/"
seg_dir = os.path.join(BASE_DIR, "segmentation_data")
image_dir = os.path.join(BASE_DIR, "image_data/samples")
SUPPLEMENTARY_FIG_DIR = os.path.join(BASE_DIR, "supplementary_figs")

# Segmentation channels and overlays
save_dir = Path(SUPPLEMENTARY_FIG_DIR) / "supp_figure_4_tiles"
save_dir.mkdir(exist_ok=True, parents=True)

membrane_channels = ["CD14", "CD38", "CD45", "ECAD", "CK17"]
overlay_channels = ["membrane_channel", "nuclear_channel"]

fovs_mem_markers = [
    "TONIC_TMA23_R10C2",
]
for fov in fovs_mem_markers:
    p = supplementary_plot_helpers.MembraneMarkersSegmentationPlot(
        fov=fov,
        image_data=image_dir,
        segmentation_dir=seg_dir,
        membrane_channels=membrane_channels,
        overlay_channels=overlay_channels,
        q=(0.05, 0.95),
        clip=False,
        figsize=(8,4),
        layout="constrained",
        image_type="pdf"
    )
    p.make_plot(save_dir=save_dir)

fovs_seg = [
    "TONIC_TMA8_R1C2",
    "TONIC_TMA9_R4C4",
    "TONIC_TMA12_R7C6",
    "TONIC_TMA21_R2C1",
    "TONIC_TMA24_R2C6",
]

# generate overlay with selected FOVs with segmentation mask and channel overlay
for fov in fovs_seg:
    p = supplementary_plot_helpers.SegmentationOverlayPlot(
        fov=fov,
        segmentation_dir=seg_dir,
        overlay_channels=overlay_channels,
        q=(0.05, 0.95),
        figsize=(8, 4),
        clip=False,
        layout="constrained",
        image_type="pdf",
    )
    p.make_plot(save_dir = save_dir)

# fov cell counts
cell_table = pd.read_csv(os.path.join(BASE_DIR, 'analysis_files/cell_table_clusters.csv'))

cluster_counts = np.unique(cell_table.fov, return_counts=True)[1]
plt.figure(figsize=(8, 6))
g = sns.histplot(data=cluster_counts, kde=True)
sns.despine()
plt.title("Histogram of Cell Counts per Image")
plt.xlabel("Number of Cells in an Image")
plt.tight_layout()
plt.savefig(os.path.join(SUPPLEMENTARY_FIG_DIR, "supp_figure_4d.pdf"), dpi=300)
plt.close()
