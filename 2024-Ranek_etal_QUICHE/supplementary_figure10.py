import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import numpy as np
import os
import pandas as pd
import supplementary_plot_helpers

BASE_DIR = "/Volumes/Shared/Noah Greenwald/TNBC_Cohorts/SPAIN"
raw_dir = "/Volumes/Shared/Noah Greenwald/TNBC_Acquisition/"
SUPPLEMENTARY_FIG_DIR = os.path.join("publications", "supplementary_figures", "supplementary_figure10")
seg_dir = os.path.join(BASE_DIR, "segmentation", "samples")
image_dir = os.path.join(BASE_DIR, "image_data", "samples")

# Segmentation channels and overlays
save_dir = Path(SUPPLEMENTARY_FIG_DIR) / "supp_figure10_tiles"

if not os.path.exists(save_dir):
    os.makedirs(save_dir)

membrane_channels = ["CD14", "CD38", "CD45", "ECAD", "CK17"]
overlay_channels = ["membrane_channel", "nuclear_channel"]

fovs_mem_markers = ["TMA31_R6C6"]

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
        image_type="pdf")
    p.make_plot(save_dir=save_dir)

fovs_seg = ["TMA31_R6C6",
    "TMA32_R1C3",
    "TMA32_R3C9",
    "TMA32_R5C8",
    "TMA33_R1C4"]

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
cell_table = pd.read_csv(os.path.join(BASE_DIR, 'intermediate_files', 'post_processing', 'cell_table_clusters.csv'))

cluster_counts = np.unique(cell_table.fov, return_counts=True)[1]
plt.figure(figsize=(8, 6))
sns.set_style('ticks')
g = sns.histplot(data=cluster_counts, kde=True)
g.tick_params(labelsize=16)
sns.despine()
plt.xlabel("Cells Per Image", fontsize = 16)
plt.ylabel("Count", fontsize = 16)
plt.tight_layout()
plt.savefig(os.path.join(SUPPLEMENTARY_FIG_DIR, "supp_figure_10c.pdf"), dpi=300)
plt.close()