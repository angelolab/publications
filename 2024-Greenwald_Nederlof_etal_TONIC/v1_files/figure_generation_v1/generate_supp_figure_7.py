import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

import os
import pandas as pd

import python_files.supplementary_plot_helpers as supplementary_plot_helpers

BASE_DIR = "/Volumes/Shared/Noah Greenwald/TONIC_Cohort/"
SUPPLEMENTARY_FIG_DIR = os.path.join(BASE_DIR, "supplementary_figs")
seg_dir = os.path.join(BASE_DIR, "segmentation_data/deepcell_output")
image_dir = os.path.join(BASE_DIR, "image_data/samples/")

# Feature parameter tuning
extraction_pipeline_tuning_dir = os.path.join(SUPPLEMENTARY_FIG_DIR, "supp_figure_7_robustness")
if not os.path.exists(extraction_pipeline_tuning_dir):
    os.makedirs(extraction_pipeline_tuning_dir)

# vary the features for each marker threshold
cell_table_full = pd.read_csv(
    os.path.join(BASE_DIR, "analysis_files/combined_cell_table_normalized_cell_labels_updated.csv")
)
supplementary_plot_helpers.run_functional_marker_positivity_tuning_tests(
    cell_table_full, extraction_pipeline_tuning_dir, supplementary_plot_helpers.MARKER_INFO,
    threshold_mults=[1/4, 1/2, 3/4, 7/8, 1, 8/7, 4/3, 2, 4]
)

# vary min cell param to see how many FOVs get kept or not
cluster_broad_df = pd.read_csv(os.path.join(BASE_DIR, "output_files/cluster_df_per_core.csv"))
supplementary_plot_helpers.run_min_cell_feature_gen_fovs_dropped_tests(
    cluster_broad_df, min_cell_params=[1, 3, 5, 10, 20], compartments=["all"],
    metrics=["cluster_broad_count"], save_dir=extraction_pipeline_tuning_dir
)

# vary params for cancer mask and boundary definition inclusion
cell_table_clusters = pd.read_csv(os.path.join(BASE_DIR, "analysis_files/cell_table_clusters.csv"))
supplementary_plot_helpers.run_cancer_mask_inclusion_tests(
    cell_table_clusters, channel_dir=image_dir, seg_dir=seg_dir,
    threshold_mults=[1/4, 1/2, 3/4, 7/8, 1, 8/7, 4/3, 2, 4],
    save_dir=extraction_pipeline_tuning_dir, base_sigma=10, base_channel_thresh=0.0015,
    base_min_mask_size=7000, base_max_hole_size=1000, base_border_size=50
)
