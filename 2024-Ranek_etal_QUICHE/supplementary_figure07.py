import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import os
from alpineer import io_utils
import supplementary_plot_helpers

BASE_DIR = "/Volumes/Shared/Noah Greenwald/TNBC_Cohorts/SPAIN"
SUPPLEMENTARY_FIG_DIR = os.path.join("publications", "supplementary_figures", "supplementary_figure07")
# Generate overlays of entire panel across representative images
panel_validation_viz_dir = os.path.join(SUPPLEMENTARY_FIG_DIR, "supp_figure5_single_plex_overlay")
if not os.path.exists(panel_validation_viz_dir):
    os.makedirs(panel_validation_viz_dir)

exclude_chans = ["Au", "CD11c_nuc_exclude", "CK17_smoothed", "ECAD_smoothed", "FOXP3_nuc_include","LAG3",
                 "Noodle", "chan_39", "chan_45", "chan_48", "chan_115", "chan_141", "Calprotectin_old"]

# plot two control FOVs
controls_dir = os.path.join(BASE_DIR, "image_data/controls")
test_controls_fov = io_utils.list_folders(controls_dir)[0]
controls_channels = sorted(io_utils.remove_file_extensions(io_utils.list_files(os.path.join(controls_dir, test_controls_fov), substrs=".tiff")))
for ec in exclude_chans:
    if ec in controls_channels:
        controls_channels.remove(ec)

controls_fovs = ["TMA36_LN_top", "TMA40_Spleen_bottom"]
for cf in controls_fovs:
    supplementary_plot_helpers.validate_panel(controls_dir, cf, panel_validation_viz_dir, channels=controls_channels, num_rows=3)

# plot three example FOVs
samples_dir = os.path.join(BASE_DIR, "image_data", "samples")
test_samples_fov = io_utils.list_folders(samples_dir)[0]
for ec in exclude_chans:
    if ec in controls_channels:
        controls_channels.remove(ec)
sample_fovs = ["TMA32_R4C1", "TMA33_R6C4", "TMA32_R6C2"]
for sf in sample_fovs:
    supplementary_plot_helpers.validate_panel(
        samples_dir, sf, panel_validation_viz_dir,
        channels=controls_channels, num_rows=3)