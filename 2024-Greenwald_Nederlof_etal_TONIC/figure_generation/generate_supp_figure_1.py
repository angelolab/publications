import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

import os


from alpineer import io_utils
import python_files.supplementary_plot_helpers as supplementary_plot_helpers


BASE_DIR = "/Volumes/Shared/Noah Greenwald/TONIC_Cohort/"
SUPPLEMENTARY_FIG_DIR = os.path.join(BASE_DIR, "supplementary_figs")


# Generate overlays of entire panel across representative images
panel_validation_viz_dir = os.path.join(SUPPLEMENTARY_FIG_DIR, "supp_figure_1_overlays_text")
if not os.path.exists(panel_validation_viz_dir):
    os.makedirs(panel_validation_viz_dir)
exclude_chans = ["Au", "CD11c_nuc_exclude", "CK17_smoothed", "ECAD_smoothed", "FOXP3_nuc_include",
                 "LAG3", "Noodle", "chan_39", "chan_45", "chan_48", "chan_115", "chan_141"]

# plot two control FOVs
controls_dir = os.path.join(BASE_DIR, "image_data/controls")
test_controls_fov = io_utils.list_folders(controls_dir)[0]
controls_channels = sorted(io_utils.remove_file_extensions(
    io_utils.list_files(os.path.join(controls_dir, test_controls_fov), substrs=".tiff")
))
for ec in exclude_chans:
    if ec in controls_channels:
        controls_channels.remove(ec)
controls_fovs = ["TONIC_TMA6_ln_top", "TONIC_TMA14_NKI_Spleen1"]
for cf in controls_fovs:
    supplementary_plot_helpers.validate_panel(
        controls_dir, cf, panel_validation_viz_dir,
        channels=controls_channels, num_rows=3
    )

# plot two sample FOVs
samples_dir = os.path.join(BASE_DIR, "image_data/samples")
test_samples_fov = io_utils.list_folders(samples_dir)[0]
samples_channels = sorted(io_utils.remove_file_extensions(
    io_utils.list_files(os.path.join(samples_dir, test_samples_fov), substrs=".tiff")
))
for ec in exclude_chans:
    if ec in samples_channels:
        samples_channels.remove(ec)
sample_fovs = ["TONIC_TMA5_R1C2", "TONIC_TMA16_R1C3"]
for sf in sample_fovs:
    supplementary_plot_helpers.validate_panel(
        samples_dir, sf, panel_validation_viz_dir,
        channels=samples_channels, num_rows=3
    )
