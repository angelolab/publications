import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import os
import supplementary_plot_helpers

BASE_DIR = "/Volumes/Shared/Noah Greenwald/TNBC_Cohorts/SPAIN"
raw_dir = "/Volumes/Shared/Noah Greenwald/TNBC_Acquisition/SPAIN"
SUPPLEMENTARY_FIG_DIR = os.path.join("publications", "supplementary_figures", "supplementary_figure08")

# show a run with images pre- and post-Rosetta
rosetta_tiling = os.path.join(SUPPLEMENTARY_FIG_DIR, "supp_figure_08_tiles")
if not os.path.exists(rosetta_tiling):
    os.makedirs(rosetta_tiling)

run_name = "2023-04-13_SPAIN_TMA32_run1"
pre_rosetta_dir = os.path.join(raw_dir, "extracted")
post_rosetta_dir = os.path.join(raw_dir, "rosetta")

# NOTE: images not scaled up programmatically, this happens manually in Photoshop
supplementary_plot_helpers.stitch_before_after_rosetta(
    pre_rosetta_dir, post_rosetta_dir, rosetta_tiling, run_name,
    [22], "CD20", pre_rosetta_subdir="", post_rosetta_subdir="rescaled", padding=0, step=1,
    save_separate=True
)

supplementary_plot_helpers.stitch_before_after_rosetta(
    pre_rosetta_dir, post_rosetta_dir, rosetta_tiling, run_name,
    [21], "CD20", pre_rosetta_subdir="", post_rosetta_subdir="rescaled", padding=0, step=1,
    save_separate=True
)
supplementary_plot_helpers.stitch_before_after_rosetta(
    pre_rosetta_dir, post_rosetta_dir, rosetta_tiling, run_name,
    [22], "CD4", pre_rosetta_subdir="", post_rosetta_subdir="rescaled", padding=0, step=1,
    save_separate=True
)
supplementary_plot_helpers.stitch_before_after_rosetta(
    pre_rosetta_dir, post_rosetta_dir, rosetta_tiling, run_name,
    [14], "CD31", pre_rosetta_subdir="", post_rosetta_subdir="rescaled", padding=0, step=1,
    save_separate=True
)
supplementary_plot_helpers.stitch_before_after_rosetta(
    pre_rosetta_dir, post_rosetta_dir, rosetta_tiling, run_name,
    [13], "CD8", pre_rosetta_subdir="",  post_rosetta_subdir="rescaled", padding=0, step=1,
    save_separate=True
)

supplementary_plot_helpers.stitch_before_after_rosetta(
    pre_rosetta_dir, post_rosetta_dir, rosetta_tiling, run_name,
    [19], "CD68", pre_rosetta_subdir="",  post_rosetta_subdir="rescaled", padding=0, step=1,
    save_separate=True
)
