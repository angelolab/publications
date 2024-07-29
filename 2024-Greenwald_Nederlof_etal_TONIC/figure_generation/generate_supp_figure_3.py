import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

import os

import python_files.supplementary_plot_helpers as supplementary_plot_helpers

BASE_DIR = "/Volumes/Shared/Noah Greenwald/TONIC_Cohort/"
raw_dir = "/Volumes/Shared/Noah Greenwald/TONIC_Acquisition/"
SUPPLEMENTARY_FIG_DIR = os.path.join(BASE_DIR, "supplementary_figs")

# show a run with images pre- and post-Rosetta
rosetta_tiling = os.path.join(SUPPLEMENTARY_FIG_DIR, "supp_figure_3b_tiles")
if not os.path.exists(rosetta_tiling):
    os.makedirs(rosetta_tiling)

run_name = "2022-01-14_TONIC_TMA2_run1"
pre_rosetta_dir = os.path.join(raw_dir, "extracted")
post_rosetta_dir = os.path.join(raw_dir, "rosetta")

# NOTE: images not scaled up programmatically, this happens manually in Photoshop
supplementary_plot_helpers.stitch_before_after_rosetta(
    pre_rosetta_dir, post_rosetta_dir, rosetta_tiling, run_name,
    [11], "CD4", post_rosetta_subdir="normalized", padding=0, step=1,
    save_separate=True
)
supplementary_plot_helpers.stitch_before_after_rosetta(
    pre_rosetta_dir, post_rosetta_dir, rosetta_tiling, run_name,
    [18], "CD56", post_rosetta_subdir="normalized", padding=0, step=1,
    save_separate=True
)
supplementary_plot_helpers.stitch_before_after_rosetta(
    pre_rosetta_dir, post_rosetta_dir, rosetta_tiling, run_name,
    [45], "CD31", post_rosetta_subdir="normalized", padding=0, step=1,
    save_separate=True
)
supplementary_plot_helpers.stitch_before_after_rosetta(
    pre_rosetta_dir, post_rosetta_dir, rosetta_tiling, run_name,
    [30], "CD8", post_rosetta_subdir="normalized", padding=0, step=1,
    save_separate=True
)