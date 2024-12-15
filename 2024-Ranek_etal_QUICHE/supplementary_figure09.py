import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import os
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, Normalize
import matplotlib.gridspec as gridspec
import seaborn as sns
import numpy as np
import shutil
from alpineer import io_utils
import supplementary_plot_helpers

BASE_DIR = "/Volumes/Shared/Noah Greenwald/TNBC_Cohorts/SPAIN"
raw_dir = "/Volumes/Shared/Noah Greenwald/TNBC_Acquisition/SPAIN"
SUPPLEMENTARY_FIG_DIR = os.path.join("publications", "supplementary_figures", "supplementary_figure09")

# show a run with images stitched in acquisition order pre- and post-normalization
norm_tiling = os.path.join(SUPPLEMENTARY_FIG_DIR, "supp_figure_09_tiles")
if not os.path.exists(norm_tiling):
    os.makedirs(norm_tiling)

run_name = "2023-04-13_SPAIN_TMA32_run1"
pre_norm_dir = os.path.join(raw_dir, "rosetta")
post_norm_dir = os.path.join(raw_dir, "normalized")

# NOTE: images not scaled up programmatically, this happens manually in Photoshop
supplementary_plot_helpers.stitch_before_after_norm(
    pre_norm_dir, post_norm_dir, norm_tiling, run_name,
    [11, 12, 13, 14, 15, 17, 18, 20, 22, 23, 24, 28, 29, 30, 31, 32, 33, 34, 35,
    36, 39, 40, 41, 42, 43, 44, 45, 46, 47], "H3K9ac", pre_norm_subdir="rescaled", num_rows = 7, padding=0, step=1)

# Generate QC plots showing signal intensity over control TMAs and within each TMA
qc_metrics = ["Non-zero mean intensity"]

channel_exclude = ["Au", "CD11c_nuc_exclude", "CK17_smoothed", "ECAD_smoothed", "FOXP3_nuc_include",
                 "Noodle", "chan_39", "chan_45", "chan_48", "chan_115", "chan_141", "Calprotectin_old"]

# FOV spatial location
cohort_path = os.path.join(BASE_DIR, "image_data", "samples")
qc_tma_metrics_dir = os.path.join(raw_dir, "qc_metrics", "qc_tma_metrics")
if not os.path.exists(qc_tma_metrics_dir):
    os.makedirs(qc_tma_metrics_dir)

fovs = io_utils.list_folders(cohort_path)
tmas = list(set([fov.split('_R')[0] for fov in fovs]))

qc_tmas = supplementary_plot_helpers.QCTMA(
    qc_metrics=qc_metrics,
    cohort_path=cohort_path,
    metrics_dir=qc_tma_metrics_dir,
)

qc_tmas.compute_qc_tma_metrics(tmas=tmas)
qc_tmas.qc_tma_metrics_zscore(tmas=tmas, channel_exclude=channel_exclude)

supplementary_plot_helpers.qc_tmas_metrics_plot(qc_tmas=qc_tmas, tmas=tmas, save_figure=True, dpi=300)

shutil.copy(os.path.join(qc_tma_metrics_dir, "figures", "cross_TMA_averages_nonzero_mean_stats.pdf"),
            os.path.join(SUPPLEMENTARY_FIG_DIR, "supp_figure_09c.pdf"))

# longitudinal controls
control_path = os.path.join(BASE_DIR, "image_data", "controls")
qc_control_metrics_dir = os.path.join(raw_dir, "qc_metrics", "qc_longitudinal_control")
if not os.path.exists(qc_control_metrics_dir):
    os.makedirs(qc_control_metrics_dir)

folders = io_utils.list_folders(control_path, "TMA34_")
control_substrs = [name.split("_")[1] + '_' + name.split("_")[2] if len(name.split("_")) == 3
                   else name.split("_")[1] + '_' + name.split("_")[2]+'_' + name.split("_")[3]
                   for name in folders]

# loop over each control sample
all_folders = io_utils.list_folders(control_path)
for i, control in enumerate(control_substrs):
    control_sample_name = control
    fovs = [folder for folder in all_folders if control.lower() in folder.lower()]
    if 'duplicate1' not in control:
        fovs = [i for i in fovs if 'duplicate1' not in i]
    else:
        fovs = [i for i in fovs if 'duplicate1' in i]
        
    qc_control = supplementary_plot_helpers.QCControlMetrics(
        qc_metrics=qc_metrics,
        cohort_path=control_path,
        metrics_dir=qc_control_metrics_dir,
    )

    qc_control.compute_control_qc_metrics(
        control_sample_name=control_sample_name,
        fovs=fovs,
        channel_exclude=channel_exclude,
        channel_include=None,
    )

    supplementary_plot_helpers.longitudinal_control_heatmap(
        qc_control=qc_control, control_sample_name=control_sample_name, save_figure=True, dpi=300
    )
# aggregate data from each control sample
dfs = []
for control in control_substrs:
    df = pd.read_csv(os.path.join(qc_control_metrics_dir, f"{control}_combined_nonzero_mean_stats.csv"))
    df['fov'] = [i.replace('_' + control.lower(), '') for i in list(df['fov'].str.lower())]
    log2_norm_df: pd.DataFrame = df.pivot(
        index="channel", columns="fov", values="Non-zero mean intensity"
    ).transform(func=lambda row: np.log2(row / row.mean()), axis=1)
    
    if control != 'Tonsil_Bottom_duplicate1':
        dup_col = [col for col in log2_norm_df.columns if 'duplicate1' in col]
        log2_norm_df = log2_norm_df.drop(columns=dup_col) if dup_col else log2_norm_df

    mean_t_df: pd.DataFrame = (
        log2_norm_df.mean(axis=0)
        .to_frame(name="mean")
        .transpose()
        .sort_values(by="mean", axis=1)
    )
    transformed_df: pd.DataFrame = pd.concat(
        objs=[log2_norm_df, mean_t_df]
    ).sort_values(by="mean", axis=1, inplace=False)
    transformed_df.rename_axis("channel", axis=0, inplace=True)
    transformed_df.rename_axis("fov", axis=1, inplace=True)

    dfs.append(transformed_df)
    
all_data = pd.concat(dfs).replace([np.inf, -np.inf], 0).groupby(['channel']).mean().sort_values(by="mean", axis=1, inplace=False).round(2)
all_data.columns = all_data.columns.str.upper()

# generate heatmap
fig = plt.figure(figsize=(10,12), dpi=300)
fig.set_layout_engine(layout="constrained")
gs = gridspec.GridSpec(nrows=2, ncols=1, figure=fig, height_ratios=[len(all_data.index) - 1, 1])
_norm = Normalize(vmin=-1, vmax=1)
_cmap = sns.color_palette("vlag", as_cmap=True)
fig.suptitle(f"Average per TMA - QC: Non-zero Mean Intensity")

annotation_kws = {
    "horizontalalignment": "center",
    "verticalalignment": "center",
    "fontsize": 8,
}

ax_heatmap = fig.add_subplot(gs[0, 0])
sns.heatmap(
    data=all_data[~all_data.index.isin(["mean"])],
    ax=ax_heatmap,
    linewidths=1,
    linecolor="black",
    cbar_kws={"shrink": 0.5},
    annot=True,
    annot_kws=annotation_kws,
    xticklabels=False,
    norm=_norm,
    cmap=_cmap,
)

ax_heatmap.collections[0].colorbar.ax.set_title(r"$\log_2(QC)$")
ax_heatmap.set_yticks(
    ticks=ax_heatmap.get_yticks(),
    labels=ax_heatmap.get_yticklabels(),
    rotation=0,
)
ax_heatmap.set_xlabel(None)

ax_avg = fig.add_subplot(gs[1, 0])
sns.heatmap(
    data=all_data[all_data.index.isin(["mean"])],
    ax=ax_avg,
    linewidths=1,
    linecolor="black",
    annot=True,
    annot_kws=annotation_kws,
    fmt=".2f",
    cmap=ListedColormap(["white"]),
    cbar=False,
)
ax_avg.set_yticks(
    ticks=ax_avg.get_yticks(),
    labels=["Mean"],
    rotation=0,
)
ax_avg.set_xticks(
    ticks=ax_avg.get_xticks(),
    labels=ax_avg.get_xticklabels(),
    rotation=45,
    ha="right",
    rotation_mode="anchor",
)
ax_heatmap.set_ylabel("Channel")
ax_avg.set_xlabel("FOV")

fig.savefig(fname=os.path.join(qc_control_metrics_dir, "figures", "log2_avgs.pdf"), dpi=300,
            bbox_inches="tight")
shutil.copy(os.path.join(qc_control_metrics_dir, "figures", "log2_avgs.pdf"),
            os.path.join(SUPPLEMENTARY_FIG_DIR, "supp_figure_09d.pdf"))