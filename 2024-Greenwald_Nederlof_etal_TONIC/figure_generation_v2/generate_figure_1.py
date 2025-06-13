import os
import shutil
import numpy as np

import pandas as pd
import matplotlib

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

import matplotlib.pyplot as plt
import seaborn as sns
import skimage.io as io
from skimage.segmentation import find_boundaries


from ark.utils.data_utils import erode_mask
from ark.utils.plot_utils import cohort_cluster_plot
import ark.settings as settings
from venny4py.venny4py import venny4py


base_dir = '/Volumes/Shared/Noah Greenwald/TONIC_Cohort'
metadata_dir = os.path.join(base_dir, 'intermediate_files/metadata')
plot_dir = os.path.join(base_dir, 'figures')
image_dir = os.path.join(base_dir, 'image_data/samples/')
sequence_dir = os.path.join(base_dir, 'sequencing_data')
segmentation_dir = os.path.join(base_dir, 'segmentation_data/deepcell_output/')

harmonized_metadata = pd.read_csv(os.path.join(metadata_dir, 'harmonized_metadata.csv'))
timepoint_metadata = pd.read_csv(os.path.join(metadata_dir, 'TONIC_data_per_timepoint.csv'))
timepoint_metadata = timepoint_metadata.loc[timepoint_metadata.MIBI_data_generated, :]
timepoint_metadata = timepoint_metadata[timepoint_metadata.Timepoint.isin(['baseline', 'primary', 'pre_nivo', 'on_nivo'])]

outcome_data = pd.read_csv(os.path.join(base_dir, 'intermediate_files/metadata/patient_clinical_data.csv'))
timepoint_metadata = timepoint_metadata.merge(outcome_data, on='Patient_ID')
timepoint_metadata = timepoint_metadata.loc[timepoint_metadata.Clinical_benefit.isin(['Yes', 'No']), :]

# create venn diagram across modalities
wes_metadata = pd.read_csv(os.path.join(sequence_dir, 'preprocessing/TONIC_WES_meta_table.tsv'), sep='\t')
rna_metadata = pd.read_csv(os.path.join(sequence_dir, 'preprocessing/TONIC_tissue_rna_id.tsv'), sep='\t')
rna_metadata = rna_metadata.merge(harmonized_metadata[['Patient_ID', 'Tissue_ID']].drop_duplicates(), on='Tissue_ID', how='left')

# exclude SD patients
wes_metadata = wes_metadata.rename(columns={'Individual.ID': 'Patient_ID'})
wes_metadata = wes_metadata.loc[wes_metadata.Clinical_benefit.isin(['Yes', 'No']), :]
rna_metadata = rna_metadata.merge(outcome_data, on='Patient_ID')
rna_metadata = rna_metadata.loc[rna_metadata.Clinical_benefit.isin(['Yes', 'No']), :]

MIBI_ids = set(timepoint_metadata.Patient_ID.values)
WES_ids = set(wes_metadata.Patient_ID.values)
RNA_ids = set(rna_metadata.Patient_ID.values)

sets = {
    'MIBI': MIBI_ids,
    'WES': WES_ids,
    'RNA': RNA_ids}

venny4py(sets=sets)
plt.savefig(os.path.join(plot_dir, 'figure1b_venn_diagram.pdf'), dpi=300, bbox_inches='tight')
plt.close()

# identify representative images for visualization
study_fovs = harmonized_metadata.loc[harmonized_metadata.Timepoint.isin(['primary', 'baseline', 'pre_nivo', 'on_nivo']), 'fov'].values

# representative images
fov_list = ['TONIC_TMA21_R6C6', 'TONIC_TMA2_R6C6', 'TONIC_TMA18_R4C5', 'TONIC_TMA20_R12C4', 'TONIC_TMA10_R3C1', 'TONIC_TMA10_R6C2']

overlay_fovs = ['TONIC_TMA2_R6C6', 'TONIC_TMA10_R3C1']
overlay_folder = os.path.join(plot_dir, 'Figure1_representative_images')
if not os.path.exists(overlay_folder):
    os.mkdir(overlay_folder)

# copy images to new folder
for fov in overlay_fovs:
    fov_folder = os.path.join(image_dir, fov)
    save_folder = os.path.join(overlay_folder, fov)
    if not os.path.exists(save_folder):
        os.mkdir(save_folder)

    fov_files = os.listdir(fov_folder)
    fov_files = [x for x in fov_files if x.endswith('.tiff')]

    for file in fov_files:
        shutil.copy(os.path.join(fov_folder, file), os.path.join(save_folder, file))

# generate color overlays
cell_table_clusters = pd.read_csv(os.path.join(base_dir, 'analysis_files/cell_table_clusters.csv'))

# scale bars: 800 um, 1/8 = 100um
# crops = 500 pixels, 500/2048 = 0.25, one quarter the size. Scale bar = 1/2 of the crop size


crop_cmap = pd.DataFrame({'cell_cluster_broad': ['Cancer', 'Structural', 'Mono_Mac', 'T', 'Other', 'Granulocyte', 'NK', 'B'],
                          'color': ['dimgrey', 'darksalmon', 'red', 'yellow', 'yellowgreen', 'aqua', 'dodgerblue', 'darkviolet']})

crop_plot_dir = os.path.join(plot_dir, 'Figure1_crop_overlays')
if not os.path.exists(crop_plot_dir):
    os.mkdir(crop_plot_dir)

cell_table_subset = cell_table_clusters.loc[cell_table_clusters.fov.isin(overlay_fovs), :]

cohort_cluster_plot(
    fovs=overlay_fovs,
    seg_dir=segmentation_dir,
    save_dir=crop_plot_dir,
    cell_data=cell_table_subset,
    erode=True,
    fov_col=settings.FOV_ID,
    label_col=settings.CELL_LABEL,
    cluster_col='cell_cluster_broad',
    seg_suffix="_whole_cell.tiff",
    display_fig=False,
    cmap=crop_cmap,
)

# # make segmentation only overlay
# for fov in overlay_fovs:
#     seg_mask = io.imread(os.path.join(segmentation_dir, '{}_whole_cell.tiff'.format(fov)))[0, :, :]
#     bool_mask = find_boundaries(seg_mask, mode='inner')
#
#     greyscale_cells = np.full((seg_mask.shape[0], seg_mask.shape[1]), 255, dtype='uint8')
#
#     greyscale_cells[seg_mask > 0] = 0
#     greyscale_cells[bool_mask] = 160
#
#     io.imsave(os.path.join(crop_plot_dir, '{}_segmentation_overlay.png'.format(fov)), greyscale_cells)

# tumor compartment overlays
annotations_by_mask = pd.read_csv(os.path.join(base_dir, 'intermediate_files/mask_dir/cell_annotation_mask.csv'))

# generate overlays
compartment_colormap = pd.DataFrame({'mask_name': ['cancer_core', 'cancer_border', 'stroma_border', 'stroma_core'],
                         'color': ['blue', 'deepskyblue', 'lightcoral', 'firebrick']})

compartment_plot_dir = os.path.join(plot_dir, 'Figure1_compartment_overlays')
if not os.path.exists(compartment_plot_dir):
    os.mkdir(compartment_plot_dir)


cell_table_subset = annotations_by_mask.loc[(annotations_by_mask.fov.isin(overlay_fovs)), :]

cohort_cluster_plot(
    fovs=overlay_fovs,
    seg_dir=segmentation_dir,
    save_dir=compartment_plot_dir,
    cell_data=cell_table_subset,
    erode=True,
    fov_col=settings.FOV_ID,
    label_col=settings.CELL_LABEL,
    cluster_col='mask_name',
    seg_suffix="_whole_cell.tiff",
    cmap=compartment_colormap,
    display_fig=False,
)


# generate crops
crop_dict = {'TONIC_TMA2_R6C6': [[300, 400], [1300, 800], ['CD45_ECAD_Overlay.tif', 'CD14_CD31_SMA.tif', 'Ki67_CD38_CD8_overlay.tif']]}
            #'TONIC_TMA10_R3C1': [[750, 750], [1300, 1100], 'CD8_CD45RO_HLADR.tif']}

for fov, crop_info in crop_dict.items():
    save_folder = os.path.join(plot_dir, 'Figure1_representative_images/{}'.format(fov))
    crop_coords_1, crop_coords_2, files = crop_info

    # crop image
    for crop_coords in [crop_coords_1, crop_coords_2]:
        # read in segmentation mask and crop to correct size for visualization
        seg_mask = io.imread(os.path.join(segmentation_dir, '{}_whole_cell.tiff'.format(fov)))
        seg_crop = seg_mask[0, crop_coords[0]:crop_coords[0] + 500, crop_coords[1]:crop_coords[1] + 500]
        seg_crop = erode_mask(seg_crop, connectivity=2, mode="thick", background=0)
        seg_crop[seg_crop > 0] = 255
        seg_crop = seg_crop.astype('uint8')
        io.imsave(os.path.join(save_folder, 'crop_{}_mask.png'.format(crop_coords[0])), seg_crop)

        # same for overlay mask
        overlay_mask = io.imread(os.path.join(crop_plot_dir, 'cluster_masks_colored/{}.tiff'.format(fov)))
        overlay_mask = overlay_mask[crop_coords[0]:crop_coords[0] + 500, crop_coords[1]:crop_coords[1] + 500, :]
        io.imsave(os.path.join(crop_plot_dir, 'cluster_masks_colored/{}_{}_crop.tiff'.format(fov, crop_coords[0])), overlay_mask)

        # and compartment mask
        compartment_mask = io.imread(os.path.join(compartment_plot_dir, 'cluster_masks_colored/{}.tiff'.format(fov)))
        compartment_mask = compartment_mask[crop_coords[0]:crop_coords[0] + 500, crop_coords[1]:crop_coords[1] + 500, :]
        io.imsave(os.path.join(compartment_plot_dir, 'cluster_masks_colored/{}_{}_crop.tiff'.format(fov, crop_coords[0])), compartment_mask)



# cell cluster heatmap
study_fovs = np.delete(study_fovs, np.where(study_fovs == 'TONIC_TMA14_R1C1'))

# Markers to include in the heatmap
markers = ["ECAD", "CK17", "CD45", "CD3", "CD4", "CD8", "FOXP3", "CD20", "CD56", "CD14", "CD68",
           "CD163", "CD11c", "HLADR", "ChyTr", "Calprotectin", "FAP",  "Fibronectin",
           "Collagen1", "Vim", "SMA", "CD31"]

cell_ordering = ['Cancer_1', 'Cancer_2', 'Cancer_3', 'CD4T', 'CD8T', 'Treg', 'T_Other', 'B',
                 'NK', 'CD68_Mac', 'CD163_Mac', 'Mac_Other', 'Monocyte', 'APC','Mast', 'Neutrophil',
                 'CAF', 'Fibroblast', 'Smooth_Muscle','Endothelium']

# # Get average across each cell phenotype
# cell_counts = pd.read_csv(os.path.join(base_dir, "analysis_files/cell_table_counts.csv"))
# phenotype_col_name = "cell_cluster"
# cell_counts = cell_counts.loc[cell_counts.fov.isin(study_fovs), :]
# mean_counts = cell_counts.groupby(phenotype_col_name)[markers].mean()
# mean_counts.to_csv(os.path.join(plot_dir, "figure1_cell_cluster_marker_means.csv"))

# read previously generated averages
mean_counts = pd.read_csv(os.path.join(plot_dir, "figure1_cell_cluster_marker_means.csv"))
mean_counts = mean_counts.set_index('cell_cluster')
mean_counts = mean_counts.reindex(cell_ordering)

# set column order
mean_counts = mean_counts[markers]


# functional marker expression per cell type
core_df_func = pd.read_csv(os.path.join(base_dir, 'output_files/functional_df_per_core_filtered_plot.csv'))
plot_df = core_df_func.loc[core_df_func.fov.isin(study_fovs), :]
plot_df = plot_df.loc[plot_df.metric == 'cluster_freq', :]
plot_df = plot_df.loc[plot_df.subset == 'all', :]
plot_df = plot_df.loc[~plot_df.functional_marker.isin(['Vim', 'CD45RO_CD45RB_ratio', 'H3K9ac_H3K27me3_ratio', 'HLA1'])]

# average across cell types
plot_df = plot_df[['cell_type', 'functional_marker', 'value']].groupby(['cell_type', 'functional_marker']).mean().reset_index()
plot_df = pd.pivot(plot_df, index='cell_type', columns='functional_marker', values='value')

# set index based on cell_ordering
plot_df = plot_df.reindex(cell_ordering)

cols = ['Ki67','GLUT1','CD45RO','CD69', 'PD1','CD57','TBET', 'TCF1',
        'CD45RB', 'TIM3', 'PDL1','Fe','IDO','CD38']
plot_df = plot_df[cols]

# combine together
combined_df = pd.concat([mean_counts, plot_df], axis=1)

# plot heatmap
sns.clustermap(combined_df, z_score=1, cmap="vlag", center=0, vmin=-3, vmax=3, xticklabels=True, yticklabels=True,
                    row_cluster=False,col_cluster=False)
plt.tight_layout()
plt.savefig(os.path.join(plot_dir, 'figure1f_combined_heatmap.pdf'))
plt.close()