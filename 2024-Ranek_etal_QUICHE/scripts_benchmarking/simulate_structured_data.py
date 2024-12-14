
import os
import pandas as pd
import numpy as np
import quiche as qu
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Rectangle
import random
import os
import anndata
from sketchKH import *
from scipy.spatial.distance import cdist
import scanpy as sc
import shutil
import gc

base_dir = 'data/tnbc_spain'
annotations_by_mask = pd.read_csv(os.path.join(base_dir, 'cell_annotation_mask.csv'))
cell_table = pd.read_csv(os.path.join(base_dir, 'cell_table_size_normalized_samples_cell_labels_updated.csv'))
compartmentalized_tumors = ['TMA31_R3C1', 'TMA31_R3C9', 'TMA41_R4C4', 'TMA31_R4C5', 'TMA31_R5C4', 'TMA31_R5C5', 'TMA31_R7C1', 'TMA32_R5C7', 'TMA32_R8C5', 'TMA32_R10C4', 'TMA33_R5C8',
                            'TMA33_R8C4', 'TMA33_R9C4', 'TMA33_R10C5', 'TMA33_R12C2', 'TMA34_R4C2', 'TMA34_R9C8', 'TMA34_R12C3', 'TMA35_R3C2', 'TMA35_R4C3', 'TMA36_R2C7',
                            'TMA36_R9C9', 'TMA37_R3C1', 'TMA37_R4C4', 'TMA37_R4C7', 'TMA37_R7C4', 'TMA37_R10C5', 'TMA38_R5C2', 'TMA39_R5C6', 'TMA39_R1C1', 'TMA39_R2C4', 'TMA39_R3C4',
                            'TMA39_R5C4', 'TMA39_R5C6', 'TMA39_R5C8', 'TMA39_R6C1', 'TMA39_R9C2', 'TMA39_R9C6', 'TMA40_R4C7', 'TMA40_R5C2', 'TMA40_R6C3', 'TMA40_R6C6', 'TMA40_R7C6',
                            'TMA40_R7C7', 'TMA40_R8C6', 'TMA40_R10C7', 'TMA41_R1C3', 'TMA41_R2C3', 'TMA41_R4C2', 'TMA41_R4C3', 'TMA41_R4C4', 'TMA42_R2C2', 'TMA42_R3C5', 'TMA42_R4C1',
                            'TMA42_R6C1', 'TMA42_R6C5', 'TMA42_R7C4', 'TMA43_R1C3', 'TMA43_R3C3', 'TMA43_R5C7', 'TMA43_R8C7', 'TMA43_R9C8', 'TMA43_R11C5', 'TMA44_R3C3', 'TMA44_R3C7',
                            'TMA44_R7C2', 'TMA44_R7C6', 'TMA44_R8C1', 'TMA44_R8C3', 'TMA44_R9C5', 'TMA44_R10C6', 'TMA44_R12C2', 'TMA44_R12C7', 'TMA44_R13C7', 'TMA44_R14C7']
compartment_colormap = pd.DataFrame({'mask_name': ['cancer_core', 'cancer_border', 'stroma', 'immune1'], 
                                     'color': ['blue', 'deepskyblue','firebrick', 'orange']})
annotations_by_mask_merged = pd.merge(annotations_by_mask, cell_table.loc[:, ['fov', 'label', 'centroid-1', 'centroid-0']], on = ['fov', 'label'])
annotations_by_mask_merged['mask_name'].replace(['stroma_core', 'stroma_border'], ['stroma', 'stroma'], inplace = True)
annotations_by_mask_merged = annotations_by_mask_merged[np.isin(annotations_by_mask_merged.fov, compartmentalized_tumors)]
adata_expression = anndata.read_h5ad(os.path.join(base_dir, 'adata_simulated_expression_groups_large.h5ad'))
pct_change_list = [0.05, 0.08, 0.1, 0.2, 0.25]
radius_list = [100, 250, 500]
prevalence_list = [0.2, 0.4, 0.6, 0.8, 1.0]
save_directory = os.path.join('data', 'simulated', 'structured', 'adata')
qu.pp.make_directory(save_directory)
for trial in range(5):
    for pct_change in pct_change_list:
        for radius in radius_list:
            for prevalence in prevalence_list:
                attempt = 0 
                if pct_change in [0.08, 0.2] and radius != 500:
                    continue
                while attempt < 5:
                    try:
                        print(f"Trial: {trial}, Attempt: {attempt}, pct_change: {pct_change}, radius: {radius}, prevalence: {prevalence}")
                        prev = int(prevalence * 100)
                        pct = int(pct_change * 100)
                        param_id = f'prev{prev}_r{radius}_p{pct}_trial{trial}'
                        num_to_change = int(radius * pct_change)
                        count_df = annotations_by_mask_merged.groupby(['fov', 'mask_name']).count()['label'].unstack().copy()
                        ratio_df = pd.Series([400, 600, 800], index=['cancer_border', 'cancer_core', 'stroma'])
                        fov_list = (count_df > ratio_df + num_to_change).sum(axis=1)
                        fov_list = list(fov_list[fov_list == 3].index)
                        ratio_df = pd.Series([400, 600, 800, num_to_change], index=['cancer_border', 'cancer_core', 'stroma', 'immune1'])
                        annotations_by_mask_run = annotations_by_mask.copy()
                        adata = qu.tl.simulate_structured_data(
                            annotations_by_mask_run, cell_table, adata_expression, n_cond=10,
                            fov_key='fov', fov_list=fov_list, labels_key='mask_name',
                            cond1='cancer_core', cond2='cancer_border', radius=radius, p=2,
                            condition_id='immune1', num_to_change=num_to_change,
                            compartment_colormap=compartment_colormap, prevalence=prevalence,
                            ratio_df=ratio_df, cell_types=['cancer_core', 'cancer_border', 'stroma', 'immune1'],
                            sim_cell_types=['Group1', 'Group1', 'Group2', 'Group3'], group_key='group',
                            spatial_key='spatial', n_jobs=8)
                        adata.write_h5ad(os.path.join(save_directory, f'adata_{param_id}.h5ad'))
                        print(f"Completed successfully for Trial: {trial}, pct_change: {pct_change}, radius: {radius}, prevalence: {prevalence}")
                        del adata
                        gc.collect()
                        break
                    except Exception as e:
                        print(f"Error: {e}. Attempt {attempt} failed for Trial: {trial}, pct_change: {pct_change}, radius: {radius}, prevalence: {prevalence}")
                        attempt += 1
                        if attempt >= 5:
                            print(f"Failed after 5 attempts for Trial: {trial}, pct_change: {pct_change}, radius: {radius}, prevalence: {prevalence}")
