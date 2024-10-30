# Install Libraries
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import colors
import seaborn as sns

def main():
    # Set paths
    path = '/Users/aleadelmastro/Library/CloudStorage/GoogleDrive-alead@stanford.edu/.shortcut-targets-by-id/1GXGQlL1W-w8iJENEwcrx7NBGohHlnlCY/PAP/PAP_MIBI/Cohort/MIBI_data/all_points'
    cell_table = pd.read_csv(path + "/celltable_05032024.csv")
    fov_labels = pd.read_csv(path + "/fov_labels.csv")
    color_key = pd.read_csv(path + "/neighborhood_colorkey.csv")
    
    # Add disease category to cell_table
    cell_table = pd.merge(cell_table, fov_labels[["point","FOV_Category"]], on="point", how="left")
    
    # Functional Markers
    markers = ['CD45RO','IFNg','HLA.DR', 'TIM3','pS6','CD57','MMP9','Ki67','GrzB','HO.1','IDO','iNOS','H3K9Ac','H3K27me3']

    # Keep only lymphocytes within FOVs with a TLS
    points = [12, 14, 38, 43, 44]
    cell_order = ['Bcell',"CD4+_Tcell","CD8+_Tcell","CD57+_CD8+_Tcell","Treg"]
    cell_table = cell_table[cell_table['point'].isin(points)]
    cell_cell = cell_table[cell_table['name'].isin(cell_order)]
    
    # Keep only columns of marker expression and cell phenotype name
    cols = markers + ['name']
    cell_cell = cell_cell[cols]
    
    # Find average marker expression grouped by cell phenotype name
    cell_Ex = cell_cell.groupby('name').mean()
        
    # Plot Heatmap
    cell_Ex = cell_Ex.T
    cell_Ex = cell_Ex[cell_order]
    hex_codes = pd.read_csv(path + "/heatmap_colorPalette.csv") #From R colorspace package: Blue-Red 3
    cmap = colors.ListedColormap(hex_codes.x)
    g = sns.clustermap(cell_Ex, cmap = cmap, vmin = -2, square = True, vmax = 2, linewidths=1,
                linecolor='grey', col_cluster=False, cbar_pos = (1, .2, .03, .4), z_score = 0)
    plt.savefig(path + '/Sub-panels/SupplementalFigures/SuppFig7c_TLS_functionalMarkerHeatmap.pdf', dpi = 1000, transparent = True)
        
if __name__ == '__main__':
     main()