# Install Libraries
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from PIL import ImageColor
import seaborn as sns
import umap.umap_ as umap
import math

def main():
    ## USER INPUT:
    cell = 'M2_Mac'
    chosen_markers = ['IFNg', 'HLA.DR', 'CD45RO', 'TIM3', 'pS6', 'CD57', 'MMP9', 'Ki67', 'GrzB', 'HO.1', 'IDO', 'iNOS', 'H3K9Ac', 'H3K27me3']

    # Set paths
    path = '/Users/aleadelmastro/Library/CloudStorage/GoogleDrive-alead@stanford.edu/.shortcut-targets-by-id/1GXGQlL1W-w8iJENEwcrx7NBGohHlnlCY/PAP/PAP_MIBI/Cohort/MIBI_data/all_points'
    cell_table = pd.read_csv(path + "/celltable_05032024.csv")
    fov_labels = pd.read_csv(path + "/fov_labels.csv")
    color_map = pd.read_csv(path + "/cellpheno_num_colorkey.csv")
    
    # Add disease category to cell_table
    cell_table = pd.merge(cell_table, fov_labels[["point","FOV_Category"]], on="point", how="left")
    
    # Define Color Map for Disease categories
    mycolors = {"sJIA-daPAP": "#1f77b4", "non-sJIA-PAP": "#ff7f0e", "uninvolved": "#2ca02c"} 
    
    #Define markers
    markers = ["CD11c", "CD14", "CD16", "CD163", "CD20", "CD206", "CD209", "CD3", "CD31", "CD4", "CD45", "CD45RO",
                  "CD57", "CD68", "CD8", "Calprotectin", "EPOX", "Foxp3", "GrzB", "HH3", "H3K27me3", "H3K9Ac", "HLA.DR",
                  "HO.1", "IDO", "IFNg", "Ki67", "MMP9", "NaKATPase", "PanCK", "SMA", "TIM3", "Tryptase", "VIM", "iNOS",
                  "pS6"]

    #Subset Data to only include PAP involved and uninvolved cells of interest
    cell_table_sub = cell_table[(cell_table['FOV_Category'] != "non-PAP") & (cell_table['name'] == cell)]
    cell_table_sub_umap = cell_table[(cell_table['FOV_Category'] != "non-PAP") & (cell_table['name'] == cell)][markers]

    #Perform UMAP on subsetted data
    u = umap.UMAP(n_neighbors=15, min_dist=0.5).fit_transform(cell_table_sub_umap)
    
    #Plot UMAP for all cells of interest
    plt.figure(figsize=(4, 3))
    plt.tick_params(axis='x', which='both', bottom=False,
                    top=False, labelbottom=False)
    plt.tick_params(axis='y', which='both', right=False,
                    left=False, labelleft=False)
    plt.scatter(u[:, 0], u[:, 1], c=cell_table_sub['FOV_Category'].map(mycolors))
    for pos in ['right', 'top', 'bottom', 'left']:
        plt.gca().spines[pos].set_visible(False)
    plt.savefig(path + '/Sub-panels/Figure3/' + cell + '/CategoryUMAP.png',dpi=1000, transparent = True)
    plt.show()

    #Plot UMAP with color expression of each marker per cell
    for marker in chosen_markers:
        marker_expr = cell_table_sub[marker].values
        scale_expr = (marker_expr - min(marker_expr))/(max(marker_expr) - min(marker_expr))
        plt.figure(figsize=(4, 3))
        plt.scatter(u[:, 0], u[:, 1], c=scale_expr)
        plt.tick_params(axis='x', which='both', bottom=False,
                    top=False, labelbottom=False)
        plt.tick_params(axis='y', which='both', right=False,
                    left=False, labelleft=False)
        for pos in ['right', 'top', 'bottom', 'left']:
            plt.gca().spines[pos].set_visible(False)
        plt.savefig(path + '/Sub-panels/Figure3/' + cell + '/UMAP_' + marker + '.png',dpi=1000, transparent = True)
        plt.show()

    #Plot UMAP colored by disease category (outline) and lipoprotein status (fill)
    colors = {1: 'white', 0: 'black'}
    edgecolors = mycolors
    plt.figure(figsize=(4, 3))
    plt.scatter(u[:, 0], u[:, 1], c=cell_table_sub['lipoprotein'].map(colors), edgecolors = cell_table_sub['FOV_Category'].map(edgecolors))
    plt.tick_params(axis='x', which='both', bottom=False,
                    top=False, labelbottom=False)
    plt.tick_params(axis='y', which='both', right=False,
                    left=False, labelleft=False)
    for pos in ['right', 'top', 'bottom', 'left']:
        plt.gca().spines[pos].set_visible(False)
    plt.savefig(path + '/Sub-panels/Figure3/' + cell + '/LipoUMAP.png',dpi=1000, transparent = True)
    plt.show()

    # Save UMAP embeddings for future use
    umap_df = pd.DataFrame(u, columns = ['UMAP1', 'UMAP2'])
    umap_df.to_csv(path + '/Sub-panels/Figure3/' + cell + '/umapEmbeddings.csv')

if __name__ == '__main__':
    main()