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
    
    # Functional Markers of Interest
    fun_markers = ["CD45RO","CD57","GrzB","H3K27me3","H3K9Ac","HLA.DR","HO.1","IDO","IFNg","iNOS","Ki67","MMP9","pS6", "TIM3","topic"]
    
    # Find average expression of each marker by ME
    cell_table = cell_table[fun_markers].groupby(["topic"]).mean().reset_index()
    cell_table = cell_table.set_index("topic")
    
    # Scale values
    for col in cell_table.columns:
        cell_table[col] = (cell_table[col].values - np.mean(cell_table[col].values)) / np.std(cell_table[col].values)
    
    # Plot heatmap
    fig, ax = plt.subplots()
    hex_codes = pd.read_csv(path + "/heatmap_colorPalette.csv") #From R colorspace package: Blue-Red 3
    cmap = colors.ListedColormap(hex_codes.x)
    ax = sns.heatmap(cell_table, square=True, cmap=cmap,  linewidths=1,
                linecolor='grey', vmin = -2, vmax = 2)
    plt.tight_layout()
    plt.savefig(path + '/Sub-panels/Figure4/Fig4b_ME_heatmap_functional.pdf', dpi = 1000, transparent = True)
    
if __name__ == '__main__':
     main()