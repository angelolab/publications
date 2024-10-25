# Install Libraries
import numpy as np
import pandas as pd
from PIL import Image
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
    
    # Keep Only PAP-involved and lipoproteinosis+ cells
    cell_table = cell_table[cell_table['FOV_Category'].isin(['sJIA-daPAP','non-sJIA-PAP'])]
    cell_table = cell_table[cell_table['lipoprotein'] == 1]

    lipo_analysis = cell_table[["topic","lipoprotein"]]
    count = lipo_analysis.groupby(['topic']).size().reset_index()

    # Count number of cells in each ME that are in lipoproteinosis regions
    count = count.sort_values(by=["topic"]).reset_index()
    plot_new = count
    plot_new = plot_new[0].values / np.sum(plot_new[0].values)

    # Generate Bar Plot
    plot = pd.DataFrame()
    labels = [1,2,3,4,5,6,7,8,9]
    plot["labels"] = color_key.Topic.values
    plot["values"] = plot_new * 100 #For Percent
    plot["mycolors"] = color_key.Hex.values
    plot = plot.sort_values('values', ascending=False)
    
    fig1, ax = plt.subplots()
    ax.bar("labels", "values",width = 0.7, color = "mycolors",  data = plot, edgecolor = "black")
    ax.set_ylabel("% of All Lipoproteinosis+ Cells")
    ax.set_xlabel("Microenvironments")
    ax.set_xticks(np.arange(min(labels), max(labels)+1, 1.0))
    plt.savefig(path + '/Sub-panels/Figure4/Fig4g_LipoCells_ME_Distribution.png', dpi = 1000, transparent = True)

if __name__ == '__main__':
     main()