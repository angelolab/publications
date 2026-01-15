# Install Libraries
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import colors
import seaborn as sns
from scipy import stats

def main():
    # Set paths
    path = '/Users/aleadelmastro/Library/CloudStorage/GoogleDrive-alead@stanford.edu/.shortcut-targets-by-id/1GXGQlL1W-w8iJENEwcrx7NBGohHlnlCY/PAP/PAP_MIBI/Cohort/MIBI_data/all_points'
    cell_table = pd.read_csv(path + "/celltable_05032024.csv")
    fov_labels = pd.read_csv(path + "/fov_labels.csv")
    color_map = pd.read_csv(path + "/cellpheno_num_colorkey.csv")
    
    # Define Color Map for Cell Phenotypes
    labels = color_map.Pheno
    mycolors = color_map.Hex
    cmap = dict(zip(labels, mycolors))
    
    # Add disease category to cell_table
    cell_table = pd.merge(cell_table, fov_labels[["point","FOV_Category"]], on="point", how="left")

    #List of ME numbers
    topics = [2,4,5]
    data = pd.DataFrame(columns = ['point','cell','Frequency'])
    cells = ['CD4+_Tcell','Bcell', 'CD8+_Tcell','Treg', 'CD57+_CD8+_Tcell']
    points = [12, 14, 38, 43, 44] # Points with visible TLS
    
    #Calculate proportion of cells associated with listed MEs above
    cell_ME = cell_table[cell_table['topic'].isin(topics)]
    for point in points:
        cell_point = cell_ME[cell_ME['point'] == point]
        cell_lymph = cell_point[cell_point['name'].isin(cells)]
        for cell in cells:
            cell_cell = len(cell_lymph[cell_lymph['name'] == cell])/len(cell_lymph)
            data = data.append({'point': point, 'cell': cell, 'Frequency': cell_cell}, ignore_index=True)
        
    # Plot Bar Graph
    fig1, ax = plt.subplots()
    ax = sns.barplot(x="cell", y="Frequency",
                   data=data, palette=cmap, capsize=.2)

    ax.set_ylim(0,0.7)
    plt.legend([],[], frameon=False)
    plt.rcParams["figure.figsize"] = (5,5)
    plt.savefig(path + '/Sub-panels/SupplementalFigures/SuppFig7b_TLS_frequency.png', dpi = 1000)
    
if __name__ == '__main__':
     main()