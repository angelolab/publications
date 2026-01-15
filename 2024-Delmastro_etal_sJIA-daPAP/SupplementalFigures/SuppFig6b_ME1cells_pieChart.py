#Install Libraries
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
    
    # Add disease category to cell_table
    cell_table = pd.merge(cell_table, fov_labels[["point","FOV_Category"]], on="point", how="left")
    
    # Define Color Map for Cell Phenotypes
    labels = color_map.Pheno
    mycolors = color_map.Hex

    # Keep only PAP involved points
    cats = ['sJIA-daPAP','non-sJIA-PAP']
    cell_table = cell_table[cell_table['FOV_Category'].isin(cats)]
       
    # Keep only ME1 Cells
    cell_table = cell_table[cell_table['topic'] == 1]
    
    #Calculate proportion of cells associated with a given cell phenotype
    values = []
    for cell in labels:
        cell_topic = len(cell_table[cell_table["name"] == cell]) / len(cell_table)
        values.append(cell_topic)
        
    plot_df = pd.DataFrame()
    plot_df["name"] = labels
    plot_df["value"] = values
    plot_df["colors"] = mycolors
    
    #Sort Proportions
    plot_df = plot_df.sort_values(by=['value'])
    
    #Plot Pie Chart
    fig1, ax = plt.subplots()
    ax.pie(plot_df["value"], colors = plot_df["colors"], counterclock=False, startangle = 90, wedgeprops = {"edgecolor" : "white",
                      'linewidth': 2,
                      'antialiased': True})
    my_circle=plt.Circle( (0,0), 0.5, color='white')
    p=plt.gcf()
    p.gca().add_artist(my_circle)
    plt.savefig(path + '/Sub-panels/SupplementalFigures/SuppFig6b_ME1cells_pieChart.pdf', dpi = 1000, transparent = True)
    
if __name__ == '__main__':
     main()