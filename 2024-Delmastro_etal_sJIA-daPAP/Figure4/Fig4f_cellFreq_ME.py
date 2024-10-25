#Install Libraries
import numpy as np
import pandas as pd
from PIL import Image
import matplotlib.pyplot as plt
from matplotlib import colors
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

def main():
    # Set paths
    path = '/Users/aleadelmastro/Library/CloudStorage/GoogleDrive-alead@stanford.edu/.shortcut-targets-by-id/1GXGQlL1W-w8iJENEwcrx7NBGohHlnlCY/PAP/PAP_MIBI/Cohort/MIBI_data/all_points'
    cell_table = pd.read_csv(path + "/celltable_05032024.csv")
    fov_labels = pd.read_csv(path + "/fov_labels.csv")
    color_map = pd.read_csv(path + "/neighborhood_colorkey.csv")
    
    # Add disease category to cell_table
    cell_table = pd.merge(cell_table, fov_labels[["point","FOV_Category"]], on="point", how="left")
    
    # Define Topic of Interest
    topicInterest = 1
    
    # Keep on FOVs from regions of interests and define cell types to compare 
    data = pd.DataFrame(columns = ['point','cell','category','Frequency'])
    cats = ['sJIA-daPAP','non-sJIA-PAP']
    cells = ['Eosinophil', 'Neutrophil', 'M2_Mac']

    # Define Color Map for Microenvironments
    topics = color_map.Topic
    mycolors = color_map.Hex
    
    #Calculate proportion of cells associated with a given ME
    for cat in cats:
        cell_cat = cell_table[cell_table['FOV_Category'] == cat]
        points = cell_cat["point"].unique()
            
        for point in points:
            cell_point = cell_cat[cell_cat['point'] == point]
            cell_topic = cell_point[cell_point["topic"] == topicInterest]
            if len(cell_topic) != 0:
                for cell in cells:
                    if len(cell_point[cell_point["name"] == cell]) != 0:
                        cell_cell = len(cell_topic[cell_topic['name'] == cell])/len(cell_point[cell_point["name"] == cell])
                        data = data.append({'point': point, 'cell': cell, 'category': cat, 'Frequency': cell_cell}, ignore_index=True)

    # Perform statistical comparison between disease categories for each cell phenotype
    for cell in cells:
        data_test = data[data['cell'] == cell]
        data_SJIA = data_test[data_test['category'] == 'sJIA-daPAP']['Frequency'].values
        data_nonSJIA = data_test[data_test['category'] == 'non-sJIA-PAP']['Frequency'].values
        _, p = stats.ttest_ind(data_SJIA, data_nonSJIA)
        print(cell,": ", p)
        
    # Plot box and whisker
    fig1, ax = plt.subplots()
    ax = sns.boxplot(x="cell", y="Frequency", hue="category",
                   data=data, palette="tab10", dodge=True, showfliers = False)
    ax.set_ylim(-0.05,1.13)
    plt.legend([],[], frameon=False)
    plt.rcParams["figure.figsize"] = (5,5)
    plt.savefig(path + '/Sub-panels/Figure4/Fig4f_cellFreq_ME' + str(topicInterest) +'.png', dpi = 1000, transparent = True)
    
if __name__ == '__main__':
     main()