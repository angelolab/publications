#Install Libraries
import numpy as np
import pandas as pd
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

    # Define Color Map for Microenvironments
    topics = color_map.Topic
    mycolors = color_map.Hex
            
    data = pd.DataFrame(columns = ['point','topic','FOV_Category','Frequency'])
    cats = ['sJIA-daPAP','non-sJIA-PAP','uninvolved']
    
    #Calculate proportion of cells associated with a given ME
    for cat in cats:
        cell_cat = cell_table[cell_table['FOV_Category'] == cat]
        points = cell_cat.point.unique()
            
        for point in points:
            cell_point = cell_cat[cell_cat['point'] == point]
            for topic in topics:
                cell_topic = len(cell_point[cell_point["topic"] == topic]) / len(cell_point)
                data = data.append({'point': point, 'topic': topic, 'FOV_Category': cat, 'Frequency': cell_topic}, ignore_index=True)
    
    #Perform Statistical Comparisons between proportions by group and ME
    for topic in topics:
        data_test = data[data['topic'] == topic]
        data_SJIA = data_test[data_test['FOV_Category'] == 'sJIA-daPAP']['Frequency'].values
        data_nonSJIA = data_test[data_test['FOV_Category'] == 'non-sJIA-PAP']['Frequency'].values
        data_un = data_test[data_test['FOV_Category'] == 'uninvolved']['Frequency'].values
        _, p = stats.ttest_ind(data_SJIA, data_nonSJIA)
        print("SJIA-daPAP vs. non-sJIA-PAP - Topic ",topic,": ", p)
        _, p = stats.ttest_ind(data_SJIA, data_un)
        print("SJIA-daPAP vs. Uninvolved - Topic ",topic,": ", p)
        _, p = stats.ttest_ind(data_nonSJIA, data_un)
        print("non-sJIA-PAP vs. Uninvolved - Topic ",topic,": ", p)
        
    # Generate Grouped Box & Whisker Plots
    fig1, ax = plt.subplots()
    plt.rcParams["figure.figsize"] = (10,5)
    ax = sns.boxplot(x="topic", y="Frequency", hue="FOV_Category",
                   data=data, palette="tab10", dodge=True, showfliers = False)
    ax.set_ylim(0,1.0)
    plt.legend([],[], frameon=False)
    plt.savefig(path + '/Sub-panels/Figure4/Fig4d_Frequency_ME_byPAPCat.png', dpi = 1000, transparent = True)
    
if __name__ == '__main__':
     main()