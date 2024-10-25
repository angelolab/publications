# Install Libraries
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import colors
import seaborn as sns
import statistics

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
    cell_table_PAP = cell_table[cell_table['FOV_Category'].isin(["sJIA-daPAP","non-sJIA-PAP"])]
    points = cell_table_PAP.point.unique()
    points.sort()

    # Differentiate lipoprotein+ and Lipoprotein- cells and groupby cell phenotype
    pos_cells = cell_table_PAP[cell_table_PAP["lipoprotein"] == 1]
    neg_cells = cell_table_PAP[cell_table_PAP["lipoprotein"] == 0]
    pos_cells = pos_cells.groupby(["point","pixelfreq_hclust_cap"]).size().reset_index()
    neg_cells = neg_cells.groupby(["point","pixelfreq_hclust_cap"]).size().reset_index()
    datas = pd.DataFrame()
    
    # Calculate enrichment of Lipoprotein+/Lipoprotein- per cell phenotype
    for i in points:
        pos = pos_cells[pos_cells['point'] == i]
        neg = neg_cells[neg_cells['point'] == i]
        pos = pos[["pixelfreq_hclust_cap",0]]
        neg = neg[["pixelfreq_hclust_cap",0]]
        
        for j in range(1,21):
            if j not in pos["pixelfreq_hclust_cap"].values:
                pos.loc[len(pos)] = [j, 0]
            if j not in neg["pixelfreq_hclust_cap"].values:
                neg.loc[len(neg)] = [j, 0]
            
        pos = pos.sort_values("pixelfreq_hclust_cap", axis = 0, ascending = True).reset_index()
        neg = neg.sort_values("pixelfreq_hclust_cap", axis = 0, ascending = True).reset_index()
        pos[0] = pos[0]/np.sum(pos[0].values)
        neg[0] = neg[0]/np.sum(neg[0].values)
        ratio = np.log2(pos[0]/neg[0])
        df = pd.DataFrame()
        df[i] = list(ratio.values)
        df = df.T
        
        datas = pd.concat([datas, df], ignore_index=True)
    
    # Remove Giant cells
    datas.columns = labels
    datas = datas.drop(["Giant_cell"], axis = 1)
    
    # Remove non-finite numbers
    data_list = datas.transpose().values.tolist()
    final = []
    order = []
    count = 0
    
    for lst in data_list:
        lst_new = []
        for i in lst:
            if i != float('-inf'):
                if i != float('inf'): # 0 lipoprotein- cell phenotype in an FOV
                    if i != float('nan'): # 0 of cell phenotype in an FOV
                        lst_new.append(i)
                        
        lst_2 = [j for j in lst_new if str(j) != 'nan'] #Confirm no 'NaN's remain
        final.append(lst_2)
        order.append(count)
        count = count + 1

    # Sort enrichment results by median
    med = []
    
    for lst in final:
        med.append(statistics.median(lst))
        
    res = {med[i]: labels[i] for i in range(len(med))}
    res2 = {med[i]: mycolors[i] for i in range(len(med))}
    res3 = {med[i]: order[i] for i in range(len(med))}
    labels_ordered = []
    mycolors_ordered = []
    list_ordered = []
    for k in sorted (res.keys()):
        labels_ordered.append(res[k])
        
    for k in sorted (res2.keys()):
        mycolors_ordered.append(res2[k])
        
    for k in sorted (res3.keys()):
        list_ordered.append(final[res3[k]])

    # Plot Enrichment Score across all FOVs per Cell Phenotype
    fig, ax = plt.subplots()
    bplot = plt.boxplot(list_ordered, patch_artist=True, labels = labels_ordered)
    plt.xticks(rotation = 90)
    
    # Fill with colors
    for patch, color in zip(bplot['boxes'], mycolors_ordered):
        patch.set_facecolor(color)
        
    plt.savefig(path + '/Sub-panels/Figure2/Fig2e_CP_LipoEnrichment.png', dpi = 1000)
        
if __name__ == '__main__':
     main()
