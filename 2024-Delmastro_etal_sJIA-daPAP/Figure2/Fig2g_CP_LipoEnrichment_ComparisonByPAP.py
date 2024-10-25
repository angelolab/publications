# Install Libraries
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import colors
import seaborn as sns
import statistics
from scipy import stats

def main():
    # Set paths
    path = '/Users/aleadelmastro/Library/CloudStorage/GoogleDrive-alead@stanford.edu/.shortcut-targets-by-id/1GXGQlL1W-w8iJENEwcrx7NBGohHlnlCY/PAP/PAP_MIBI/Cohort/MIBI_data/all_points'
    cell_table = pd.read_csv(path + "/celltable_05032024.csv")
    fov_labels = pd.read_csv(path + "/fov_labels.csv")
    color_map = pd.read_csv(path + "/cellpheno_num_colorkey.csv")
    
    # Cells of Interest for Comparison
    cellInterest = ['M2_Mac','Neutrophil','Eosinophil']
    
    # Add disease category to cell_table
    cell_table = pd.merge(cell_table, fov_labels[["point","FOV_Category"]], on="point", how="left")
    
    # Define Color Map for Cell Phenotypes
    labels = color_map.Pheno
    mycolors = color_map.Hex

    # Keep only PAP involved points
    cats = ['sJIA-daPAP','non-sJIA-PAP']
    cell_table_PAP = cell_table[cell_table['FOV_Category'].isin(cats)]
    
    # Differentiate lipoprotein+ and Lipoprotein- cells and groupby cell phenotype
    pos_cells = cell_table_PAP[cell_table_PAP["lipoprotein"] == 1]
    neg_cells = cell_table_PAP[cell_table_PAP["lipoprotein"] == 0]
    pos_cells = pos_cells.groupby(["point","name",'FOV_Category']).size().reset_index()
    neg_cells = neg_cells.groupby(["point","name",'FOV_Category']).size().reset_index()
    
    data = pd.DataFrame(columns = ['point','name','FOV_Category','Enrichment'])
    
    # Calculate enrichment of Lipoprotein+/Lipoprotein- per cell phenotype
    for cat in cats:
        points = fov_labels[fov_labels["FOV_Category"] == cat].point.values
            
        for i in points:
            pos = pos_cells[pos_cells['point'] == i]
            neg = neg_cells[neg_cells['point'] == i]
            pos = pos[["name",0]]
            neg = neg[["name",0]]

            for cell in labels:
                if cell not in pos["name"].values:
                    df = pd.DataFrame([[cell, 0]], columns = ["name",0])
                    pos = pd.concat([pos,df], ignore_index=True)
                if cell not in neg["name"].values:
                    df = pd.DataFrame([[cell, 0]], columns = ["name",0])
                    neg = pd.concat([neg,df], ignore_index=True)

            #Remove Giant cells
            pos = pos[pos["name"] != "Giant_cell"]
            neg = neg[neg["name"] != "Giant_cell"]
                
            pos = pos.sort_values("name", axis = 0, ascending = True).reset_index()
            neg = neg.sort_values("name", axis = 0, ascending = True).reset_index()

            pos[0] = pos[0]/np.sum(pos[0].values)
            neg[0] = neg[0]/np.sum(neg[0].values)
            ratio = np.log2(pos[0]/neg[0])
            
            for k in range(len(ratio)):
                name = pos.iloc[k]['name']
                enr = (ratio.values)[k]
                df = pd.DataFrame([[i,name,cat,enr]], columns = ["point","name","FOV_Category","Enrichment"])
                data = pd.concat([data,df], ignore_index=True)

    data = data.dropna(axis=0)
    data = data[~data["Enrichment"].isin([float('inf'),float('-inf')])]
    
    # Compare and Plot Only the Cell Types of Interest
    for cell in cellInterest:
        data_test = data[data['name'] == cell]
        data_sJIA = data_test[data_test['FOV_Category'] == 'sJIA-daPAP']['Enrichment'].values
        data_nonSJIA = data_test[data_test['FOV_Category'] == 'non-sJIA-PAP']['Enrichment'].values
        _, p = stats.ttest_ind(data_sJIA, data_nonSJIA)
        print("sJIA-daPAP vs. non-sJIA-PAP: ",cell,": ", p)
        
    data = data[data["name"].isin(cellInterest)]
        
    fig1, ax = plt.subplots()
    ax = sns.boxplot(x="name", y="Enrichment", hue="FOV_Category",
                   data=data, palette="tab10", dodge=True, showfliers = False)
    plt.legend([],[], frameon=False)
    plt.rcParams["figure.figsize"] = (5,5)
    ax.set_ylim(-3,8)
    plt.savefig(path + '/Sub-panels/Figure2/Fig2g_CP_LipoEnrichment_ComparisonByPAP.png', dpi = 1000)
    
if __name__ == '__main__':
     main()