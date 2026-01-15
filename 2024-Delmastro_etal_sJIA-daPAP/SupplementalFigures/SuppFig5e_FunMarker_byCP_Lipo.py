# Install Libraries
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import colors
import seaborn as sns
import statistics 
import scipy.stats as stats

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
    
    # Cell of interest
    cell = 'Lung_Epithelium'
    cell_table = cell_table[cell_table["name"] == cell].drop(columns = ["name"])
    
    # Functional markers of interest
    markers = ['IFNg']
       
    # Keep only Lipoproteinosis+ Cells
    cell_table = cell_table[cell_table['lipoprotein'] == 1]
    
    #Generate Violin Plots for each functional marker of interest
    for marker in markers:
        df = cell_table[[marker,"FOV_Category"]]
        lst_nonSJIA = np.array(df[df["FOV_Category"] == "non-sJIA-PAP"].drop(columns=["FOV_Category"])[marker].values)
        lst_SJIA = np.array(df[df["FOV_Category"] == "sJIA-daPAP"].drop(columns=["FOV_Category"])[marker].values)

        #Plot violin plot
        fig = plt.figure()
        fig.set_size_inches(4, 6)
        ax = sns.violinplot(data = [lst_SJIA, lst_nonSJIA], palette = 'tab10', scale = "width", inner = "quartile")
        #ax = sns.boxplot(data = [lst_SJIA, lst_PAP], palette = 'tab10')
        #ax = sns.swarmplot(data = [lst_SJIA, lst_PAP], palette = 'tab10')
        plt.title(marker  + " Expression in Lipoproteinosis+ " + cell + "s across Tissue Categories ")
        plt.ylabel("Normalized " + marker + " Expression")
        plt.ylim(-0.3,1.5)
        plt.show()
    
        #Perform t-test between distribution of intensities
        print(marker)
        _, p = stats.ttest_ind(lst_SJIA, lst_nonSJIA)
        print("sJIA-daPAP vs. non-sJIA-PAP: ", p)
        
        fig.savefig(path + "/Sub-panels/SupplementalFigures/SuppFig5e_Lipo+_" + str(cell) + "_"+ str(marker) + ".png", dpi = 1000)

if __name__ == '__main__':
     main()