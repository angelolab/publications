# Install Libraries
import seaborn as sns
from scipy import stats
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import colors

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

    #List of all cell types in dataset for reference
    #all_cells = ["Endothelial", "CD11c+_mDC", "CD57+_ImmuneOther", "iNOS+_Mac", "iNOS+_Pneumocyte", "Fibroblast", "CD209+_Mac", "CD8+_Tcell", "Neutrophil", "CD14+_Mono", "Bcell", "CD4+_Tcell", "M2_Mac", "Treg", "Mesenchymal", "Lung_Epithelium", "Mast_cell", "CD57+_CD8+_Tcell", "Eosinophil", "CD16+_ImmuneOther","Giant_cell"]
    
    # Cells of Interest for Comparison
    all_cells = ['Neutrophil']
    
    # Functional markers of interest
    markers = ['MMP9']
            
    #Perform comparison for each cell type and marker of interest
    for cell in all_cells:
        cell_table_name = cell_table[cell_table["name"] == cell].drop(columns = ["name"])
        for marker in markers:
            df = cell_table_name[[marker,"FOV_Category","lipoprotein"]]

            # Compare between disease category and lipoproteinosis status
            lst_lipo_nonSJIA = np.array(df[(df["lipoprotein"] == 1) & (df["FOV_Category"] == "non-sJIA-PAP")].drop(columns=["lipoprotein","FOV_Category"])[marker].values)
            lst_lipo_SJIA = np.array(df[(df["lipoprotein"] == 1) & (df["FOV_Category"] == "sJIA-daPAP")].drop(columns=["lipoprotein","FOV_Category"])[marker].values)
            lst_nolipo_nonSJIA = np.array(df[(df["lipoprotein"] == 0) & (df["FOV_Category"] == "non-sJIA-PAP")].drop(columns=["lipoprotein","FOV_Category"])[marker].values)
            lst_nolipo_SJIA = np.array(df[(df["lipoprotein"] == 0) & (df["FOV_Category"] == "sJIA-daPAP")].drop(columns=["lipoprotein","FOV_Category"])[marker].values)

            #Plot violin plot
            cmap = ["#000000","#FFFFFF"]
            fig = plt.figure()
            fig.set_size_inches(4, 6)
            ax = sns.violinplot(x="FOV_Category", y=marker, hue="lipoprotein", data=df, palette=cmap, scale = "width", inner = "quartile", legend = False)
            plt.legend([], [], frameon=False)
            ax.set_yticks(np.arange(-0.2, 1.8, step=0.2))
            plt.ylim(-0.3,1.7)
            plt.show()

            #Perform t-test between distribution of intensities
            print(marker)
            _, p = stats.ttest_ind(lst_lipo_SJIA, lst_lipo_nonSJIA)
            print("Lipoprotein: sJIA-daPAP vs. non-sJIA-PAP: ", p)
            _, p = stats.ttest_ind(lst_nolipo_SJIA, lst_nolipo_nonSJIA)
            print("No Lipoprotein: sJIA-daPAP vs. non-sJIA-PAP: ", p)
            _, p = stats.ttest_ind(lst_lipo_SJIA, lst_nolipo_SJIA)
            print("sJIA-daPAP: Lipoprotein vs. no Lipoprotein: ", p)
            _, p = stats.ttest_ind(lst_lipo_nonSJIA, lst_nolipo_nonSJIA)
            print("non-sJIA-PAP: Lipoprotein vs. no Lipoprotein: ", p)

            fig.savefig(path + "/Sub-panels/SupplementalFigures/SuppFig5f_" + str(cell) + "_"+ str(marker) + ".png", dpi = 1000, transparent = True)

if __name__ == '__main__':
     main()