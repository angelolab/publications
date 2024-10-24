import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

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
    cmap = dict(zip(labels, mycolors))

    ##USER INPUT HERE:
    cat = 'uninvolved'

    cell_table = cell_table[cell_table['FOV_Category'] == cat]

    # Define immune cells and cell order for pie chart
    immune_cells = ['Bcell',"CD4+_Tcell","CD8+_Tcell","CD57+_CD8+_Tcell","Treg","CD11c+_mDC","CD14+_Mono",
                  "CD209+_Mac","iNOS+_Mac","M2_Mac","Giant_cell", "Eosinophil","Mast_cell","Neutrophil","CD16+_ImmuneOther",
                 "CD57+_ImmuneOther"]
    cells_ordered = ['Bcell',"CD4+_Tcell","CD8+_Tcell","CD57+_CD8+_Tcell","Treg","CD11c+_mDC","CD14+_Mono",
                  "CD209+_Mac","iNOS+_Mac","M2_Mac","Giant_cell", "Eosinophil","Mast_cell","Neutrophil","CD16+_ImmuneOther",
                 "CD57+_ImmuneOther","Endothelial", "Fibroblast","iNOS+_Pneumocyte","Lung_Epithelium","Mesenchymal"]

    # Calculate proportion of cells associated with a given phenotype and whether immune or non-immune
    nonimmune_values = []
    outer_nonimmune_colors = []
    immune_values = []
    outer_immune_colors = []
    for cell in cells_ordered:
        cell_prop = len(cell_table[cell_table["name"] == cell]) / len(cell_table)
        if cell in immune_cells:
            immune_values.append(cell_prop)
            outer_immune_colors.append(cmap[cell])
        else:
            nonimmune_values.append(cell_prop)
            outer_nonimmune_colors.append(cmap[cell])

    inner_colors = ["#808080", "#D3D3D3"]

    size = 0.35
    vals_all = np.array(immune_values + nonimmune_values)
    outer_colors = outer_immune_colors + outer_nonimmune_colors
    vals_sum = np.array([np.sum(np.array(immune_values)) / np.sum(vals_all), np.sum(np.array(nonimmune_values)) / np.sum(vals_all)])

    # Plot the pie chart
    fig, ax = plt.subplots()
    
    ax.pie(vals_all,
           radius=1,
           colors=outer_colors,
           startangle=90,
           counterclock=False,
           wedgeprops={"edgecolor": "white",
                       'linewidth': 1,
                       'width': size})

    ax.pie(vals_sum,
           radius=1 - size,
           colors=inner_colors,
           startangle=90,
           counterclock=False,
           wedgeprops={"edgecolor": "white",
                       'linewidth': 1,
                       'width': size})

    ax.set(aspect="equal")
    plt.savefig(path + '/Sub-panels/Figure2/Fig2c_pie_' + cat + '.pdf', dpi=1000, bbox_inches='tight', transparent=True)

if __name__ == '__main__':
    main()
