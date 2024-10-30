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
    cell_table = pd.merge(cell_table, fov_labels[["point","PatientID"]], on="point", how="left")
    
    # Define Color Map for Cell Phenotypes
    labels = color_map.Pheno
    mycolors = color_map.Hex
    cmap = dict(zip(labels, mycolors))

    ##USER INPUT HERE:
    patient = 13
    cell_table = cell_table[cell_table['PatientID'] == patient]

    # Define immune cells and cell order for pie chart
    immune_cells = ['Bcell',"CD4+_Tcell","CD8+_Tcell","CD57+_CD8+_Tcell","Treg","CD11c+_mDC","CD14+_Mono",
                  "CD209+_Mac","iNOS+_Mac","M2_Mac","Giant_cell", "Eosinophil","Mast_cell","Neutrophil","CD16+_ImmuneOther",
                 "CD57+_ImmuneOther"]
    cells_ordered = ['Bcell',"CD4+_Tcell","CD8+_Tcell","CD57+_CD8+_Tcell","Treg","CD11c+_mDC","CD14+_Mono",
                  "CD209+_Mac","iNOS+_Mac","M2_Mac","Giant_cell", "Eosinophil","Mast_cell","Neutrophil","CD16+_ImmuneOther",
                 "CD57+_ImmuneOther","Endothelial", "Fibroblast","iNOS+_Pneumocyte","Lung_Epithelium","Mesenchymal"]

    # Calculate proportion of cells associated with a given phenotype and whether immune or non-immune
    colors = []
    values = []
    for cell in cells_ordered:
           values.append(len(cell_table[cell_table["name"] == cell]) / len(cell_table))
           colors.append(cmap[cell])

    # Plot the pie chart
    fig, ax = plt.subplots()
    
    ax.pie(values,
           radius=1,
           colors=colors,
           startangle=90,
           counterclock=False,
           wedgeprops = {"edgecolor" : "white",
                      'linewidth': 2,
                      'antialiased': True})
    my_circle=plt.Circle( (0,0), 0.5, color='white')
    p=plt.gcf()
    p.gca().add_artist(my_circle)
    ax.set(aspect="equal")
    plt.savefig(path + '/Sub-panels/SupplementalFigures/SuppTable3_patient' + str(patient) + '_CPM.pdf', dpi=1000, bbox_inches='tight', transparent=True)

if __name__ == '__main__':
    main()
