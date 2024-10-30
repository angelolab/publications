#Install Libraries
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import colors
import seaborn as sns
from scipy import stats
from natsort import index_natsorted

def main():
    # Set paths
    path = '/Users/aleadelmastro/Library/CloudStorage/GoogleDrive-alead@stanford.edu/.shortcut-targets-by-id/1GXGQlL1W-w8iJENEwcrx7NBGohHlnlCY/PAP/PAP_MIBI/Cohort/MIBI_data/all_points'
    cell_table = pd.read_csv(path + "/celltable_05032024.csv")
    fov_labels = pd.read_csv(path + "/fov_labels.csv")
    color_map = pd.read_csv(path + "/cellpheno_num_colorkey.csv")
    
    # Add patient ids to cell_table
    cell_table = pd.merge(cell_table, fov_labels[["point","PatientID.Block-FOV"]], on="point", how="left")

    # Define Color Map for Cell Phenotypes
    phenos = color_map.Pheno
    mycolors = color_map.Hex   
    cmap = dict(zip(phenos, mycolors))
   
    # Calculate proportion of cells associated with each cell phenotype per FOV
    patient_df = pd.DataFrame(columns = sorted(phenos, key=str.lower))
    IDs = cell_table["PatientID.Block-FOV"].unique()
    for ID in IDs:
        cell_patient = cell_table[cell_table['PatientID.Block-FOV'] == ID]     
        total_cell = len(cell_patient)
        cell_pheno = cell_patient.groupby(['name']).size().reset_index()
        for cell in phenos:
            if cell not in cell_pheno["name"].values:
                cell_pheno.loc[len(cell_pheno)] = [cell, 0]
        
        cell_pheno.sort_values("name", axis = 0, ascending = True,
                 inplace = True, key=lambda col: col.str.lower())
        
        proportion = cell_pheno[0].values/ total_cell
        patient_df.loc[len(patient_df)] = proportion
     
    # Add patient ids to dataframe
    patient_df['patient'] = IDs
        
    # Sort table by patient ID
    patient_df = patient_df.sort_values(by=['patient'], key = lambda x: np.argsort(index_natsorted(patient_df["patient"])))
    cell_order = list(patient_df.columns)
    cell_order.pop()
    
    #Plot stacked bar plot of cell phenotype proportions for each FOV
    plt.rcParams["figure.figsize"] = (15,5)
    patient_df.plot(x='patient', kind='bar', stacked=True,
    title='Proportion of Cells', legend = None, color = [cmap[key] for key in cell_order], width = 0.9, ylim=(0,1)).get_figure().savefig(path + "/Sub-panels/SupplementalFigures/SuppFig3b_Proportion_CellPhenotypes_all.png", dpi = 1000)
    plt.ylim(0,1)
    
if __name__ == '__main__':
    main()