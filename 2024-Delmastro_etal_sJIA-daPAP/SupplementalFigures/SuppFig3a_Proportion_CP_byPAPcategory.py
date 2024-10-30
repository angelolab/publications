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
    
    # Order cells in alphabetical order
    cellOrder = ['Bcell','CD11c+_mDC','CD14+_Mono','CD16+_ImmuneOther','CD209+_Mac','CD4+_Tcell','CD57+_CD8+_Tcell',
                'CD57+_ImmuneOther','CD8+_Tcell','Endothelial','Eosinophil','Fibroblast','Giant_cell','iNOS+_Mac',
                'iNOS+Pneumocyte','Lung_Epithelium','M2_Mac','Mast_cell','Mesenchymal','Neutrophil','Treg']
    
    # Add disease category to cell_table
    cell_table = pd.merge(cell_table, fov_labels[["point","FOV_Category"]], on="point", how="left")
    
    # Define Color Map for Cell Phenotypes
    labels = color_map.Pheno
    mycolors = color_map.Hex

    # Keep only PAP involved points
    cats = ['sJIA-daPAP','non-sJIA-PAP']
    cell_table_PAP = cell_table[cell_table['FOV_Category'].isin(cats)]
    
    # Differentiate lipoprotein+ and Lipoprotein- cells and groupby cell phenotype
    data = cell_table_PAP.groupby('point')['name'].value_counts(normalize=True).unstack()
    data = data.fillna(0)
    data = data.reset_index().melt(id_vars='point', var_name='name', value_name='proportion')
    data = pd.merge(data, fov_labels[["point","FOV_Category"]], on="point", how="left")
    
    # Compare and Plot Only the Cell Types of Interest
    for cell in cellOrder:
        data_test = data[data['name'] == cell]
        data_sJIA = data_test[data_test['FOV_Category'] == 'sJIA-daPAP']['proportion'].values
        data_nonSJIA = data_test[data_test['FOV_Category'] == 'non-sJIA-PAP']['proportion'].values
        _, p = stats.ttest_ind(data_sJIA, data_nonSJIA)
        print("sJIA-daPAP vs. non-sJIA-PAP: ",cell,": ", p)
        
    fig1, ax = plt.subplots()
    ax = sns.boxplot(x="name", y="proportion", hue="FOV_Category", order = cellOrder,
                   data=data, palette="tab10", dodge=True, showfliers = False)
    plt.legend([],[], frameon=False)
    plt.rcParams["figure.figsize"] = (15,5)
    plt.xticks(rotation=45)
    ax.set_ylim(0,0.5)
    plt.savefig(path + '/Sub-panels/SupplementalFigures/SuppFig3a_CP_ComparisonByPAP.pdf', dpi = 1000)
    
if __name__ == '__main__':
     main()