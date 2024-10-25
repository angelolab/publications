# Install Libraries
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import colors
import seaborn as sns
import statsmodels
from scipy import stats
import statsmodels.stats.multitest as smt

def main():
    # Set paths
    path = '/Users/aleadelmastro/Library/CloudStorage/GoogleDrive-alead@stanford.edu/.shortcut-targets-by-id/1GXGQlL1W-w8iJENEwcrx7NBGohHlnlCY/PAP/PAP_MIBI/Cohort/MIBI_data/all_points'
    cell_table = pd.read_csv(path + "/celltable_05032024.csv")
    fov_labels = pd.read_csv(path + "/fov_labels.csv")
    
    # Add disease category to cell_table
    cell_table = pd.merge(cell_table, fov_labels[["point","FOV_Category"]], on="point", how="left")
    
    # Choose ME topic
    topic = 1
    cell_table = cell_table[cell_table['topic'] == topic]
    
    markers = ['IFNg','HLA.DR', 'CD45RO','TIM3','pS6','CD57','MMP9','Ki67','GrzB','HO.1','IDO','iNOS','H3K9Ac','H3K27me3']
    cells = ['Bcell',"CD4+_Tcell","CD8+_Tcell","CD57+_CD8+_Tcell","Treg","CD11c+_mDC","CD14+_Mono",
                  "CD209+_Mac","iNOS+_Mac","M2_Mac","Eosinophil","Mast_cell","Neutrophil","CD16+_ImmuneOther",
                 "CD57+_ImmuneOther","Giant_cell"]

    # Get Lipoproteinosis+ cells
    lst = markers + ["name","point"]
    cell_Lipo = cell_table[cell_table['lipoprotein'] == 1]
    cell_Lipo = cell_Lipo[lst]
    
    # Get Lipoproteinosis- cells
    cell_noLipo = cell_table[cell_table['lipoprotein'] == 0]
    cell_noLipo = cell_noLipo[lst]
    p_val = pd.DataFrame(columns = ['Marker','p_val'])
    
    # Average marker expression across all lipoproteinosis category cells per FOV
    cell_Lipo = cell_Lipo.groupby('point').mean()
    cell_noLipo = cell_noLipo.groupby('point').mean()
    
    # Statistically compare marker expression between lipoproteinosis categories
    for marker in markers:
        _, p = stats.ttest_ind(cell_Lipo[marker].values, cell_noLipo[marker].values)
        df = pd.DataFrame([[marker,p]], columns = ['Marker','p_val'])
        p_val = pd.concat([p_val,df], ignore_index=True)
    
    _, p_adj, _, _ = smt.multipletests(p_val['p_val'].values, alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)
    p_val['p_adj'] = p_adj
    
    # Calculate log2 fold change of marker expression in each lipoproteinosis category
    cell_Lipo = cell_Lipo.mean().to_frame('Mean')
    cell_noLipo = cell_noLipo.mean().to_frame('Mean')
    
    cell_FC = cell_Lipo/cell_noLipo
    cell_FC = np.log2(cell_FC)
    
    p_val['log2FC'] = cell_FC['Mean'].values
    p_val['-log_p'] = -np.log10(p_val['p_adj'].values)
    
    # Generate Volcano Plot
    ## Define thresholds
    logFC_threshold = 0
    logP_threshold = 1.5
    for index, row in p_val.iterrows():
        if row['log2FC'] >= logFC_threshold and row['-log_p'] >= logP_threshold:
            p_val.at[index, 'cat'] = 'green'
        elif row['log2FC'] <= -logFC_threshold and row['-log_p'] >= logP_threshold:
            p_val.at[index, 'cat'] = 'green'
        else:
            p_val.at[index, 'cat'] = 'grey'

    enmax_palette = ['green','grey']
    g = sns.scatterplot(data = p_val, x = 'log2FC', y = '-log_p', hue = 'cat', palette = sns.set_palette(palette=enmax_palette), s = 100, edgecolor="black", linewidth=1)
    plt.legend([],[], frameon=False)
    plt.xlim(-2,  2)
    plt.ylim(-0.25,7.25)
    plt.rcParams["figure.figsize"] = (5,5)
    plt.savefig(path + '/Sub-panels/Figure4/Fig4j_DE_LipoCells_ME'+ str(topic) +'.png', dpi = 1000, transparent = True)
        
if __name__ == '__main__':
     main()