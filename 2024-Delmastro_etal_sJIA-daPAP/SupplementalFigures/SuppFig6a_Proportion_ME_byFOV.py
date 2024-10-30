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
    color_map = pd.read_csv(path + "/neighborhood_colorkey.csv")
    
    # Add patient ids to cell_table
    cell_table = pd.merge(cell_table, fov_labels[["point","PatientID.Block-FOV"]], on="point", how="left")

    # Define Color Map for Microenvironments
    topics = color_map.Topic
    mycolors = color_map.Hex   
    cmap = dict(zip(topics, mycolors))
   
    # Calculate proportion of cells associated with each microenvironment per FOV
    patient_df = pd.DataFrame(columns = topics)
    IDs = cell_table["PatientID.Block-FOV"].unique()
    for ID in IDs:
        cell_patient = cell_table[cell_table['PatientID.Block-FOV'] == ID]     
        total_topic = len(cell_patient)
        cell_topic = cell_patient.groupby(['topic']).size().reset_index()
        for topic in topics:
            if topic not in cell_topic["topic"].values:
                cell_topic.loc[len(cell_topic)] = [topic, 0]
        
        cell_topic.sort_values("topic", axis = 0, ascending = True,inplace = True)
        
        proportion = cell_topic[0].values/ total_topic
        patient_df.loc[len(patient_df)] = proportion
     
    # Add patient ids to dataframe
    patient_df['patient'] = IDs
        
    # Sort table by patient ID
    patient_df = patient_df.sort_values(by=['patient'], key = lambda x: np.argsort(index_natsorted(patient_df["patient"])))
    ME_order = list(patient_df.columns)
    ME_order.pop()
    
    #Plot stacked bar plot of microenvironment proportions for each FOV
    plt.rcParams["figure.figsize"] = (15,5)
    patient_df.plot(x='patient', kind='bar', stacked=True,
    title='Proportion of Microenvironments', legend = None, color = [cmap[key] for key in ME_order], width = 0.9, ylim=(0,1)).get_figure().savefig(path + "/Sub-panels/SupplementalFigures/SuppFig6a_Proportion_ME_all.png", dpi = 1000)
    plt.ylim(0,1)
    
if __name__ == '__main__':
    main()