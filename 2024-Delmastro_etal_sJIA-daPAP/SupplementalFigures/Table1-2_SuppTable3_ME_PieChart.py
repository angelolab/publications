#Install Libraries
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import colors
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

def main():
    # Set paths
    path = '/Users/aleadelmastro/Library/CloudStorage/GoogleDrive-alead@stanford.edu/.shortcut-targets-by-id/1GXGQlL1W-w8iJENEwcrx7NBGohHlnlCY/PAP/PAP_MIBI/Cohort/MIBI_data/all_points'
    cell_table = pd.read_csv(path + "/celltable_05032024.csv")
    fov_labels = pd.read_csv(path + "/fov_labels.csv")
    color_map = pd.read_csv(path + "/neighborhood_colorkey.csv")
    
    # Add disease category to cell_table
    cell_table = pd.merge(cell_table, fov_labels[["point","PatientID"]], on="point", how="left")

    # Define Color Map for Microenvironments
    topics = color_map.Topic
    mycolors = color_map.Hex

    # Keep only patient of Interest points
    patient = 13
    cell_table = cell_table[cell_table['PatientID'] == patient]
    
    #Calculate proportion of cells associated with a given ME
    values = []
    
    for topic in topics:
        cell_topic = len(cell_table[cell_table["topic"] == topic]) / len(cell_table)
        values.append(cell_topic)
        
    # Plot Pie Chart
    fig1, ax = plt.subplots()
    ax.pie(values, radius=1, colors = mycolors, counterclock=False, startangle = 90, wedgeprops = {"edgecolor" : "white",
                      'linewidth': 2,
                      'antialiased': True})
    
    my_circle=plt.Circle( (0,0), 0.5, color='white')
    p=plt.gcf()
    p.gca().add_artist(my_circle)
    ax.set(aspect="equal")
    plt.savefig(path + '/Sub-panels/SupplementalFigures/SuppTable3_patient' + str(patient) + '_ME.pdf', dpi=1000, bbox_inches='tight', transparent=True)
    
if __name__ == '__main__':
     main()