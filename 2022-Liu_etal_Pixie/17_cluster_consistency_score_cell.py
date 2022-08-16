"""
After cluster_consistency_score_cell.R, finish calculation here
Find score for each cell, create overlays colored according to cluster consistency score
Author: Candace Liu
Date: 8/15/22

"""

import pandas as pd
import numpy as np
import os
from PIL import Image
import matplotlib.pyplot as plt


# Variables to change
name = "pixelComposition" #filename from cluster_consistency_score_cell.R
reps = list(range(1,6)) #all replicates

comparison_coln = ["rep"+str(x)+"_purity" for x in reps]
output_dir_colored = "cell_"+name+"_purity_tiff_colored"
seg_dir = "deepcell_output" #where segmentation output is stored

# Make output directory
if not os.path.exists(output_dir_colored):
    os.makedirs(output_dir_colored)

# All points
points = list(range(1,13))

# Get cluster consistency scores (from cluster_consistency_score_cell.R)
dat = pd.read_csv('cell_'+name+'_purity.csv')
dat['mean_purity'] = dat[comparison_coln].mean(axis=1)
dat.to_csv('cell_'+name+'_meanpurity.csv', columns=['sample','label','mean_purity'], index=False)

# Plot histogram of all cells
fig = plt.figure(figsize=(10,10))
plt.hist(dat['mean_purity'])
plt.savefig('cell_'+name+'_allCells_purity.tiff')
plt.close()

# cmap for colored TIFF
cmap = plt.get_cmap('inferno')
cmap.set_bad(color='gray')

for point in points:
    # Get data for one point
    samp_dat = dat.loc[dat['sample'] == point]
    # Make dictionary of cell label:purity score
    cell_dict = dict(zip(samp_dat['label'],samp_dat['mean_purity']))
   
    # Get cell segmentation overlay
    label_array = np.array(Image.open(os.path.join(seg_dir,"Point"+str(point)+"_feature_0.tif")))

    # Initialize with zeros
    purity_array = np.zeros(label_array.shape)
    # Assign labels
    for lab in cell_dict:
        idx = np.where(label_array == lab)
        purity_array[idx] = cell_dict[lab]

    # Change cell borders
    border_array = np.array(Image.open(os.path.join(seg_dir, "Point"+str(point)+"_overlay.tif")))
    white = (border_array==255).all(axis=2)
    border_idx = np.where(white==True)
    purity_array[border_idx] = 0

    # Change color of 0 to gray
    purity_array_mask = np.ma.masked_where(purity_array == 0, purity_array)

    # Save colored TIFF
    fig = plt.figure(figsize=(10,10))
    plt.imshow(purity_array_mask, cmap=cmap, vmin=1, vmax=5)
    plt.colorbar(fraction=0.04)
    plt.axis('off')
    fig.suptitle("Point"+str(point))
    plt.savefig(os.path.join(output_dir_colored, "Point"+str(point)+"_purity.tiff"), bbox_inches='tight')
    plt.close()
    
    print(point)
    

