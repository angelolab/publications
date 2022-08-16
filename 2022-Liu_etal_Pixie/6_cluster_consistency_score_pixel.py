"""
After cluster_consistency_score_pixel.R, finish calculation here
Find score for each pixel, create overlays colored according to cluster consistency score
Author: Candace Liu
Date: 8/15/22

"""

import pandas as pd
import numpy as np
import os
from PIL import Image
import matplotlib.pyplot as plt


# Variables to change
name = "sigma2_passes10_threshold95" #name of file created in cluster_consistency_score_pixel.R
reps = list(range(1,6)) #all replicates
comparison_coln = ["rep"+str(x)+"_purity" for x in reps]
size = 1024 #size of image
output_dir_colored = "purity_tiff_colored_"+name #directory to store overlays colored according to score
minval = 1 #minimum value for colored overlays
maxval = 5 #maximum value for colored overlays

# Make output directory
if not os.path.exists(output_dir_colored):
    os.makedirs(output_dir_colored)

# All points
points = list(range(1,13))

# Get cluster consistency scores (from cluster_consistency_score_pixel.R)
dat = pd.read_csv('purity_'+name+'.csv')
# Get mean and sd for each pixel
dat['mean_purity'] = dat[comparison_coln].mean(axis=1)
dat['sd_purity'] = dat[comparison_coln].std(axis=1)
dat.to_csv('meanpurity_'+name+'.csv', columns=['sample','x','y','mean_purity','sd_purity'], index=False)

# Plot histogram of all pixels
fig = plt.figure(figsize=(10,10))
plt.hist(dat['mean_purity'])
plt.savefig(name+'_allPixels_purity.tiff')
plt.close()

# Define function to assign values
def assign_label(clust_array,x,y,clust):
    clust_array[int(x),int(y)] = clust

# cmap for colored TIFF
cmap = plt.get_cmap('inferno')
cmap.set_bad(color='gray')

for point in points:
    # Get data for one point
    samp_dat = dat.loc[dat['sample'] == point]
    
    # Create array for cluster labels, fill with 0's
    purity_array = np.zeros((size,size))
    # Fill in array
    [assign_label(purity_array,row[0],row[1],row[2]) for row in samp_dat[['x','y','mean_purity']].values]
   
    # Change color of 0 to gray
    purity_array_mask = np.ma.masked_where(purity_array == 0, purity_array)
    # Save colored TIFF
    fig = plt.figure(figsize=(10,10))
    plt.imshow(purity_array_mask, cmap=cmap, vmin=minval, vmax=maxval)
    plt.colorbar(fraction=0.04)
    plt.axis('off')
    fig.suptitle("Point"+str(point))
    plt.savefig(os.path.join(output_dir_colored, "Point"+str(point)+"_purity.tiff"), bbox_inches='tight')
    plt.close()
    
    print(point)
    

