"""
Save pixel phenotype map TIFs with cell border overlaid
Save TIF were pixel values correspond to pixel cluster id and TIF colored according to pixel cluster id

Author: Candace Liu
Date: 8/15/22

"""

import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib import colors
from PIL import Image
import scipy.ndimage as ndimage

name = "sigma2_passes10_rep1" #name of file output from pixelClustering.R
clusters_path = "pixelClustering_"+name+"_clusters.csv"
clust_coln = "hCluster_cap" #column name of pixel cluster ids
output_dir = "overlays_pxClusterIDs_"+name #output directory to store overlays
output_dir_colored = "overlays_pxClusterIDs_colored_"+name #output directory to store colored overlays
clust_to_pheno_path = "pixelClustering_"+name+"_mapping.csv" #csv mapping each cluster to its phenotype
colors_path = "px_colors.csv" #file specifying the colors for each cluster
size = 1024 #size of image

# Make output directory
if not os.path.exists(output_dir):
    os.makedirs(output_dir)
if not os.path.exists(output_dir_colored):
    os.makedirs(output_dir_colored)

# Read data
clusters = pd.read_csv(clusters_path, sep=',')
# Point names
points = list(range(1,13))
# Total number of clusters
maxk = max(np.unique(clusters[clust_coln]))

# Get phenotype mapping
clust_to_pheno = pd.read_csv(clust_to_pheno_path)
colors_tab = pd.read_csv(colors_path)
clust_to_color = pd.merge(clust_to_pheno, colors_tab, on='phenotype')
clust_to_color = clust_to_color.sort_values(['hCluster_cap'])

## Create custom cmap
mycols = list(clust_to_color["color"])
mycols.insert(0,'#000000') #first color is black
# Make bounds
bounds = [i-0.5 for i in np.linspace(0,maxk+1,maxk+2)]
colmap = colors.ListedColormap(mycols)
norm = colors.BoundaryNorm(bounds, colmap.N)

# Define function to assign labels
def assign_label(clust_array,x,y,clust):
    clust_array[x,y] = clust

for i in range(len(points)):
    point = points[i]

    # Get cluster data for one point
    samp_clust = clusters.loc[clusters['sample'] == point]

    # Create array for cluster labels, fill with 0's
    clust_array = np.full((size,size), 0, dtype=int)
    # Fill in array
    [assign_label(clust_array,row[0],row[1],row[2]) for row in samp_clust[['x','y',clust_coln]].values]

    # Save overlay as TIF
    im = Image.fromarray(clust_array.astype(np.int32))
    im.save(os.path.join(output_dir, "Point"+str(point)+"_pxClusterIDs.tiff"))

    # Save colored overlay
    image = colmap(norm(clust_array))
    plt.imsave(os.path.join(output_dir_colored, "Point"+str(point)+"_pxClusterIDs_colored.tiff"), image)

    print(point)

