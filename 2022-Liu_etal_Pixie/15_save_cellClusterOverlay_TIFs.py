"""
Save cell phenotype map TIFs
Save TIF where pixel values correspond to cell cluster id and TIF colored according to cell cluster id

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

name = "pixelComposition" #name of output
clusters_path = "cellClustering_kcell14_passes10_cellRep1.csv" #output of cell clustering
clust_to_pheno_path = "cellClustering_kcell14_passes10_cellRep1_mapping_manual.csv" #csv mapping each cluster to its phenotype

clust_coln = "phenotype_num" #column name of cell cluster ids
seg_dir = "deepcell_output" #segmentation output
output_dir = "overlays_pxCellClusterIDs_"+name #output directory to store overlays
output_dir_colored = "overlays_pxCellClusterIDs_colored_"+name #output directory to store colored overlays
colors_path = "cell_colors.csv" #file specifying the colors for each cluster

# Make output directory
if not os.path.exists(output_dir):
    os.makedirs(output_dir)
if not os.path.exists(output_dir_colored):
    os.makedirs(output_dir_colored)

# Read data
clusters = pd.read_csv(clusters_path)
# Get manual mapping
clust_to_pheno = pd.read_csv(clust_to_pheno_path)
clusters = pd.merge(clusters, clust_to_pheno, on='pixelfreq_cluster')
# Get colors
colors_tab = pd.read_csv(colors_path)
clusters = pd.merge(clusters, colors_tab, on='phenotype')

# Get points
points = np.unique(clusters['sample'])
maxk = max(np.unique(clusters[clust_coln]))

## Create custom cmap
colors_only = colors_tab[colors_tab['phenotype_num'] <= maxk]
mycols = list(colors_only["color"])
mycols.insert(0,'#000000') #first color is black
mycols.append('#646464') #gray for cells that were not assigned
# Make bounds
bounds = [i-0.5 for i in np.linspace(0,maxk+2,maxk+3)]
colmap = colors.ListedColormap(mycols)
norm = colors.BoundaryNorm(bounds, colmap.N)

for i in range(len(points)):
    point = points[i]

    # Make dictionary of cell label:cluster id
    samp_clust = clusters.loc[clusters['sample'] == point]
    clust_dict = dict(zip(samp_clust['label'],samp_clust[clust_coln]))

    # Get cell segmentation overlay
    label_array = np.array(Image.open(os.path.join(seg_dir, "Point"+str(point)+"_feature_0.tif")))

    # Initialize array with number for unassigned
    pxfreq_cell_array = np.full(label_array.shape,maxk+1)
    # Assign 0s
    pxfreq_cell_array[np.where(label_array==0)] = 0
    # Assign labels
    for lab in clust_dict:
        idx = np.where(label_array == lab)
        pxfreq_cell_array[idx] = clust_dict[lab]

    # Change color of cell borders
    border_array = np.array(Image.open(os.path.join(seg_dir, "Point"+str(point)+"_overlay.tif")))
    white = (border_array==255).all(axis=2)
    border_idx = np.where(white==True)
    pxfreq_cell_array[border_idx] = 0

    # Save overlay as TIF
    im = Image.fromarray(pxfreq_cell_array.astype(np.int32))
    im.save(os.path.join(output_dir, "Point"+str(point)+"_pxCellClusterIDs.tiff"))

    # Save colored overlay
    image = colmap(norm(pxfreq_cell_array))
    plt.imsave(os.path.join(output_dir_colored, "Point"+str(point)+"_pxCellClusterIDs_colored.tiff"), image)

    print(point)

