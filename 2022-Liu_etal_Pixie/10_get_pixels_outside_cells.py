# Get percentage of pixels outside of cells
# Author: Candace Liu
# Date: 8/15/22

import pandas as pd
import numpy as np
from PIL import Image
import os

px_dir = "overlays_pxClusterIDs_sigma2_passes10_rep1" #name of directory that contains overlay TIFs where numbers are cluster numbers
seg_dir = "deepcell_output" #segmentation output directory
size = 1024 #size of images (nxn)

# All points
points = list(range(1,13))

all_counts = {}
for i in range(len(points)):
    point = points[i]

    # Get px cluster overlay
    px = np.array(Image.open(os.path.join(px_dir, "Point"+str(point)+"_pxClusterIDs.tiff")))
    # Get seg overlay
    seg = np.array(Image.open(os.path.join(seg_dir, "Point"+str(point)+"_feature_0.tif")))

    # Dummy matrix
    counts = np.ones(px.shape)
    # Set where pixels are 0 to 0
    counts[px==0] = 0
    # Set inside cells to 0
    counts[seg!=0] = 0

    # Get number of pixels outside of cells
    all_counts[point] = int(np.sum(counts))

# Turn into data frame
df = pd.DataFrame(all_counts.items(), columns=['point','counts'])
df['total'] = size*size
df['perc'] = df['counts']/df['total']
df['group'] = "ln"

df.to_csv('pixels_outside_cells.csv', index=False)

