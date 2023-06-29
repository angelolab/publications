"""
Save preprocessed data as image
Author: Candace Liu
Date: 1/4/23

"""

import numpy as np
import pandas as pd
import os
from PIL import Image

output_dir = "preprocessed_single_channels"
markers = ["CD14","CD209","HLA-DR-DQ-DP","CD4","MPO","CD3","SMA","CD11c","CD68","CD8","CD45","CD21","CD20","CD163","CD206","CD31"]
size = 1024

# Get normalization values
norm_vals = pd.read_csv("avg999_sigma2.csv")

# Function to assign pixel values
def assign_pixels(im_array,x,y,val):
    im_array[x,y] = val

# All points
points = list(range(1,13))

# Make output directories
for p in points:
    out_path = os.path.join(output_dir,"Point"+str(p))
    if not os.path.exists(out_path):
        os.makedirs(out_path)

# Process data
for point in points:
    # Get blurred data
    dat = pd.read_csv(os.path.join("pixel_mats","Point"+str(point)+"_sigma2.csv"))

    # Split up data
    dat_meta = dat[["sample","x","y","label"]]
    dat_markers = dat[markers]

    # Frequency normalization
    dat_markers = dat_markers.div(dat_markers.sum(axis=1), axis=0)

    # 99.% normalization
    dat_markers = dat_markers.div(norm_vals[markers].iloc[0], axis=1)
    save_dat = pd.concat([dat_meta,dat_markers], axis=1)

    # Save as image
    for marker in markers:
        new_im = np.zeros((size,size))
        # Assign new values
        [assign_pixels(new_im,int(row[0]),int(row[1]),row[2]) for row in save_dat[['x','y',marker]].values]
        im = Image.fromarray(new_im)
        im.save(os.path.join(output_dir, "Point"+str(point), marker+".tiff"))

    print(point) 
    

