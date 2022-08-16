"""
Convert single-channel TIFs into pixel x markers matrix
Author: Candace Liu
Date: 8/15/22

"""

import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
from PIL import Image
import scipy.ndimage as ndimage

# Sigma for Gaussian blur
s = 2

# Markers for clustering
markers = ["CD14","CD209","HLA-DR-DQ-DP","CD4","MPO","CD3","SMA","CD11c","CD68","CD8","CD45","CD21","CD20","CD163","CD206","CD31"]
num_markers = len(markers)

# Number of pixels (nxn)
num_pixels = 1024

# Paths to single-channel TIFs and segmentation masks
base_dir = "mibi-final"
seg_dir = "deepcell_output"

# All points
points = list(range(1,13))
paths_all_points = [os.path.join(base_dir,"Point"+str(x),"TIFs") for x in points]

# Create directory to store output
output_dir = 'pixel_mats'
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# For each point, create pixel x markers matrix
for i in range(len(points)):
    # Point
    samp = points[i]
    # Path to data
    path = paths_all_points[i]

    # Initialize np array
    # 1 column for each marker, sample #, x, y, cell label from segmentation
    pixel_mat = np.zeros((num_pixels**2,num_markers+4))

    for j in range(len(markers)):
        mark = markers[j]
        im = Image.open(os.path.join(path,mark+".tif"))
        im_array = np.array(im)
        im_array = im_array.astype(np.float64)

        # Add Gaussian blur
        im_blur = ndimage.gaussian_filter(im_array, sigma=s)

        # Turn 1024x1024 -> 1048576x1
        img = im_blur.flatten()
        pixel_mat[:,j] = img

    # Add column for sample #
    pixel_mat[:,num_markers] = np.repeat(samp, num_pixels**2)

    # Add column for pixel x,y
    pixel_mat[:,num_markers+1] = np.repeat(range(num_pixels),num_pixels)
    pixel_mat[:,num_markers+2] = np.tile(range(num_pixels),num_pixels)

    # Add column for cell label from segmentation
    seg_labels = Image.open(os.path.join(seg_dir,"Point"+str(samp)+"_feature_0.tif"))
    seg_labels_array = np.array(seg_labels)
    pixel_mat[:,num_markers+3] = seg_labels_array.flatten()

    # Remove all 0s
    pixel_mat_nozeros = pixel_mat[np.sum(pixel_mat[:,:num_markers], axis=1)!=0,:]

    # Write to file
    header = ','.join(markers)
    header = header+',sample,x,y,label'
    np.savetxt(os.path.join(output_dir,"Point"+str(samp)+"_sigma"+str(s)+".csv"), pixel_mat_nozeros, delimiter=",", header=header, comments='')

    print(samp)

