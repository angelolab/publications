"""
Apply Otsu thresholding to images
Author: Candace Liu
Date: 11/30/22

"""

import os
import numpy as np
import pandas as pd
import skimage.io as io
from skimage.filters import threshold_otsu
import scipy.ndimage as ndimage
from PIL import Image
import feather

# sigma for Gaussian blur
s = 2

# Markers for clustering
markers = ["CD14","CD209","HLA-DR-DQ-DP","CD4","MPO","CD3","SMA","CD11c","CD68","CD8","CD45","CD21","CD20","CD163","CD206","CD31"]
num_markers = len(markers)

# Number of pixels (nxn)
num_pixels = 1024

# Paths to single-channel TIFs and segmentation masks
base_dir = "mibi-final"

# All points
points = list(range(1,13))
paths_all_points = [os.path.join(base_dir,"Point"+str(x),"TIFs") for x in points]

# Create directory to store pixel matrices
output_dir = 'pixel_mats_thresholded'
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# Create directories to store thresholded images
image_output_dir = 'images_thresholded'
output_paths_all_points = [os.path.join(image_output_dir,"Point"+str(x)) for x in points]
for p in np.unique(output_paths_all_points):
    if not os.path.exists(p):
        os.makedirs(p)

# For each point, create pixel x markers matrix
for i in range(len(points)):
    # Point
    samp = points[i]
    # Path to data
    path = paths_all_points[i]
    # Path to image output
    output_path = output_paths_all_points[i]

    # Initialize np array
    # 1 column for each marker, sample #, x, y
    pixel_mat = np.zeros((num_pixels**2,num_markers+3))

    for j in range(len(markers)):
        mark = markers[j]
        im = io.imread(os.path.join(path,mark+".tif"))
        im_array = np.array(im)
        im_array = im_array.astype(np.float64)

        # Add Gaussian blur
        im_blur = ndimage.gaussian_filter(im_array, sigma=s)

        # Otsu threshold
        thresh = threshold_otsu(im_blur)
        im_thresh = (im_blur > thresh).astype(np.int16)

        # Save image
        im_save = Image.fromarray(im_thresh)
        im_save.save(os.path.join(output_path, mark+".tif"))

        # Turn 1024x1024 -> 1048576x1
        img = im_thresh.flatten()
        pixel_mat[:,j] = img

    # Add column for sample #
    pixel_mat[:,num_markers] = np.repeat(samp, num_pixels**2)

    # Add column for pixel x,y
    pixel_mat[:,num_markers+1] = np.repeat(range(num_pixels),num_pixels)
    pixel_mat[:,num_markers+2] = np.tile(range(num_pixels),num_pixels)

    # Remove all 0s
    pixel_mat_nozeros = pixel_mat[np.sum(pixel_mat[:,:num_markers], axis=1)!=0,:]

    # Write to file
    df = pd.DataFrame(pixel_mat_nozeros)
    df.columns = markers+['sample','x','y']
    feather.write_dataframe(df, os.path.join(output_dir,"Point"+str(samp)+"_sigma"+str(s)+".feather"), compression='uncompressed')

    print(samp)

