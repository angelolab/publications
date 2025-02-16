"""This script takes the output of the buffer polygons from ArcGIS and does two things:

1. Removes the padding
2. Standardizes the buffer numbers/colors

The resulting greyscale masks are exported for downstream processing.

"""

import matplotlib.pyplot as plt
import os
import skimage.io as io
import cv2
import scipy.spatial.distance
from skimage.measure import label, regionprops_table
import pandas as pd
import numpy as np

# Define directories
base_dir = '/Volumes/T7 Shield/MIBI_data/NHP_TB_Cohort/Panel2'
polygon_dir = os.path.join(base_dir, 'masks', 'polygons')
buffer_dir = '/Users/erinmccaffrey/Library/CloudStorage/GoogleDrive-erinmcc@stanford.edu/My ' \
             'Drive/Grad_School/AngeloLab/MIBIProjects/NHP_TB/Cohort/downstream_analysis/TB_Granulomas_GIS/output' \
             '/buffer_polygons/rescaled'
results_dir = os.path.join(base_dir, 'masks', 'buffers')

# Get list of samples to process
# samples = io_utils.list_files(buffer_dir, substrs=['sample'])
# cohort_samples = [sub.replace('_buffers.tif', '') for sub in samples]

# define the padding parameters for trimming
max_x = 8192
max_y = 10240

# define minimum buffer size
buffer_area_min = 15000
buffer_area_max = 4262421

# iterate over samples

sample = 'sample52'

print("Working on: " + sample)

# read in the buffer mask and the original polygon
poly_mask = io.imread(os.path.join(polygon_dir, sample + '_poly.tif'))
poly_padded_mask = (io.imread(os.path.join(polygon_dir, sample + '_poly_padded.tif'))).astype(np.uint8)
buffers_padded = io.imread(os.path.join(buffer_dir, sample + '_buffers.tif'))

# convert the buffer image to hsv
hsv = cv2.cvtColor(buffers_padded, cv2.COLOR_BGR2HSV)
h, s, v = cv2.split(hsv)
plt.imshow(v)

# threshold value image and invert
buffer_mask = cv2.threshold(v, 128, 255, cv2.THRESH_BINARY)[1]
plt.imshow(buffer_mask)

# # mask the h image to get the buffers
# buffer_mask = s > 0
# plt.imshow(buffer_mask)

# zero out necrotic core and non-buffered areas
zero_mask = poly_padded_mask == 0
buffer_mask = buffer_mask * (zero_mask.astype(int))
buffer_mask = buffer_mask.astype(np.uint8)
plt.imshow(buffer_mask)

# convert to label image
buffers_label_img = label(buffer_mask)
plt.imshow(buffers_label_img, 'jet')
plt.colorbar()

# remove objects outside of buffer size parameters
buffer_props = pd.DataFrame(regionprops_table(buffers_label_img, properties=('label', 'area')))
buffer_props
large_objects = buffer_props.loc[buffer_props['area'] > buffer_area_max, 'label'].tolist()
small_objects = buffer_props.loc[buffer_props['area'] < buffer_area_min, 'label'].tolist()
other_objects = [47]
remove_objects = large_objects + small_objects + other_objects
for remove_label in remove_objects:
    buffers_label_img[np.where(buffers_label_img == remove_label)] = 0
plt.imshow(buffers_label_img)

# re-label to buffer number

# first get the perimeter of the necrosis
contours1, hierarchy1 = cv2.findContours(poly_padded_mask, cv2.RETR_TREE, cv2.CHAIN_APPROX_SIMPLE)
output1 = np.zeros(poly_padded_mask.shape, dtype=np.uint8)
cv2.drawContours(output1, contours1, -1, (255, 255, 255), 3)
binary_perim1 = (output1 / 255).astype(int)
necrosis_props = regionprops_table(binary_perim1, properties=('centroid', 'coords'))
perim_coords = necrosis_props['coords'][0]

# get buffer labels
buffer_labels = np.unique(buffers_label_img)
buffer_labels = list(filter(lambda num: num != 0, buffer_labels))

# create array to track distances for relabeling
distance_data = pd.DataFrame(columns=['label', 'distance'])

# iterate to get a table of the buffer labels and min distance
for label in buffer_labels:

    # generate mask of buffer
    label_mask = (buffers_label_img == label).astype(np.uint8)

    # get the area coordinates
    buffer_props = regionprops_table(label_mask, properties=('centroid', 'coords'))
    buffer_coords = buffer_props['coords'][0]

    # subsample just 1% of border coords to speed computation
    N = int(0.01*(buffer_coords.shape[0]))
    buffer_coords_sub = buffer_coords[np.random.choice(buffer_coords.shape[0], N, replace=False)]

    # get the shortest distance
    distances = scipy.spatial.distance.cdist(perim_coords, buffer_coords_sub)
    min_distance = np.min(distances)

    # store results
    label_dist_data = {'label': label, 'distance': min_distance}
    distance_data = pd.concat([distance_data, pd.DataFrame(label_dist_data, index=[0])])

# sort the values and add rank
distance_data = distance_data.sort_values(by=['distance'])
distance_data['rank'] = distance_data['distance'].rank(ascending=True, method='first').astype('int32')

distance_data['buffer'] = distance_data['rank']

# for multifocal: group further
distance_data.loc[distance_data['rank'] <= 5, 'buffer'] = 1
distance_data.loc[(distance_data['rank'] >= 6) & (distance_data['rank'] <= 10), 'buffer'] = 2
distance_data.loc[(distance_data['rank'] >= 11) & (distance_data['rank'] <= 15), 'buffer'] = 3
distance_data.loc[(distance_data['rank'] >= 16) & (distance_data['rank'] <= 20), 'buffer'] = 4
distance_data.loc[(distance_data['rank'] >= 21) & (distance_data['rank'] <= 24), 'buffer'] = 5
distance_data.loc[(distance_data['rank'] >= 25) & (distance_data['rank'] <= 28), 'buffer'] = 6
distance_data.loc[(distance_data['rank'] >= 29) & (distance_data['rank'] <= 32), 'buffer'] = 7
distance_data.loc[(distance_data['rank'] >= 33) & (distance_data['rank'] <= 36), 'buffer'] = 8
distance_data.loc[(distance_data['rank'] >= 37) & (distance_data['rank'] <= 40), 'buffer'] = 9
distance_data.loc[(distance_data['rank'] >= 41) & (distance_data['rank'] <= 43), 'buffer'] = 10
distance_data.loc[(distance_data['rank'] >= 44) & (distance_data['rank'] <= 46), 'buffer'] = 11
distance_data.loc[distance_data['rank'] == 47, 'buffer'] = 12
distance_data.loc[distance_data['rank'] == 48, 'buffer'] = 13
distance_data.loc[distance_data['rank'] == 49, 'buffer'] = 14
distance_data.loc[distance_data['rank'] == 50, 'buffer'] = 15
distance_data.loc[distance_data['rank'] == 51, 'buffer'] = 16

# re-label
relabeled = np.copy(buffers_label_img)
for label in buffer_labels:
    buffer_num = distance_data.loc[distance_data['label'] == label, 'buffer'].item()
    relabeled[np.where(buffers_label_img == label)] = buffer_num
plt.imshow(relabeled)

# for cropping get the actual size and the difference
sample_x = poly_mask.shape[0]
sample_y = poly_mask.shape[1]
diff_x = (max_x - sample_x) / 2
diff_y = (max_y - sample_y) / 2

# define the crop parameters
x_start = int(diff_x)
x_end = int(x_start + sample_x)
y_start = int(diff_y)
y_end = int(y_start + sample_y)

# crop
label_im_cropped = relabeled[x_start:x_end, y_start:y_end]
plt.imshow(label_im_cropped)

# save
io.imsave(os.path.join(results_dir, sample + '_buffers.png'), label_im_cropped.astype(np.uint8), check_contrast=False)







