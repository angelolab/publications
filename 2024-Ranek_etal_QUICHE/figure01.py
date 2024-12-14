import imageio as io
from skimage.segmentation import find_boundaries
import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

## figure 1 example overlay
rgb_to_cmyk = np.array([[0.0, 1.0, 1.0],
                         [1.0, 0.0, 1.0],
                         [1.0, 1.0, 0.0]])
cmyk_from_rgb = np.linalg.inv(rgb_to_cmyk)

def rgb_to_cmyk(rgb):
    out = np.dot(rgb, cmyk_from_rgb)
    return out

def cmyk_to_rgb(cmyk):
    return np.dot(cmyk, rgb_to_cmyk)
data_dir = r'/Volumes/Shared/Noah Greenwald/TNBC_Cohorts/SPAIN/image_data/samples'
seg_dir = r'/Volumes/Shared/Noah Greenwald/TNBC_Cohorts/SPAIN/segmentation/samples/deepcell_output'
save_directory = os.path.join('publications', 'figures', 'figure1')
os.makedirs(save_directory)
channels = ['ECAD', 'CD45', 'Fibronectin']
nuclei_channels = ["H3K27me3", "H3K9ac"]
cyan = np.linspace(0, 100, 256)
magenta = np.linspace(0, 100, 256)
yellow = np.linspace(0, 100, 256)
dimgrey = np.linspace(105, 105, 105)
key = np.linspace(0, 100, 256)
cmyk_colormap = np.stack([cyan, magenta, yellow, key], axis=-1) / 256.0
cmyk_cmap = ListedColormap(cmyk_colormap)
for fov in ['TMA32_R2C9', 'TMA32_R8C3', 'TMA33_R2C1', 'TMA33_R2C2']:
    img = None
    pred = None
    inst_seg = np.squeeze(io.imread(os.path.join(seg_dir, fov +'_whole_cell.tiff'))).astype(np.float32)
    boundaries = find_boundaries(inst_seg, mode='inner')
    inst_seg[boundaries] = 0
    inst_seg = (inst_seg > 0).astype(np.float32)
    for channel in channels:
        data_path = os.path.join(data_dir, fov, channel + '.tiff')
        if img is None:
            img = [np.squeeze(io.imread(data_path))]
        else:
            img += [np.squeeze(io.imread(data_path))]

    nuclei_img = []
    for nuclei_chan in nuclei_channels:
        data_path = os.path.join(data_dir, fov, nuclei_chan + '.tiff')
        nuclei_img += [np.squeeze(io.imread(data_path))]
    nuclei_img = np.stack(nuclei_img, axis=-1)
    nuclei_img = nuclei_img / np.quantile(nuclei_img, 0.99)
    nuclei_img = nuclei_img.mean(axis=-1)
    
    img = np.stack(img, axis=-1)
    img = img / np.quantile(img, 0.99, axis=(0,1))
    img = np.clip(img, 0, 1)
    fig, ax = plt.subplots(1,1, figsize=(5,5))
    img = rgb_to_cmyk(img)
    ax.imshow(nuclei_img, cmap="gray", vmin=0, vmax=1.0, interpolation="none")
    ax.imshow(img, vmin=0, vmax=1.0, interpolation="none", alpha=0.8)
    ax.axis("off")
    channel_names = "".join(channels)
    plt.title(fov)
    plt.savefig(os.path.join(save_directory, f'{fov}_overlay.svg'), transparent=True, dpi=300)