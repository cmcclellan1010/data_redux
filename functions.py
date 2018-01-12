import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gs
from glob import glob
from astropy.io import fits


def plot_grid(datacube, imagenames):
    n_images = len(datacube)
    xplots = int(np.around(np.sqrt(n_images)))
    yplots = xplots + 1
    gridspec = gs.GridSpec(yplots, xplots)
    plt.figure(figsize=(15, 15))
    for i in range(n_images):
        image = datacube[i]
        plt.subplot(gridspec[i])
        plt.imshow(np.log10(image), origin='lower', cmap='gray')
        plt.title(imagenames[i])


def search_names(img_type='', prompt1='', prompt2=''):
    id1 = prompt1
    id2 = prompt2
    if prompt1 != '':
        id1 = raw_input(prompt1)
    if prompt2 != '':
        id2 = raw_input(prompt2)
    filename_list = []
    for filename in list(set().union(glob('*.FIT'), glob('*.fit'))):
        if img_type in filename and '_'+id1 in filename and id2 in filename:
            filename_list.append(filename)
    return filename_list, id1, id2


def load_fits(filename):
    hdu_list = fits.open(filename)
    data = hdu_list[0].data
    hdu_list.close()
    return data

