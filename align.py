from functions import plot_grid, search_names
from astropy.io import fits
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import interpolation as interp
from skimage.feature.register_translation import (register_translation, _upsampled_dft)
import warnings
warnings.filterwarnings('ignore')

input = search_names('stacked', 'Input the object name: ')

stacked_list = input[0]
name = input[1]

image_data = {}
for image_name in stacked_list:
    image_data[image_name] = fits.getdata(image_name)

stacked_cube = np.stack([image_data[frame] for frame in stacked_list], axis=0)

print "\nRecord the index (first image is 1) of the reference image to use for alignment, and close the figure."
# show the images:
plot_grid(stacked_cube, stacked_list)
plt.show()

selection = raw_input('Type the number of the desired reference image for alignment. (1-'+str(len(stacked_list))+'): ')
zero_shift_image = stacked_list[int(selection)-1]

print "\nShifting images..."
imshifts = {}
for image in stacked_list:
    result, error, diffphase = register_translation(image_data[zero_shift_image], image_data[image], 1000)
    imshifts[image] = result

shifted_science_data = {}
for i in range(len(stacked_list)):
    shifted_science_data[stacked_list[i]] = interp.shift(image_data[stacked_list[i]], imshifts[stacked_list[i]])

if not os.path.exists("./aligned/"):
    os.makedirs("./aligned/")

for filename in stacked_list:
    hdu = fits.PrimaryHDU(shifted_science_data[filename])
    hdu.writeto("./aligned/"+filename)
    print "Aligned science image saved to ./aligned/"+filename+"."
