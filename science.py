from functions import search_names, plot_grid
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from skimage.feature.register_translation import (register_translation, _upsampled_dft)
from scipy.ndimage import interpolation as interp

input = search_names(prompt1="Input the exposure time identifier for these science images: ",
                     prompt2="Input the filter identifier for these science images: ")

science_list = input[0]
exposure_time = input[1]
filter_id = input[2]

hdu_list = fits.open('master_dark_'+exposure_time+'.fit')
master_dark = hdu_list[0].data
hdu_list.close()

hdu_list = fits.open('master_flat_'+filter_id+'.fit')
master_flat = hdu_list[0].data
hdu_list.close()

image_data = {}
for image_name in science_list:
    image_data[image_name] = (fits.getdata(image_name) - master_dark)/master_flat
science_cube = np.stack([image_data[image_name] for image_name in science_list], axis=0)

print "\nRecord the number of the reference image to use for alignment, and close the figure."
plot_grid(science_cube, science_list)
plt.show()

selection = raw_input('Type the number of the desired reference image for alignment. (1-'+str(len(science_list))+'): ')
zero_shift_image = science_list[int(selection)-1]

print "\nShifting and stacking images..."
imshifts = {}
for image in science_list:
    result, error, diffphase = register_translation(
        science_cube[zero_shift_image],
        science_cube[image], 750)
    imshifts[image] = result

shifted_science_data = {}
for i in range(len(science_list)):
    shifted_science_data[science_list[i]] = interp.shift(science_cube[science_list[i]], imshifts[science_list[i]])

science_cube = np.stack(shifted_science_data.values(), axis=0)
science_stacked = np.average(science_cube, axis=0)
# science_stacked = np.median(science_cube, axis=0)

plt.figure(1)
plt.figure(figsize=(15, 15));
plt.title('Aligned and Stacked Science image')
plt.imshow(np.log10(science_stacked), origin='lower', cmap='gray', vmin=1.5, vmax=3)
plt.show()

name = raw_input("Input the name of this object (for saving the final image): ")

hdu = fits.PrimaryHDU(science_stacked)
hdu.writeto('stacked_'+name+'_'+filter_id+'_'+exposure_time+'.fit')

print "\nStacked science image saved as 'stacked_"+name+"_"+filter_id+"_"+exposure_time+".fit'."
