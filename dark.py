import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from functions import plot_grid, search_names

input = search_names("dark", "Input the exposure time identifier for this set of darks: ")
dark_list, exposure_time = input[0], input[1]

raw_image_data = {}
for image_name in dark_list:
    raw_image_data[image_name] = fits.getdata(image_name)

dark_cube = np.stack([raw_image_data[dark_frame] for dark_frame in dark_list], axis=0)

# Preview images
plot_grid(dark_cube, dark_list)
plt.show()

# master_dark = np.average(dark_cube, axis=0)    # to combine with an average
master_dark = np.median(dark_cube, axis=0)    # to combine with a median

# Preview master dark
plt.figure(figsize=(15, 15))
plt.imshow(np.log10(master_dark), origin='lower', cmap='gray')
plt.title('Master Dark')
plt.show()

hdu = fits.PrimaryHDU(master_dark)
hdu.writeto('master_dark_'+exposure_time+'.fit')
print "\nMaster dark image saved as 'master_dark_"+exposure_time+".fit'."
