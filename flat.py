import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from functions import plot_grid, search_names, load_fits

input = search_names("flat", "Input the exposure time identifier for this set of flats: ",
                     "Input the filter identifier for this flat frame: ")

flat_list = input[0]
exposure_time = input[1]
filter_id = input[2]

master_dark = load_fits('master_dark_'+str(exposure_time)+'.fit')

image_data = {}
for image_name in flat_list:
    image_data[image_name] = fits.getdata(image_name) - master_dark

flat_cube = np.stack([image_data[image_name] for image_name in flat_list], axis=0)

# Preview images
plot_grid(flat_cube, flat_list)
plt.show()

master_flat = np.median(flat_cube, axis=0)

print '\nmaster flat median: %.3f counts' % np.median(master_flat)
print 'master flat mean: %.3f counts' % np.mean(master_flat)
print 'master flat max value: %.3f counts' % np.max(master_flat)
print 'master flat min value: %.3f counts' % np.min(master_flat)

normalized_master_flat = master_flat/np.median(master_flat)

# PREVIEW MASTER FLAT
plt.figure(figsize=(15, 15))
plt.imshow(normalized_master_flat, origin='lower', cmap='gray')
plt.title('Normalized Master Flat')
plt.show()

print '\nnormalized master flat median: %.3f' % np.median(normalized_master_flat)
print 'normalized master flat mean: %.3f' % np.mean(normalized_master_flat)
print 'normalized master flat max value: %.3f' % np.max(normalized_master_flat)
print 'normalized master flat min value: %.3f' % np.min(normalized_master_flat)

hdu = fits.PrimaryHDU(normalized_master_flat)
hdu.writeto('master_flat_'+filter_id+'.fit')
print "\nNormalized master flat image saved as 'master_flat_"+filter_id+".fit'."



