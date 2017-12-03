# data_redux
Code for dark-subtracting, flat-dividing, aligning, and stacking single-filter fits files for color image synthesis in DS9.

### Known Bugs
* Faint images, especially in blue filter, have trouble aligning properly. Making the upsampling factor lower appears to fix this issue, with some downsides. Adaptive upsampling should be incorporated as a permanent fix.
* When displaying the image grid, some science images display pure static, or mangled versions of the images. This is a display issue and does not mean anything is wrong with the files themselves. The final image shows up fine, without any static.

# Using the Software

## File Naming Conventions
For this software to work properly, each image must contain certain tags to uniquely specify the file's exposure time, filter, etc. The exact format of each of the tags does not matter, as long as they are consistent throughout an entire observing run and they are separated by underscores.

**Darks** must contain the word "dark", and an exposure time tag.

**Flats** must contain the word "flat", an exposure time tag, and a filter tag.

**Science images** must contain an exposure time tag and a filter tag. The target name is not required for searching the directory, but you will have to specify one when the stacked image is saved.

Image type | Example of proper tagging
-------------------|-------------------
Dark | dark_30s_001.fit
Flat | flat_15s_r_005.FIT
Science | albireo_700ms_b_004.fit

## 1. Setup
Put your science images in each filter, flats in each filter, and darks matching the exposure times of both in the same folder. Download the .ZIP file from this repository, and extract it into the same folder as your images. The code should be scattered around in the directory now, among all your flats, darks, and science images.

## 2. Master Darks
Each exposure time you use must have a corresponding master dark of the same duration. You also need master darks for each of your flat frames, which have their own exposure times.

1. Open terminal in the directory where you put all your images
2. Type "python dark.py" in the terminal.
3. Type in the exposure time tag for this set of darks. Ex: 2s
4. Look at the figures. Make sure the number of images you were expecting show up.
5. Check to make sure the "master_dark_[exposure time]" file saved properly in the directory.
6. Repeat steps 2-5 for all the different exposure times you'll need master darks for.

## 3. Master Flats
For every filter you take images in, you will need a flat field in that filter. Flat fields must be dark-subtracted, so you will also need a master dark matching the exposure time of each flat. 

1. Type “python flat.py” in the terminal.
2. Type in the exposure time tag for this set of flats (ex: 2s)
3. Type in the filter tag (ex: r)
4. Look at the figures, close them if everything looks good.
5. Master flat will be saved as “master_flat_[filter].fit”
6. Repeat steps 1-5 for each filter you have.

## 4. Science Images
You'll need to line up all the images you have in each filter so you can combine the light by stacking. Each of these images must be dark-subtracted and flat-divided before stacking.

1. Type "python science.py" in the terminal.
2. Type in the exposure time tag for this set of science images (ex: 30s)
3. Type in the filter tag (ex: r)
4. Look at the figure with all the science frames, and pick one where the stars are bright and the field looks relatively even. Count from the top left across the first row, 1, 2, 3, 4... continuing with the next row until you reach the image you picked.
5. Close the figure, and type in the number of the image you picked for alignment. Press enter.
6. A final stacked image will be presented. Check that there aren't any diagonal patterns of repeated stars--if there are, the registration algorithm failed.
7. Name the target something simple, ex: “cig”
8. Close the figure. The stacked science image will be saved as “stacked_[filename].fit”. 
9. Repeat steps 1-8 for each filter.

## 5. Aligning Color Images
Before combining your stacked images in each filter into a final color image, you need to align them. This part of the program translates your stacked images into alignment, and then saves them with the prefix "aligned_...". These images can be loaded into an RGB frame in DS9 to produce a color image.

1. Type “python align.py”
2. Type the object name you used for all three stacked images
3. Choose an image to use for alignment
4. Type in the number of the image you chose, and press enter.
5. Three aligned images will be saved as "aligned_[filename]".
6. Load images into an RGB frame in DS9 to produce a color image



