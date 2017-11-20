# data_redux
Code for dark-subtracting, flat-dividing, aligning, and stacking single-filter fits files for color image synthesis in DS9.

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
2. Type "python dark.py"
3. Type in the exposure time tag for this set of darks. Ex: 2s
4. Look at the figures. Make sure the number of images you were expecting show up.
5. Check to make sure the "master_dark_[exposure time]" file saved properly in the directory.
6. Repeat steps 2-5 for all the different exposure times you'll need master darks for.
