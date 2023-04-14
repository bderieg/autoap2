from astropy.io import fits
from astropy.wcs import WCS
from reproject import reproject_exact
import matplotlib.pyplot as plt
import numpy as np
from reproject.mosaicking import find_optimal_celestial_wcs, reproject_and_coadd
import os

############################################
# Get user input for which files to stitch #
############################################

# Prompt to change working directory
print(' ')
workdir = input("Enter working directory: ")
if workdir[-1] != '/':
    workdir += '/'

print(' ')
print('Found the following files:')
for f in os.listdir(workdir):
    print('\t'+f)
print(' ')
file1name = input('Enter the first image name:')
print(' ')
file2name = input('Enter the second image name:')
print(' ')
outputfilename = input('Enter the desired name of the stitched image file:')
print(' ')
print('Stitching . . . ')
print(' ')

#######################
# Open two data files #
#######################

# Open file 1
im1 = fits.open(workdir+file1name)
# Loop through extensions to get the one with image data
im1data = im1[0].data  # Default to first extension just in case
wcs1 = WCS(im1[0].header, naxis=[1,2])  # ^^
for ext in im1:
    if 'IMAGE' in ext.header.get('EXTNAME','').upper():
        im1data = ext.data
        # Determine WCS for sky-to-pixel aperture conversion
        wcs1 = WCS(ext.header, naxis=[1,2])
## Get rid of redundant dimensions in data if necessary
im1data = im1data.squeeze()
# Assume desired header is always in the first extension
im1hdr = im1[0].header

# Open file 2
im2 = fits.open(workdir+file2name)
# Loop through extensions to get the one with image data
im2data = im2[0].data  # Default to first extension just in case
wcs2 = WCS(im2[0].header, naxis=[1,2])  # ^^
for ext in im2:
    if 'IMAGE' in ext.header.get('EXTNAME','').upper():
        im2data = ext.data
        # Determine WCS for sky-to-pixel aperture conversion
        wcs2 = WCS(ext.header, naxis=[1,2])
## Get rid of redundant dimensions in data if necessary
im2data = im2data.squeeze()
# Assume desired header is always in the first extension
im2hdr = im2[0].header

####################
# Stitch and write #
####################

# Find new WCS
wcs_out, shape_out = find_optimal_celestial_wcs([(im1data,wcs1),(im2data,wcs2)])

# Make mosaic
mosaic, _ = reproject_and_coadd([(im1data,wcs1),(im2data,wcs2)], wcs_out, shape_out=shape_out, reproject_function=reproject_exact)

# Edit the header
outhdr_bare = wcs_out.to_header()
outhdr = im1hdr.copy()
outhdr.update(outhdr_bare)


# Write
fits.writeto(workdir+outputfilename, mosaic, outhdr, overwrite=True)
