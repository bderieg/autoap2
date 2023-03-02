import sys
import os
sys.path.insert(0, './lib')

import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
from regions import Regions

# import image_functions as imf
# import get_ned_data as gnd

import logging

#################
# Set constants #
#################

workdir = '/home/ben/Desktop/research/research_boizelle_working/ap_phot_data/'
subdir = 'GAMA177186/'

##################################################
# Do photometry on each aperture and add to file #
##################################################

# Get aperture file paths
ap_file_paths = []
filter_names = []
for r, d, f in os.walk(workdir+subdir+'apertures/'):
    for file in f:
        full_path = workdir+subdir+'apertures/'+file
        ap_file_paths.append(full_path)
        filter_names.append(file)

# Iterate through fits files
for fltr, ap_path in zip(filter_names, ap_file_paths):

    # Get path for fits file corresponding to aperture
    fits_path = workdir+subdir+'fits/'+fltr+'.fits'

    # Attempt to open fits file
    try:
        img = fits.open(fits_path)[0].data
        # Determine WCS for sky-to-pixel aperture conversion
        wcs = WCS(fits.open(fits_path)[0].header, naxis=[1,2])
    except FileNotFoundError:
        logging.warning('\''+fltr+'.fits\' not found, but there exists a corresponding aperture file')
        continue

    # Get apertures from file
    sky_aps = Regions.read(ap_path, format='ds9')
    ## Convert to pixel aperture
    pix_aps = [i.to_pixel(wcs) for i in sky_aps]
    ## Get aperture colors (green=main; red=subtraction)
    ap_colors = [i.visual['edgecolor'] for i in pix_aps]
