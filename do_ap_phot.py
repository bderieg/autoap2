import sys
import os
sys.path.insert(0, './lib')
sys.path.insert(0, './param_files')

import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
from regions import Regions
import json

# import image_functions as imf
# import get_ned_data as gnd
import flux_conversion

import logging

#################
# Set constants #
#################

workdir = '/home/ben/Desktop/research/research_boizelle_working/ap_phot_data/'
subdir = 'GAMA177186/'
target_name = 'GAMA177186'
sed_data_loc = '/home/ben/Desktop/research/research_boizelle_working/sed_data.json'

############################
# Flux integration routine #
############################

def integrate_flux(cutout, background):
    flux = 0

    for row in cutout:
        for pix_val in row:
            if pix_val==pix_val and pix_val!=0.0:
                flux += pix_val
                flux -= background

    return flux

##################################
# Do photometry on each aperture #
##################################

# Import necessary static data
tele_wl = json.load(open('./param_files/telescope_wavelengths.json'))
tele_filters = json.load(open('./param_files/telescope_filter_names.json'))

# Set up SED data structure (dict with wavelength/flux pairs)
sed_data = {}
sed_flags = {}
sed_filternames = {}
sed_telescopenames = {}

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
        img_hdr = fits.open(fits_path)[0].header
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

    # Get cutouts of all apertures
    pix_ap_masks = [i.to_mask() for i in pix_aps]
    pix_ap_cutouts = [i.multiply(img) for i in pix_ap_masks]

    # Measure flux in all apertures
    background = 0
    eff_radius = 1.0
    flux_final = 0
    for cutout, color in zip(pix_ap_cutouts, ap_colors):
        if color=='green':
            flux_final += integrate_flux(cutout, background)
        elif color=='red':
            flux_final -= integrate_flux(cutout, background)
        else:
            logging.warning('Ambiguity warning ... an aperture was an ambiguous color, so it was not used')

    # Convert to units of Jy
    flux_final = flux_final *\
            flux_conversion.brightness_conversion[fltr](img_hdr) *\
            (1/flux_conversion.beam_size[fltr](img_hdr)) *\
            flux_conversion.pix_size[fltr](img_hdr) *\
            flux_conversion.color_correction[fltr](img_hdr) *\
            flux_conversion.other_correction[fltr](img_hdr, eff_radius)

    sed_data[tele_wl[fltr]] = flux_final
    sed_filternames[tele_wl[fltr]] = fltr
    sed_telescopenames[tele_wl[fltr]] = tele_filters[fltr]

#####################
# Write SED to file #
#####################

sed_full = {target_name : {
        "sed_data" : sed_data,
        "sed_telescopenames" : sed_telescopenames,
        "sed_filternames" : sed_filternames,
        "sed_flags" : sed_flags
}}

try:
    all_sed_data = json.load(open(sed_data_loc))
    sed_outfile = open(sed_data_loc, 'w')
    all_sed_data |= sed_full
    json.dump(all_sed_data, sed_outfile, indent=5)
except FileNotFoundError:
    print("\nSED data file not found ... creating a new one here")
    sed_outfile = open(sed_data_loc, 'w')
    json.dump(sed_full, sed_outfile, indent=5)
