import sys
import os
sys.path.insert(0, './lib')
sys.path.insert(0, './param_files')

import matplotlib.pyplot as plt
from matplotlib import colors
from astropy.io import fits
from astropy.wcs import WCS
from regions import Regions
import json

import flux_conversion

import logging

#################
# Set constants #
#################

workdir = ''
target_names = []
sed_data_loc = ''

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

####################################
# Full aperture photometry routine #
####################################

def full_photometry(target_name):
    subdir = target_name + "/"

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
    bg_ap_file_paths = []
    filter_names = []
    for r, d, f in os.walk(workdir+subdir+'apertures/'):
        for file in f:
            if "background" not in file:
                full_path = workdir+subdir+'apertures/'+file
                bg_path = workdir+subdir+'apertures/background'+file
                ap_file_paths.append(full_path)
                bg_ap_file_paths.append(bg_path)
                filter_names.append(file.replace('.reg',''))

    # Iterate through fits files
    for fltr, ap_path, bg_ap_path in zip(filter_names, ap_file_paths, bg_ap_file_paths):
        # Get path for fits file corresponding to aperture
        fits_path = workdir+subdir+'fits/'+fltr+'.fits'

        # Attempt to open fits file
        try:
            # Open
            fitsfile = fits.open(fits_path)
            # Loop through extensions to get the one with image data
            img = fitsfile[0].data  # Default to first extension just in case
            wcs = WCS(fitsfile[0].header, naxis=[1,2])  # ^^
            for ext in fitsfile:
                if 'IMAGE' in ext.header.get('EXTNAME','').upper():
                    img = ext.data
                    # Determine WCS for sky-to-pixel aperture conversion
                    wcs = WCS(ext.header, naxis=[1,2])
            # Assume desired header is always in the first extension
            img_hdr = fitsfile[0].header
        except FileNotFoundError:
            logging.warning('For '+target_name+', \''+fltr+'.fits\' not found, but there exists a corresponding aperture file')
            continue

        # Get background level
        ## Get region file
        bg_ap = Regions.read(bg_ap_path, format='ds9')[0]
        ## Get background cutout
        bg_cutout = (bg_ap.to_mask()).multiply(img)
        ## Integrate background flux
        background = integrate_flux(bg_cutout, 0.0) / (len(bg_cutout)-2)**2

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
        eff_radius = 1.0
        flux_final = 0
        for cutout, color, ap in zip(pix_ap_cutouts, ap_colors, sky_aps):
            if color=='green':
                eff_radius = (ap.width.value+ap.height.value)/2
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

        # Add data to SED structure
        sed_data[tele_wl[fltr]] = flux_final
        sed_filternames[tele_wl[fltr]] = fltr
        sed_telescopenames[tele_wl[fltr]] = tele_filters[fltr]

    # Write data to SED file
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

#######################
#######################
## Start main script ##
#######################
#######################

#################################
# Set constants with user input #
#################################

# Prompt to change working directory
print(' ')
workdir = input("Enter working directory: ")
if workdir[-1] != '/':
    workdir += '/'
sed_data_loc = workdir + "sed_data.json"

# Prompt to specify targets for photometry
print(' ')
curinp = input("Do photometry on ALL subfolders in this directory? (y/n) : ")
subdirs = [f.name for f in os.scandir(workdir) if f.is_dir()]
if curinp=='y' or curinp=='Y':
    target_names = subdirs
else:
    print(' ')
    target_names = input("Enter (as a comma-separated list) the desired targets : ")
    target_names = target_names.split(',')

for target in target_names:
    if target in subdirs:
        full_photometry(target)
    else:
        logging.warning("A target was specified, but no corresponding folder was found, so no photometry was done")
