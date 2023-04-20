import sys
import os
sys.path.insert(0, './lib')
sys.path.insert(0, './param_files')

import matplotlib.pyplot as plt
from matplotlib import colors
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from regions import Regions
import json

import flux_conversion
import image_functions as imf

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

    # Set up SED data structure (dict with wavelength/flux(/other) pairs)
    sed_data = {}
    sed_unc_upper = {}
    sed_unc_lower = {}
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
            ## Get rid of redundant dimensions in data if necessary
            img = img.squeeze()
            # Assume desired header is always in the first extension
            img_hdr = fitsfile[0].header
        except FileNotFoundError:
            logging.warning('For '+target_name+', \''+fltr+'.fits\' not found, but there exists a corresponding aperture file')
            continue

        # Get flags if possible
        sed_flags = {}
        try:
            all_sed_data_temp = json.load(open(sed_data_loc))
            if target in all_sed_data_temp:
                if 'sed_flags' in all_sed_data_temp[target]:
                    sed_flags = all_sed_data_temp[target]['sed_flags']
        except FileNotFoundError:
            pass

        # Get background level
        ## Get region file
        bg_ap = Regions.read(bg_ap_path, format='ds9')[0]
        ## Get background cutout
        bg_cutout = (bg_ap.to_mask()).multiply(img)
        ## Integrate background flux
        background = integrate_flux(bg_cutout, 0.0) / (len(bg_cutout)-2)**2
        ### If background flux nan, just set to 0.0
        if background != background:
            background = 0.0

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

        # Find uncertainties
        stat_unc_upper = imf.calc_unc_background(img, bg_ap)
        stat_unc_lower = imf.calc_unc_background(img, bg_ap)
        if "SPIRE" in fltr:
            # TODO: Make sure pix_aps[0] here is actually the desired main aperture
            stat_unc_upper = imf.calc_unc_apcopy(img, pix_aps[0], background)
            stat_unc_lower = imf.calc_unc_apcopy(img, pix_aps[0], background)
        abs_unc_upper = flux_final * flux_conversion.abs_unc[fltr]
        abs_unc_lower = flux_final * flux_conversion.abs_unc[fltr]
        total_unc_upper = (stat_unc_upper**2 + abs_unc_upper**2)**0.5
        total_unc_lower = (stat_unc_lower**2 + abs_unc_lower**2)**0.5

        # If upper limit to be found, just do that and continue
        # TODO: Set a flag for "no background" and set the background to 0
        if fltr in sed_flags:
            if 'u' in sed_flags[fltr]:
                correction = 1.0
                if "PACS" in fltr:
                    correction = flux_conversion.beam_size[fltr](header)
                upper_limit = correction * 4.5*np.sqrt(np.mean(pix_ap_cutouts[0]**2))
                total_unc_upper = -1
                total_unc_lower = 0

        # Convert to units of Jy
        flux_final = flux_final *\
                flux_conversion.brightness_conversion[fltr](img_hdr) *\
                (1/flux_conversion.beam_size[fltr](img_hdr)) *\
                flux_conversion.pix_size[fltr](img_hdr) *\
                flux_conversion.color_correction[fltr](img_hdr) *\
                flux_conversion.other_correction[fltr](img_hdr, eff_radius)
        total_unc_upper = total_unc_upper *\
                flux_conversion.brightness_conversion[fltr](img_hdr) *\
                (1/flux_conversion.beam_size[fltr](img_hdr)) *\
                flux_conversion.pix_size[fltr](img_hdr) *\
                flux_conversion.color_correction[fltr](img_hdr) *\
                flux_conversion.other_correction[fltr](img_hdr, eff_radius)
        total_unc_lower = total_unc_lower *\
                flux_conversion.brightness_conversion[fltr](img_hdr) *\
                (1/flux_conversion.beam_size[fltr](img_hdr)) *\
                flux_conversion.pix_size[fltr](img_hdr) *\
                flux_conversion.color_correction[fltr](img_hdr) *\
                flux_conversion.other_correction[fltr](img_hdr, eff_radius)

        # Add data to SED structure
        ## If ALMA image, get frequency from image
        if 'ALMA' in fltr:
            tele_wl[fltr] = img_hdr['RESTFRQ']
        ## Check if a *distinct* measurement already exists at this wavelength
        while tele_wl[fltr] in sed_filternames and sed_filternames[tele_wl[fltr]] != fltr:
            tele_wl[fltr] += 1  # If so, increment the wavelength by one to distinguish (but keep essentially the same)
        sed_data[tele_wl[fltr]] = flux_final
        sed_unc_upper[tele_wl[fltr]] = total_unc_upper
        sed_unc_lower[tele_wl[fltr]] = total_unc_lower
        sed_filternames[tele_wl[fltr]] = fltr
        sed_telescopenames[tele_wl[fltr]] = tele_filters[fltr]

    # Write data to SED file
    sed_full = {target_name : {
            "sed_data" : sed_data,
            "sed_unc_upper" : sed_unc_upper,
            "sed_unc_lower" : sed_unc_lower,
            "sed_telescopenames" : sed_telescopenames,
            "sed_filternames" : sed_filternames,
            "sed_flags" : sed_flags
    }}

    # TODO: There's a problem with the logic here . . . just redo it probably
    try:
        all_sed_data = json.load(open(sed_data_loc))
        sed_outfile = open(sed_data_loc, 'w')
        if target in all_sed_data:
            if "sed_flags" in all_sed_data[target]:
                sed_full[target]["sed_flags"] |= all_sed_data[target]["sed_flags"]
            nfreqs = []
            for freq in all_sed_data[target]["sed_filternames"]:
                if "n" in all_sed_data[target]["sed_flags"][all_sed_data[target]["sed_filternames"][freq]]:
                    nfreqs.append(float(freq))
                    sed_full[target]["sed_flags"][all_sed_data[target]["sed_filternames"][freq]] = "n"
            for freq in nfreqs:
                origfreq = freq
                while freq in sed_full[target]["sed_data"]:
                    freq += 1.0
                sed_full[target]["sed_data"][str(freq)] = all_sed_data[target]["sed_data"][str(origfreq)]
                sed_full[target]["sed_unc_upper"][str(freq)] = all_sed_data[target]["sed_unc_upper"][str(origfreq)]
                sed_full[target]["sed_unc_lower"][str(freq)] = all_sed_data[target]["sed_unc_lower"][str(origfreq)]
                sed_full[target]["sed_telescopenames"][str(freq)] = all_sed_data[target]["sed_telescopenames"][str(origfreq)]
                sed_full[target]["sed_filternames"][str(freq)] = all_sed_data[target]["sed_filternames"][str(origfreq)]
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
    print(' ')
    print('Doing photometry on : ' + target)
    print(' ')
    if target in subdirs:
        full_photometry(target)
    else:
        logging.warning("A target was specified, but no corresponding folder was found, so no photometry was done")
