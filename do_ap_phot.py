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

from termcolor import colored
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

    # Set up data structure (retrieve if it exists)
    sed_data = {}
    try:
        sed_data = json.load(open(sed_data_loc))
        if target not in sed_data:
            sed_data[target] = {
                        "sed_flux" : {},
                        "sed_freq" : {},
                        "sed_unc_lower" : {},
                        "sed_unc_upper" : {},
                        "sed_telescopenames" : {},
                        "sed_flags" : {}
                }
    except FileNotFoundError:
        sed_data = {target : {
                        "sed_flux" : {},
                        "sed_freq" : {},
                        "sed_unc_lower" : {},
                        "sed_unc_upper" : {},
                        "sed_telescopenames" : {},
                        "sed_flags" : {}
                }
            }
    ## Remove non-NED measurements (to be remeasured shortly)
    bands_to_remove = []
    for band in sed_data[target]["sed_flux"]:
        if band in sed_data[target]["sed_flags"]:
            if "n" not in sed_data[target]["sed_flags"][band]:
                bands_to_remove.append(band)
        else:
            bands_to_remove.append(band)
    for band in bands_to_remove:
        sed_data[target]["sed_flux"].pop(band, None)
        sed_data[target]["sed_freq"].pop(band, None)
        sed_data[target]["sed_unc_lower"].pop(band, None)
        sed_data[target]["sed_unc_upper"].pop(band, None)
        sed_data[target]["sed_telescopenames"].pop(band, None)

    # Get aperture file paths
    ap_file_paths = []
    bg_ap_file_paths = []
    filter_names = []
    for r, d, f in os.walk(workdir+subdir+'apertures/'):
        for file in f:
            if "background" not in file and "Lower" not in file:
                full_path = workdir+subdir+'apertures/'+file
                bg_path = workdir+subdir+'apertures/background'+file
                ap_file_paths.append(full_path)
                bg_ap_file_paths.append(bg_path)
                filter_names.append(file.replace('.reg',''))

    # Iterate through fits files
    for fltr, ap_path, bg_ap_path in zip(filter_names, ap_file_paths, bg_ap_file_paths):
        print(colored(' . . . working on band \'' + fltr + '\' . . . ','green'))
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
        main_ind = 0
        for cutout, color, ap, ind in zip(pix_ap_cutouts, ap_colors, sky_aps, range(0, len(sky_aps))):
            if color=='green':
                eff_radius = (ap.width.value+ap.height.value)/2
                flux_final += integrate_flux(cutout, background)
                main_ind = ind
            elif color=='red':
                flux_final -= integrate_flux(cutout, background)
            else:
                logging.warning('Ambiguity warning ... an aperture was an ambiguous color, so it was not used')

        # Find uncertainties
        stat_unc_upper = imf.calc_unc_background(img, bg_ap)
        stat_unc_lower = imf.calc_unc_background(img, bg_ap)
        if "SPIRE" in fltr\
                or "PACS" in fltr:
            stat_unc_upper = imf.calc_unc_apcopy(img, pix_aps[main_ind], background)
            stat_unc_lower = imf.calc_unc_apcopy(img, pix_aps[main_ind], background)
        abs_unc_upper = flux_final * flux_conversion.abs_unc[fltr]
        abs_unc_lower = flux_final * flux_conversion.abs_unc[fltr]
        total_unc_upper = (stat_unc_upper**2 + abs_unc_upper**2)**0.5
        total_unc_lower = (stat_unc_lower**2 + abs_unc_lower**2)**0.5

        # If upper limit to be found, just do that and continue
        for key in sed_data[target]["sed_flags"]:
            if fltr in key:
                if 'u' in sed_data[target]["sed_flags"][key]:
                    correction = 1.0
                    if "PACS" in fltr:
                        correction = flux_conversion.beam_size[fltr](img_hdr)
                    flux_final = correction * 4.5*np.sqrt(np.mean(pix_ap_cutouts[main_ind]**2))
                    total_unc_upper = -1
                    total_unc_lower = 0
                elif 'a' in sed_data[target]["sed_flags"][key]:
                    flux_final = flux_conversion.beam_size[fltr](img_hdr) * 4.5*np.sqrt(np.mean(pix_ap_cutouts[main_ind]**2))
                    total_unc_upper = -1
                    total_unc_lower = 0

        # If 'r' flag
        for key in sed_data[target]["sed_flags"]:
            if fltr in key:
                if 'r' in sed_data[target]["sed_flags"][key]:
                    try:

                        # Get lower apertures from file
                        lower_sky_aps = Regions.read(workdir+subdir+'apertures/'+fltr+'Lower.reg', format='ds9')
                        ## Convert to pixel aperture
                        lower_pix_aps = [i.to_pixel(wcs) for i in lower_sky_aps]
                        ## Get aperture colors (green=main; red=subtraction)
                        lower_ap_colors = [i.visual['edgecolor'] for i in lower_pix_aps]

                        # Get cutouts of all apertures
                        lower_pix_ap_masks = [i.to_mask() for i in lower_pix_aps]
                        lower_pix_ap_cutouts = [i.multiply(img) for i in lower_pix_ap_masks]

                        # Measure flux in all apertures
                        lower_eff_radius = 1.0
                        lower_flux_final = 0
                        lower_main_ind = 0
                        for cutout,\
                                color,\
                                ap,\
                                ind in zip(lower_pix_ap_cutouts, lower_ap_colors, lower_sky_aps, range(0, len(lower_sky_aps))):
                            if color=='green':
                                lower_eff_radius = (ap.width.value+ap.height.value)/2
                                lower_flux_final += integrate_flux(cutout, background)
                                lower_main_ind = ind
                            elif color=='red':
                                lower_flux_final -= integrate_flux(cutout, background)
                            else:
                                logging.warning('Ambiguity warning ... an aperture was an ambiguous color, so it was not used')

                        # Find lower uncertainties
                        lower_stat_unc_lower = imf.calc_unc_background(img, bg_ap)
                        if "SPIRE" in fltr:
                            lower_stat_unc_lower = imf.calc_unc_apcopy(img, lower_pix_aps[lower_main_ind], background)
                        lower_abs_unc_lower = flux_final * flux_conversion.abs_unc[fltr]
                        lower_total_unc_lower = (lower_stat_unc_lower**2 + lower_abs_unc_lower**2)**0.5

                        # Adjust overall limits
                        total_unc_lower = lower_total_unc_lower + (flux_final - lower_flux_final)

                    except FileNotFoundError:
                        logging.warning('\'r\' keyword was specified, but no \'lower\' fits file is present; keyword has no effect')

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
        fltr_itr = 0
        if 'ALMA' in fltr:
            tele_wl[fltr] = img_hdr['RESTFRQ']
        ## Check if a *distinct* measurement already exists for this filter
        while (fltr+" "+str(fltr_itr)) in sed_data[target]["sed_flux"]:
            fltr_itr += 1
        sed_data[target]["sed_flux"][fltr+" "+str(fltr_itr)] = flux_final
        sed_data[target]["sed_freq"][fltr+" "+str(fltr_itr)] = tele_wl[fltr]
        sed_data[target]["sed_unc_upper"][fltr+" "+str(fltr_itr)] = total_unc_upper
        sed_data[target]["sed_unc_lower"][fltr+" "+str(fltr_itr)] = total_unc_lower
        sed_data[target]["sed_telescopenames"][fltr+" "+str(fltr_itr)] = tele_filters[fltr]

    # Write data to file
    try:  # If the file already exists
        # Open the file for writing
        sed_outfile = open(sed_data_loc, 'w')
        json.dump(sed_data, sed_outfile, indent=5)
    except FileNotFoundError:  # Create a new file if needed
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
    print(colored('Doing photometry on : ' + target,'green'))
    print(' ')
    if target in subdirs:
        full_photometry(target)
    else:
        logging.warning("A target was specified, but no corresponding folder was found, so no photometry was done")
