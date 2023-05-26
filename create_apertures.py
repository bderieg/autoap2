import sys
import os
sys.path.insert(0, './lib')
sys.path.insert(0, './param_files')

import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
from regions import Regions
from termcolor import colored
import json

import flux_conversion
import get_ned_data as gnd
import image_functions as imf

import logging

#################
# Set constants #
#################

workdir = ""

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
workdir = input("Enter working directory (relative or absolute) : ")
if workdir[-1] != '/' : workdir += '/'

# Prompt to specify targets for photometry
print(' ')
curinp = input("Create apertures for ALL subfolders in this directory? (y/n) : ")
subdirs = [f.name for f in os.scandir(workdir) if f.is_dir()]
if curinp=='y' or curinp=='Y':
    target_names = subdirs
else:
    print(' ')
    target_names = input("Enter (as a comma-separated list) the desired targets : ")
    target_names = target_names.split(',')

# Prompt to only make apertures for previously-non-existent apertures
only_new = True
print(' ')
curinp = input("Create only previously-non-existent apertures? (Y/n) : ")
if curinp.lower() == 'n' : only_new = False

########################
# Loop through targets #
########################

for target in target_names:
    if target in subdirs:
        filters = [(f.name).replace('.fits','') for f in os.scandir(workdir+target+"/fits/")]
        apertures = [(f.name).replace('.reg','') for f in os.scandir(workdir+target+"/apertures/")]
        if only_new : filters = [ fltr for fltr in filters if fltr not in apertures ]

        # Loop through images for this target
        for fltr in filters:

            # Describe to terminal
            print(' ')
            print(' ')
            print(colored('########## ' + target + ' ' + fltr + ' ##########','green'))
            print(' ')

            # Prompt to proceed or not
            print(' ')
            tempinp = input('Proceed with this image? (Y/n) : ')
            if tempinp.lower() == 'n' : continue
            print(' ')

            # Open image
            img_path = workdir+target+"/fits/"+fltr+".fits"
            try:
                # Open
                fitsfile = fits.open(img_path)
                # Loop through extensions to get the one with image data
                img = fitsfile[0].data  # Default to first extension just in case
                wcs = WCS(fitsfile[0].header, naxis=[1,2])  # ^^
                for ext in fitsfile:
                    if 'IMAGE' in ext.header.get('EXTNAME','').upper():
                        # Get image data
                        img = ext.data
                        # Get WCS data
                        wcs = WCS(ext.header, naxis=[1,2])
            except FileNotFoundError:
                logging.warning('For '+target_name+', \''+fltr+'.fits\' not found, but there exists a corresponding aperture file')
                continue

            # Algorithmically determine main aperture for image
            main_ap_sky, bg_ap_pix = imf.fit_ellipse_with_coordinates(
                    img_path,
                    gnd.get_coords(target),
                    gnd.get_ellipse_parameters(target)["axis_ratio"],
                    imf.find_background_flux(img_path)
                    )
            main_ap_pix = main_ap_sky.to_pixel(wcs)

            # Algorithmically determine subtraction apertures
            sub_aps_pix = imf.find_blobs(img, main_ap_pix.to_mask(), imf.integrate_flux((bg_ap_pix.to_mask()).multiply(img),0)/100)
            sub_aps_sky = [ap.to_sky(wcs) for ap in sub_aps_pix]

            # Write out apertures to files
            main_sub_aps = Regions([main_ap_sky] + sub_aps_sky)
            bg_ap_pix.write(workdir+target+"/apertures/background"+fltr+".reg", format='ds9', overwrite=True)
            main_sub_aps.write(workdir+target+"/apertures/"+fltr+".reg", format='ds9', overwrite=True)


    else:  # If a given image doesn't exist
        logging.warning("A target was specified, but no corresponding folder was found, so nothing was done.")
