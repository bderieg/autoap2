import sys
import os
sys.path.insert(0, './lib')
sys.path.insert(0, './param_files')

import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
from regions import Regions
import json

import flux_conversion
import get_ned_data as gnd
import image_functions as imf

import logging

#################
# Set constants #
#################

workdir = ""

############################
# Flux integration routine #
############################


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
        filters = [(f.name).replace('.fits','') for f in os.scandir(workdir+target+"/fits/")]
        for fltr in filters:
            img_path = workdir+target+"/fits/"+fltr+".fits"
            main_ap_sky, bg_ap_pix = imf.fit_ellipse_with_coordinates(
                    img_path,
                    gnd.get_coords(target),
                    gnd.get_ellipse_parameters(target)["axis_ratio"],
                    gnd.get_ellipse_parameters(target)["position_angle"],
                    imf.find_background_flux(img_path)
                    )
            # Write out apertures to files
            bg_ap_pix.write(workdir+target+"/apertures/background"+fltr+".reg", format='ds9', overwrite=True)
            main_ap_sky.write(workdir+target+"/apertures/"+fltr+".reg", format='ds9', overwrite=True)
    else:
        logging.warning("A target was specified, but no corresponding folder was found, so no photometry was done")

