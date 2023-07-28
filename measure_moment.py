import sys
sys.path.insert(0, './lib')

import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from photutils.isophote import EllipseGeometry, Ellipse
from astropy.io import fits
from astropy.wcs import WCS
import astropy.units as u
from regions import Regions, PixCoord, EllipsePixelRegion
import json
import mgefit.find_galaxy as fgal

import image_functions as imf
import get_ned_data as gnd

import logging

###############################
# Get user input for filepath #
###############################

# Prompt to change working directory
print(' ')
workdir = input("Enter working directory: ")
if workdir[-1] != '/':
    workdir += '/'
pa_data_loc = workdir + 'moment_pa_data.json'

# Prompt to specify target
print(' ')
target = input("Enter the desired target : ")

# Prompt to specify blob number
print(' ')
nblob = 1
nblob = input("Enter blob number (default 1) : ")
if nblob == '' : nblob=1
nblob = int(nblob)

###############
# Import data #
###############

fltr = "W1"

img_path = workdir+target+"/fits/"+fltr+".fits"
ap_path = workdir+target+"/apertures/"+fltr+"_mom.reg"
bg_ap_path = workdir+target+"/apertures/background"+fltr+".reg"

# Open image
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
    ## Get rid of redundant dimensions in data if necessary
    img = img.squeeze()
except FileNotFoundError:
    logging.warning('For '+target_name+', \''+fltr+'.fits\' not found, but there exists a corresponding aperture file')
    pass

# Get apertures
apertures = Regions.read(ap_path, format='ds9')
green_ind = [itr for itr in range(len(apertures)) if apertures[itr].visual['edgecolor']=='green']
red_ind = [itr for itr in range(len(apertures)) if apertures[itr].visual['edgecolor']=='red']
main_ap = apertures[green_ind[0]]

#####################
# Mask contaminants #
#####################

# Mask red apertures
sub_aps_pix = [apertures[ii].to_pixel(wcs) for ii in red_ind]
sub_aps_pix_g = sub_aps_pix.copy()
# for aa in sub_aps_pix_g:
#     aa.width += 200
#     aa.height += 200
mask = sum([aa.to_mask().to_image(img.shape) for aa in sub_aps_pix_g])
img = ma.masked_array(img, mask=mask)
img = img.filled(0)
img = np.nan_to_num(img, nan=0.0)

# Adjust aperture center to flux peak
main_ap_pix = main_ap.to_pixel(wcs)
flux_center = np.unravel_index(np.argmax(main_ap_pix.to_mask().to_image(img.shape)*img), img.shape)
main_ap_pix.center = PixCoord(flux_center[1], flux_center[0])

###################
# Get target flux #
###################

# Get background level
## Get region file
bg_ap = Regions.read(bg_ap_path, format='ds9')[0]
## Get background cutout
bg_cutout = (bg_ap.to_mask()).multiply(img)
## Integrate background flux
background = imf.integrate_flux(bg_cutout, 0.0) / (len(bg_cutout)-2)**2
### If background flux nan, just set to 0.0
if background != background:
    background = 0.0

# Get main flux
total_flux = imf.integrate_flux(main_ap_pix.to_mask().multiply(img), background)

##########################
# Find half-light radius #
##########################

hlap = main_ap_pix.copy()
relflux = 1.0
while relflux >= 0.5:
    hlap.width *= 0.99
    hlap.height *= 0.99
    relflux = imf.integrate_flux(hlap.to_mask().multiply(img), background) / total_flux
hlrad = hlap.width
fullap = hlap.copy()
fullap.width *= 2.5
fullap.height *= 2.5

########################
# Find level at fullap #
########################

# Make a slightly smaller region
fullap_sub = fullap.copy()
fullap_sub.width *= 0.95
fullap_sub.height *= 0.95

# Make cutouts
mask_outer = fullap.to_mask().to_image(np.shape(img))
mask_inner = fullap_sub.to_mask().to_image(np.shape(img))
mask_annulus = mask_outer - mask_inner
annulus_data = mask_annulus*img
annulus_data[annulus_data==0.0] = np.nan
level = np.nanmean(annulus_data)

#######################
# Find PA with moment #
#######################

momfit = fgal.find_galaxy(img, plot=True, quiet=True, level=level, nblob=nblob)
plt.show()

################
# Save PA data #
################

# Open the existing data
pa_data = {}
try:
    pa_data = json.load(open(pa_data_loc))
except FileNotFoundError:
    print("\nPA data file not found ... creating a new one here")
    pass

# Add new data
pa_data[target] = {}
pa_data[target]["pa"] = momfit.pa
pa_data[target]["pa_unc"] = 0.5

# Write
pa_outfile = open(pa_data_loc, 'w')
json.dump(pa_data, pa_outfile, indent=5)
