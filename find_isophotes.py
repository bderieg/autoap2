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

import image_functions as imf

import logging

###############################
# Get user input for filepath #
###############################

# Prompt to change working directory
print(' ')
workdir = input("Enter working directory: ")
if workdir[-1] != '/':
    workdir += '/'

# Prompt to specify target
print(' ')
target = input("Enter the desired target : ")

###############
# Import data #
###############

fltr = "W1"

img_path = workdir+target+"/fits/"+fltr+".fits"
# img_path = '/home/ben/Desktop/research/research_boizelle_working/ap_phot_data/NGC4261/fits/W1.fits'
ap_path = workdir+target+"/apertures/"+fltr+".reg"
# ap_path = '/home/ben/Desktop/research/research_boizelle_working/ap_phot_data/NGC4261/apertures/W1.reg'
bg_ap_path = workdir+target+"/apertures/background"+fltr+".reg"
# bg_ap_path = '/home/ben/Desktop/research/research_boizelle_working/ap_phot_data/NGC4261/apertures/backgroundW1.reg'

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
mask = sum([apertures[ii].to_pixel(wcs).to_mask().to_image(img.shape) for ii in red_ind])
img = ma.masked_array(img, mask=mask)

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
img_sub = img - background

# Get main flux
total_flux = imf.integrate_flux(main_ap.to_pixel(wcs).to_mask().multiply(img), background)

##########################
# Find half-light radius #
##########################

hlap = main_ap.copy().to_pixel(wcs)
relflux = 1.0
while relflux >= 0.5:
    hlap.width *= 0.99
    hlap.height *= 0.99
    relflux = imf.integrate_flux(hlap.to_mask().multiply(img), background) / total_flux
hlrad = hlap.width

################
# Fit isophote #
################

# Create fitting geometry
geom = EllipseGeometry(
        x0=hlap.center.xy[0],
        y0=hlap.center.xy[1],
        sma=hlap.width,
        eps=(1-hlap.height/hlap.width),
        pa=hlap.angle.value * np.pi/180
        )

# Create a bunch of isophotes
fitobj = Ellipse(img_sub, geom)
isofitlist = fitobj.fit_image(sma0=1.5*hlrad, minsma=hlrad, step=0.4*hlrad, linear=True, conver=1e-3, nclip=12, sclip=5.0)

isoaplist = []
for aa in isofitlist:
    try:
        apcolor = 'white'
        if aa.sma/hlrad > 2.5 and aa.sma/hlrad < 3.0:
            apcolor = 'red'
        isoaplist.append(EllipsePixelRegion(
                        PixCoord(aa.x0, aa.y0),
                        aa.sma,
                        aa.sma * (1-aa.eps),
                        angle=aa.pa * 180/np.pi * u.deg,
                        visual={'edgecolor':apcolor}
                    ))
        print(aa.sma/hlrad)
    except ValueError:
        continue

##############
# Show image #
##############

fig,ax = plt.subplots()
cm = plt.get_cmap('cividis').copy()
if np.amax(img) > 0.1:
    cm.set_under('black')
    cm.set_over('black')
    ax.imshow(img_sub, cmap=cm, norm=colors.SymLogNorm(linthresh=0.01, linscale=0.5, vmin=0.1))
else:
    ax.imshow(img, cmap=cm)
# main_ap.to_pixel(wcs).plot()
for aa in isoaplist:
    aa.plot()
plt.show()
