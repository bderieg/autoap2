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

###############
# Import data #
###############

# img_path = workdir+target+"/fits/"+fltr+".fits"
img_path = '/home/ben/Desktop/research/research_boizelle_working/ap_phot_data/NGC4261/fits/W1.fits'
# ap_path = workdir+target+"/apertures/"+fltr+".reg"
ap_path = '/home/ben/Desktop/research/research_boizelle_working/ap_phot_data/NGC4261/apertures/W1.reg'
bg_ap_path = '/home/ben/Desktop/research/research_boizelle_working/ap_phot_data/NGC4261/apertures/backgroundW1.reg'

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

# Get main flux
total_flux = imf.integrate_flux(main_ap.to_pixel(wcs).to_mask().multiply(img), background)

################
# Fit isophote #
################

# Initial guess (arbitrary)
hlap = main_ap.copy()
hlap.height /= 10
hlap.width /= 10

# Create fitting geometry
hlap_pix = hlap.to_pixel(wcs)
geom = EllipseGeometry(
        x0=hlap_pix.center.xy[0],
        y0=hlap_pix.center.xy[1],
        sma=hlap_pix.height,
        eps=1-hlap_pix.width/hlap_pix.height,
        pa=hlap_pix.angle.value * np.pi/180
        )

# Create fitting object and fit
fitobj = Ellipse(img, geom)
isofitlist = fitobj.fit_image()

isoaplist = []
hlhit = False
endhit = False
hlrad = 1e9  # arbitrarily high initial value
for ap in isofitlist:
    try:
        curap = EllipsePixelRegion(
                        PixCoord(ap.x0, ap.y0),
                        ap.sma,
                        ap.sma * (1-ap.eps),
                        angle=ap.pa * 180/np.pi * u.deg,
                        visual={'edgecolor':'white'}
                    )
        relflux = imf.integrate_flux(curap.to_mask().multiply(img), background) / total_flux
        if relflux >= 0.5 and hlhit is False:
            hlrad = curap.width
            curap.visual['edgecolor'] = 'orange'
            hlhit = True
        if curap.width >= 2.5*hlrad and hlhit and endhit is False:
            curap.visual['edgecolor'] = 'red'
            endhit = True
        isoaplist.append(curap)
    except ValueError:
        continue
    except TypeError:
        continue

##############
# Show image #
##############

fig,ax = plt.subplots()
cm = plt.get_cmap('cividis').copy()
if np.amax(img) > 1.0:
    cm.set_under('black')
    cm.set_over('black')
    ax.imshow(img, cmap=cm, norm=colors.SymLogNorm(linthresh=0.01, linscale=0.5, vmin=1))
else:
    ax.imshow(img, cmap=cm)
for ap in isoaplist:
    ap.plot()
plt.show()
