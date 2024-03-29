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
pa_data_loc = workdir + 'isophote_pa_data.json'

# Prompt to specify target
print(' ')
target = input("Enter the desired target : ")

###############
# Import data #
###############

fltr = "W1"

img_path = workdir+target+"/fits/"+fltr+".fits"
ap_path = workdir+target+"/apertures/"+fltr+".reg"
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

################
# Fit isophote #
################

# Create fitting geometry
geom = EllipseGeometry(
            x0=hlap.center.xy[0],
            y0=hlap.center.xy[1],
            sma=1.5*hlap.width,
            eps=(1-hlap.height/hlap.width),
            pa=hlap.angle.value * np.pi/180
        )

# Fit
fitobj = Ellipse(img, geom)
isofitlist = fitobj.fit_image(
                    sma0=hlrad, 
                    minsma=0.1*hlrad, 
                    maxsma=2.0*hlrad, 
                    step=0.2*hlrad, 
                    linear=True, 
                    fix_center=True,
                    conver=1e-3, 
                    minit=20,
                    nclip=12, 
                    sclip=2.0,
                    maxgerr=2.0, 
                    fflag=0.1 
                )

# Save as apertures
isoaplist = []
pa_list = []
pa_unc_list = []
rad_list = []
flat_list = []
factor = True
for aa in isofitlist:
    try:
        if factor:
            pa_list.append(aa.pa*180/np.pi)
            pa_unc_list.append(aa.pa_err*180/np.pi)
            rad_list.append(aa.sma)
            flat_list.append(1.-aa.eps)
        apcolor = 'white'
        if aa.sma/hlrad > 0.9 and aa.sma/hlrad < 1.1:
            apcolor = 'red'
            factor = False
        isoaplist.append(EllipsePixelRegion(
                        PixCoord(aa.x0, aa.y0),
                        aa.sma,
                        aa.sma * (1-aa.eps),
                        angle=aa.pa * 180/np.pi * u.deg,
                        visual={'edgecolor':apcolor}
                    ))
    except ValueError:
        continue

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

# Get PA from ellipse fit
pa_list = [pa for pa,pau in zip(pa_list,pa_unc_list) if pau!=0.0]
pa_unc_list = [pau for pau in pa_unc_list if pau!=0.0]
avg_pa_sig = sum([(1/punc**2) * p for p,punc in zip(pa_list, pa_unc_list)]) / sum([(1/punc**2) for punc in pa_unc_list])
avg_pa_sig_unc = np.sqrt(
                        sum(
                            [
                                (1/pauu**2)*ss**2 for pauu,ss in\
                                zip(
                                        pa_unc_list,\
                                        [np.sqrt((1/pau**2)/sum([(1/punc**2) for punc in pa_unc_list])) for pau in pa_unc_list]
                                    )
                            ]
                            ) 
                            / (len(pa_list)*(len(pa_list)-1))
                        )
avg_pa_sig += 90
if avg_pa_sig > 180.0 : avg_pa_sig -= 180.0
pa_data[target] = {}
pa_data[target]["pa"] = avg_pa_sig
pa_data[target]["pa_unc"] = avg_pa_sig_unc

# Write
pa_outfile = open(pa_data_loc, 'w')
json.dump(pa_data, pa_outfile, indent=5)

##############
# Show image #
##############

plt.figure(0)
fig,ax = plt.subplots()
cm = plt.get_cmap('cividis').copy()
if np.amax(img) > 0.1:
    cm.set_under('black')
    cm.set_over('black')
    ax.imshow(img, cmap=cm, norm=colors.SymLogNorm(linthresh=0.01, linscale=0.5, vmin=0.1))
else:
    ax.imshow(img, cmap=cm)
main_ap_pix.plot()
for aa in [jeff for jeff in isoaplist if jeff.visual['edgecolor']=='red']:
    aa.plot()
for aa in [jeff for jeff in isoaplist if jeff.visual['edgecolor']=='white']:
    aa.plot()
plt.legend(['Integration Region','Isophote Limit','Other Isophotes'])
ax.invert_yaxis()

plt.figure(1)
fig,ax = plt.subplots()
cm = plt.get_cmap('cividis').copy()
if np.amax(img) > 0.1:
    cm.set_under('black')
    cm.set_over('black')
    ax.imshow(img, cmap=cm, norm=colors.SymLogNorm(linthresh=0.01, linscale=0.5, vmin=0.1))
else:
    ax.imshow(img, cmap=cm)
ax.invert_yaxis()

plt.show()
