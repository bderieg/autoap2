import matplotlib.pyplot as plt
import get_ned_data as gnd
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy import units as u
from regions import Regions, PixCoord, EllipseSkyRegion, EllipsePixelRegion, RectanglePixelRegion
import numpy as np
from skimage.morphology import erosion, disk
from skimage import feature, measure
import matplotlib.colors as colors
from matplotlib.patches import Ellipse
import subprocess as sp
from multiprocessing import Process

import logging


def integrate_flux(cutout, background):
    flux = 0

    for row in cutout:
        for pix_val in row:
            if pix_val==pix_val and pix_val!=0.0:
                flux += pix_val
                flux -= background

    return flux


def find_background_flux(img_filename):

    ##############
    # Open image #
    ##############

    # Attempt to open fits file
    try:
        # Open
        fitsfile = fits.open(img_filename)
        # Loop through extensions to get the one with image data
        img = fitsfile[0].data  # Default to first extension just in case
        for ext in fitsfile:
            if 'IMAGE' in ext.header.get('EXTNAME','').upper():
                # Get image data
                img = ext.data
        ## Get rid of redundant dimensions in data if necessary
        img = img.squeeze()
    except FileNotFoundError:
        return 0.0

    #############################
    # Get some image parameters #
    #############################

    flux_sigma = np.nanstd(img.flatten())
    flux_median = np.nanmedian(img.flatten())

    ##############
    # Clip image #
    ##############

    img_sigma_clipped = img
    img_sigma_clipped[np.where(img > flux_median)] = 0

    ##########################################
    # Get background flux from clipped image #
    ##########################################
    
    background_flux = np.nanmedian((img_sigma_clipped.flatten())[img_sigma_clipped.flatten() != 0])
    return background_flux


def find_blobs(full_img, main_mask, background_flux):

    ##############
    # Get cutout #
    ##############

    # Cutout
    img = main_mask.multiply(full_img)
    # Indices of removed pixels to be readded later
    xcut = ((main_mask.get_overlap_slices(full_img.shape)[0])[1]).start
    ycut = ((main_mask.get_overlap_slices(full_img.shape)[0])[0]).start

    ######################
    # Set blob threshold #
    ######################

    bth = 0.1 * (np.nanmax(img) - background_flux)

    ###########################
    # Create dummy WCS object #
    ###########################

    header = fits.Header()
    header['SIMPLE'] = 'T'
    header['BITPIX'] = -32
    header['NAXIS'] = 2
    header['NAXIS1'] = len(img)
    header['NAXIS2'] = len(img[0])
    header['CTYPE1'] = 'RA---TAN'
    header['CTYPE2'] = 'DEC--TAN'
    header['CRPIX1'] = 0.0
    header['CRPIX2'] = 0.0
    header['CRVAL1'] = 0.0
    header['CRVAL2'] = 0.0
    header['CDELT1'] = 1.0
    header['CDELT2'] = 1.0
    header['CUNIT1'] = 'deg'
    header['CUNIT2'] = 'deg'
    header['OBJECT'] = 'none'
    header['EXPTIME'] = 1.0
    wcs = WCS(header)

    #############################################
    # Prompt user to input background apertures #
    #############################################

    # Temporarily save image as fits
    fits.writeto('./_temp_main_ap_autoap2.fits', img, header=header, overwrite=True)

    # Open DS9 for user to edit aperture
    def open_ds9():
        sp.run(["ds9", './_temp_main_ap_autoap2.fits'])
    ds9proc = Process(target=open_ds9)
    ds9proc.start()

    # Prompt user to save the new aperture when ready
    print(' ')
    input('A DS9 window should open; edit the apertures as desired, then type anything here when done to save (do not close DS9 manually) : ')

    # Read in newly-saved aperture file as the new region
    sp.run(["xpaset","-p","ds9","region","system","wcs"])
    sp.run(["xpaset","-p","ds9","region","sky","icrs"])
    sp.run(["xpaset","-p","ds9","region","skyformat","degrees"])
    sp.run(["xpaset","-p","ds9","region","save","./_temp_subaps_autoap2.reg"])
    sub_aps = Regions.read('./_temp_subaps_autoap2.reg', format='ds9')
    sub_aps = [ap.to_pixel(wcs) for ap in sub_aps]

    # Close DS9
    sp.run(["xpaset","-p","ds9","exit"])

    # Delete temporary region file
    sp.run(["rm", "./_temp_subaps_autoap2.reg"])
    sp.run(["rm", "./_temp_main_ap_autoap2.fits"])

    ############################################
    # Re-add full image coordinates and return #
    ############################################

    for ap in sub_aps:
        oldx = ap.center.x
        oldy = ap.center.y
        ap.center = PixCoord(oldx+xcut, oldy+ycut)

    return sub_aps


def fit_ellipse_with_coordinates(img_filename, icrs_coord, axis_ratio, background_flux):

    #####################
    # Set some stuff up #
    #####################

    # Define some constants
    major_ax = 1e-3
    rad_increment = 5e-3
    threshold = 0.0005
    pa = 0  # Arbitrary position angle; will be fixed later
    sb_res = 30

    # Open image as 2D numpy array
    try:
        # Open
        fitsfile = fits.open(img_filename)
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
        return False, False

    ################################
    # Prompt user to draw aperture #
    ################################

    # Open DS9 for user to edit aperture
    def open_ds9(img_path):
        sp.run(["ds9", img_path])
    ds9proc = Process(target=open_ds9, args=(img_filename,))
    ds9proc.start()

    # Prompt user to save the new aperture when ready
    print(' ')
    input('A DS9 window should open; edit the aperture as desired, then type anything here when done to save (do not close DS9 manually) : ')

    # Read in newly-saved aperture file as the new region
    sp.run(["xpaset","-p","ds9","region","system","wcs"])
    sp.run(["xpaset","-p","ds9","region","sky","icrs"])
    sp.run(["xpaset","-p","ds9","region","skyformat","degrees"])
    sp.run(["xpaset","-p","ds9","region","save","./_temp_mainap_autoap2.reg"])
    main_ap_sky = (Regions.read('./_temp_mainap_autoap2.reg', format='ds9'))[0]

    # Close DS9
    sp.run(["xpaset","-p","ds9","exit"])

    # Delete temporary region file
    sp.run(["rm", "./_temp_mainap_autoap2.reg"])

    ####################################
    # Ask user for background location #
    ####################################

    print(' ')
    print('Click to select the center of a 10x10 background region; close window after final selection')

    # Define click event
    global bg_x, bg_y
    bg_x = 0.0
    bg_y = 0.0
    def onclick(event):
        global bg_x, bg_y
        bg_x = event.xdata
        bg_y = event.ydata
        print(' ')
        print('Last selection: (', event.xdata, ',', event.ydata, ')')

    # Display image
    fig, ax = plt.subplots()
    cm = plt.get_cmap('cividis').copy()
    if np.amax(img) > 1.0:
        cm.set_under('black')
        cm.set_over('black')
        ax.imshow(img, cmap=cm, norm=colors.SymLogNorm(linthresh=0.01, linscale=0.5, vmin=1))
    else:
        ax.imshow(img, cmap=cm)
    (main_ap_sky.to_pixel(wcs)).plot(ax=ax, lw=2.0)
    cid = fig.canvas.mpl_connect('button_press_event', onclick)
    plt.show()

    # Define background aperture based on final click location
    bg_ap_pix = RectanglePixelRegion(PixCoord(bg_x,bg_y), width=10, height=10, visual={'edgecolor':'blue'})

    return main_ap_sky, bg_ap_pix


def calc_unc_background(img_data, bg_ap):

    bg_cutout = (bg_ap.to_mask()).multiply(img_data)
    bg_cutout = bg_cutout[1:-1,1:-1]  # Trim 0s

    return np.sqrt( np.mean( bg_cutout**2 ) ) * len(bg_cutout)**2

def calc_unc_apcopy(img_data, main_ap, bg=0.0):
    # Set constants
    bg_peak_tol = 0.15
    bg_avg_tol = 0.005
    iters = 0

    # Mask image
    masked_img = (-1*main_ap.to_mask().to_image(np.shape(img_data)) + 1) * img_data
    masked_img = np.where(masked_img == 0.0, np.nan, masked_img)

    # Copy apertures and accept if they're good
    fluxes = []
    for i in range(20):
        tempcutout = [[0.0]]
        while np.isnan(tempcutout).any()\
                or (np.nanmax(tempcutout)-bg)>bg_peak_tol*(np.nanmax(img_data)-bg)\
                or (np.average(tempcutout)-bg)>bg_avg_tol*(np.nanmax(img_data)-bg):
            while tempcutout is None:
                main_ap.center = PixCoord(np.random.randint(0, len(img_data)), np.random.randint(0, len(img_data)))
                tempcutout = main_ap.to_mask().multiply(img_data)
            iters += 1
            if iters > 20:
                bg_avg_tol += 0.001
        fluxes.append(integrate_flux(tempcutout, bg))

    return np.std(fluxes)
