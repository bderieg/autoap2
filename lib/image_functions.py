import matplotlib.pyplot as plt
import get_ned_data as gnd
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy import units as u
from regions import Regions, PixCoord, CirclePixelRegion, EllipseSkyRegion, EllipsePixelRegion, RectanglePixelRegion
import numpy as np
from skimage.morphology import erosion, disk
from skimage import feature, measure
import matplotlib.colors as colors
import PySimpleGUI as sg
from tkinter import ttk
import tkinter as tk
import threading
from matplotlib.patches import Ellipse
import subprocess as sp
from multiprocessing import Process


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
                img = ext.data
    except FileNotFoundError:
        logging.warning('For '+target_name+', \''+fltr+'.fits\' not found, but there exists a corresponding aperture file')
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


def find_blobs(img, background_flux):

    ##########################
    # Perform blob detection #
    ##########################

    # Run Laplacian of Gaussian method for detection
    blob_log = feature.blob_log(
            img, 
            min_sigma=1.5, 
            max_sigma=0.5*len(img), 
            threshold=background_flux+0.01*(np.nanmax(img)-background_flux),
            overlap=0.95
            )

    # Populate list of blobs' coordinates and radii
    blobs = []
    for blob in blob_log:
        y, x, r = blob
        blobs.append((x,y,r))

    ###########################################
    # Use moments to get blob characteristics #
    ###########################################

    # Define some empty arrays to fill
    centroids_x = []
    centroids_y = []
    radii = []
    position_angles = []
    flattenings = []

    # Iterate through all blobs
    for x, y, r in blobs:
        if (r == 0) or (x == 0) or (y == 0):
            continue

        # Blow up to a reasonably sized radius
        r *= 4

        # Define region for use in moment calculation
        source_region = img[int(y-r):int(y+r), int(x-r):int(x+r)].copy()

        # Arbitrarily threshold image to make moment calculation more effective
        thresh_bin_mask = source_region < (0.05*(np.max(source_region)-np.min(source_region))-np.min(source_region))
        source_region[thresh_bin_mask] = 0

        # Calculate raw moment and centroid
        raw_moment = measure.moments(source_region)
        centroids_x.append(x - (len(source_region[0])/2 - raw_moment[1,0]/raw_moment[0,0]))
        centroids_y.append(y - (len(source_region)/2 - raw_moment[0,1]/raw_moment[0,0]))

        # Calculate central moment and covariance matrix
        central_moment = measure.moments_central(source_region)
        central_moment_norm = central_moment / raw_moment[0,0]
        cov_mat = np.asarray([[central_moment_norm[2,0],central_moment_norm[1,1]],[central_moment_norm[1,1],central_moment_norm[0,2]]])

        # Get eigenvalues of covariance matrix
        eigenvalues = np.linalg.eigvals(cov_mat)

        # Calculate position angle
        pa = 0.5 * np.arctan2(2 * central_moment_norm[1,1], central_moment_norm[2,0] - central_moment_norm[0,2])
        position_angles.append(pa)

        # Calculate flattening
        flat = np.sqrt(np.min(eigenvalues) / np.max(eigenvalues))
        flattenings.append(flat)

        # Add radius to list
        radii.append(r)

    ########################################################
    # Remove blob if it's the target (i.e., in the center) #
    ########################################################

    target_blob_mask = [
            PixCoord(len(img[0])/2, len(img)/2)
            not in 
            EllipsePixelRegion(center=PixCoord(x=cx,y=cy), width=2*r*flat, height=2*r, angle=pa*u.rad) 
            for cx, cy, r, pa, flat in zip(centroids_x, centroids_y, radii, position_angles, flattenings)
            ]
    centroids_x = np.compress(target_blob_mask, centroids_x)
    centroids_y = np.compress(target_blob_mask, centroids_y)
    radii = np.compress(target_blob_mask, radii)
    position_angles = np.compress(target_blob_mask, position_angles)
    flattenings = np.compress(target_blob_mask, flattenings)

    ###################
    # Store apertures #
    ###################

    sub_aps = []
    for cx, cy, r, pa, flat in zip(centroids_x, centroids_y, radii, position_angles, flattenings):
        sub_aps.append(EllipsePixelRegion(center=PixCoord(x=cx,y=cy), width=2*r, height=2*r/flat, angle=pa*u.rad, visual={'edgecolor':'red'}))

    ###############################
    # Confirm apertures with user #
    ###############################

    # Display aperture with matplotlib for confirmation
    fig, ax = plt.subplots()
    cm = plt.get_cmap('cividis').copy()
    cm.set_under('black')
    cm.set_over('black')
    ax.imshow(img, cmap=cm, norm=colors.SymLogNorm(linthresh=0.01, linscale=0.5, vmin=1))
    for ap in sub_aps:
        ap.plot(ax=ax, lw=2.0)
    plt.show(block=False)

    # Is the user satisfied?
    print(' ')
    curinp = input('Are you satisfied with this aperture? (y/n) : ')

    # If the user is not satisfied . . . 
    if curinp == "n" or curinp == "N":
        # Temporarily save apertures
        Regions(sub_aps).write('./subaps_temp_autoap2.reg', format='ds9', overwrite=True)
        input('temp')

        # Open DS9 for user to edit aperture
        def open_ds9(img_path):
            sp.run(["ds9", img_path, "-region", './subaps_temp_autoap2.reg'])
        ds9proc = Process(target=open_ds9, args=(img_filename,))
        ds9proc.start()

        # Prompt user to save the new aperture when ready
        print(' ')
        input('A DS9 window should open; edit the apertures as desired, then type anything here when done to save (do not close DS9 manually) : ')
        ds9proc.terminate()

        # Read in newly-saved aperture file as the new region
        sp.run(["xpaset","-p","ds9","region","system","wcs"])
        sp.run(["xpaset","-p","ds9","region","sky","icrs"])
        sp.run(["xpaset","-p","ds9","region","skyformat","degrees"])
        sp.run(["xpaset","-p","ds9","region","save","./subaps_temp_autoap2.reg"])
        sub_aps = (Regions.read('./subaps_temp_autoap2.reg', format='ds9'))

        # Delete temporary region file
        sp.run(["rm", "./subaps_temp_autoap2.reg"])

    # Close pyplot if still open
    plt.close()

    ####################
    # Return apertures #
    ####################

    return sub_aps


def fit_ellipse_with_coordinates(img_filename, icrs_coord, axis_ratio, pa, background_flux):

    #####################
    # Set some stuff up #
    #####################

    # Define some constants
    major_ax = 0.01
    rad_increment = 0.005
    threshold = 0.005

    # Open image as 2D numpy array
    try:
        # Open
        fitsfile = fits.open(img_filename)
        # Loop through extensions to get the one with image data
        img = fitsfile[0].data  # Default to first extension just in case
        wcs = WCS(fitsfile[0].header, naxis=[1,2])  # ^^
        for ext in fitsfile:
            if 'IMAGE' in ext.header.get('EXTNAME','').upper():
                img = ext.data
                wcs = WCS(ext.header, naxis=[1,2])
    except FileNotFoundError:
        logging.warning('For '+target_name+', \''+fltr+'.fits\' not found, but there exists a corresponding aperture file')
        return False, False

    ##############################
    # Find correct aperture size #
    ##############################

    # Make aperture bigger, iterating until a threshold is reached
    flux_max = 1e9  # Arbitrary starting value
    flux_edge_avg = flux_max
    while (flux_edge_avg-background_flux) >= (threshold*(flux_max-background_flux)):
        # Increment aperture size
        major_ax += rad_increment

        # Define the main aperture
        main_ap_sky = EllipseSkyRegion(
                center=SkyCoord(icrs_coord[0], icrs_coord[1], unit='deg', frame='icrs'), 
                height=major_ax*axis_ratio*u.deg, 
                width=major_ax*u.deg, 
                angle=pa*u.deg,
                visual={'edgecolor':'green'}
                )
        main_ap_props = {
                    "x" : (main_ap_sky.to_pixel(wcs)).center.xy[0],
                    "y" : (main_ap_sky.to_pixel(wcs)).center.xy[1],
                    "width" : (main_ap_sky.to_pixel(wcs)).width,
                    "height" : (main_ap_sky.to_pixel(wcs)).height,
                    "angle" : (main_ap_sky.to_pixel(wcs)).angle.value
                }

        # Get cutouts
        aperture_full_mask = (main_ap_sky.to_pixel(wcs)).to_mask()
        aperture_edge_mask = aperture_full_mask - erosion(aperture_full_mask, disk(2))
        aperture_full_cutout = aperture_full_mask.multiply(img)
        aperture_edge_cutout = aperture_full_cutout * aperture_edge_mask 

        flux_max = np.amax(aperture_full_cutout)

        # Find the average flux around the outside edge of the aperture
        flux_edge_avg = np.median(aperture_edge_cutout[np.nonzero(aperture_edge_cutout)])

    ##############################
    # Confirm aperture with user #
    ##############################

    # Display aperture with matplotlib for confirmation
    fig, ax = plt.subplots()
    cm = plt.get_cmap('cividis').copy()
    cm.set_under('black')
    cm.set_over('black')
    ax.imshow(img, cmap=cm, norm=colors.SymLogNorm(linthresh=0.01, linscale=0.5, vmin=1))
    (main_ap_sky.to_pixel(wcs)).plot(ax=ax, lw=2.0)
    plt.show(block=False)

    # Is the user satisfied?
    print(' ')
    curinp = input('Are you satisfied with this aperture? (y/n) : ')

    # If the user is not satisfied . . . 
    if curinp == "n" or curinp == "N":
        # Temporarily save aperture
        main_ap_sky.write('./mainap_temp_autoap2.reg', format='ds9', overwrite=True)

        # Open DS9 for user to edit aperture
        def open_ds9(img_path):
            sp.run(["ds9", img_path, "-region", './mainap_temp_autoap2.reg'])
        ds9proc = Process(target=open_ds9, args=(img_filename,))
        ds9proc.start()

        # Prompt user to save the new aperture when ready
        print(' ')
        input('A DS9 window should open; edit the aperture as desired, then type anything here when done to save (do not close DS9 manually) : ')
        ds9proc.terminate()

        # Read in newly-saved aperture file as the new region
        sp.run(["xpaset","-p","ds9","region","system","wcs"])
        sp.run(["xpaset","-p","ds9","region","sky","icrs"])
        sp.run(["xpaset","-p","ds9","region","skyformat","degrees"])
        sp.run(["xpaset","-p","ds9","region","save","./mainap_temp_autoap2.reg"])
        main_ap_sky = (Regions.read('./mainap_temp_autoap2.reg', format='ds9'))[0]

        # Delete temporary region file
        sp.run(["rm", "./mainap_temp_autoap2.reg"])

    # Close pyplot if still open
    plt.close()

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

    fig, ax = plt.subplots()
    cm = plt.get_cmap('cividis').copy()
    cm.set_under('black')
    cm.set_over('black')
    ax.imshow(img, cmap=cm, norm=colors.SymLogNorm(linthresh=0.01, linscale=0.5, vmin=1))
    (main_ap_sky.to_pixel(wcs)).plot(ax=ax, lw=2.0)
    cid = fig.canvas.mpl_connect('button_press_event', onclick)
    plt.show()

    bg_ap_pix = RectanglePixelRegion(PixCoord(bg_x,bg_y), width=10, height=10, visual={'edgecolor':'blue'})

    return main_ap_sky, bg_ap_pix
