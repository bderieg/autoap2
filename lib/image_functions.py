import matplotlib.pyplot as plt
import get_ned_data as gnd
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy import units as u
from regions import PixCoord, CirclePixelRegion, EllipseSkyRegion, EllipsePixelRegion
import numpy as np
from skimage.morphology import erosion, disk
from skimage import feature, measure
import matplotlib.colors as colors
import PySimpleGUI as sg
from tkinter import ttk
import tkinter as tk
import threading
from matplotlib.patches import Ellipse



def find_background_flux(img_filename):

    ##############
    # Open image #
    ##############

    img = fits.open(img_filename)[0].data

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

    #####################################
    # Confirm blob detections with user #
    #####################################

    fig, ax = plt.subplots()
    cm = plt.get_cmap('cividis').copy()
    cm.set_under('black')
    cm.set_over('black')
    ax.imshow(img, cmap=cm, norm=colors.SymLogNorm(linthresh=0.01, linscale=0.5, vmin=1))
    for cx, cy, r, pa, flat in zip(centroids_x, centroids_y, radii, position_angles, flattenings):
        region = EllipsePixelRegion(center=PixCoord(x=cx,y=cy), width=2*r, height=2*r/flat, angle=pa*u.rad)
        region.plot(ax=ax, color='red', lw=2.0)
    plt.show()


def fit_ellipse_with_coordinates(img_filename, icrs_coord, axis_ratio, pa, background_flux):

    #####################
    # Set some stuff up #
    #####################

    # Define some constants
    major_ax = 0.01
    rad_increment = 0.005
    threshold = 0.005

    # Open image as 2D numpy array
    img = fits.open(img_filename)[0].data

    # Define WCS system to convert to pixel apertures
    wcs = WCS(fits.open(img_filename)[0].header, naxis=[1,2])

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
                angle=pa*u.deg
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

    class EditMainRegion:
        def __init__(self, img, ap_props):
            self.fig, self.ax = plt.subplots()
            self.cm = plt.get_cmap('cividis').copy()
            self.cm.set_under('black')
            self.cm.set_over('black')
            self.ax.imshow(img, cmap=self.cm, norm=colors.SymLogNorm(linthresh=0.01, linscale=0.5, vmin=1))
            self.region = Ellipse(xy=(main_ap_props['x'],main_ap_props['y']), width=main_ap_props['width'], height=main_ap_props['height'], angle=main_ap_props['angle'], fill=False, ec='red')
            self.ax.add_artist(self.region)

    app = EditMainRegion(img, main_ap_props)

    return aperture_full_cutout

name = "NGC 1380"

target_cutout = fit_ellipse_with_coordinates(
        "/home/ben/Desktop/research/research_boizelle_working/FITS/NGC13802MASS_H.fits",
        gnd.get_coords(name),
        gnd.get_ellipse_parameters(name)["axis_ratio"],
        gnd.get_ellipse_parameters(name)["position_angle"],
        find_background_flux("/home/ben/Desktop/research/research_boizelle_working/FITS/NGC13802MASS_H.fits")
        )

find_blobs(target_cutout, find_background_flux("/home/ben/Desktop/research/research_boizelle_working/FITS/NGC13802MASS_H.fits"))
