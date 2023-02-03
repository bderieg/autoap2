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

def find_background_coord(img_filename):
    # Open image as 2D numpy array
    img = fits.open(img_filename)[0].data

    # Perform blob detection
    blob_log = feature.blob_log(img, min_sigma=80, max_sigma=90, threshold=10.0)
    blobs = []
    for blob in blob_log:
        y, x, r = blob
        blobs.append((x,y,r))

    # Compute blob image moments
    moments = []
    for x, y, r in blobs:
        region = img[int(y-r):int(y+r), int(x-r):int(x+r)]
        moment = measure.moments_central(region, (y,x))
        moments.append(moment)

    # Use moments to get position angles and flattenings
    position_angles = []
    flattenings = []
    for moment in moments:
        # Get eigenvalues of covariance matrix from blob moment
        eigenvalues, _ = np.linalg.eig([[moment[2,0],moment[1,1]],[moment[1,1],moment[0,2]]])

        pa = 0.5 * np.arctan2(2 * moment[1,1], moment[2,0] - moment[0,2])
        position_angles.append(pa)

        flat = 1 - np.sqrt(np.min(eigenvalues) / np.max(eigenvalues))
        flattenings.append(flat)

    # Display detected blobs
    fig, ax = plt.subplots()
    ax.imshow(img, cmap='gray')
    for blob, pa, flat in zip(blobs, position_angles, flattenings):
        region = EllipsePixelRegion(center=PixCoord(x=blob[0],y=blob[1]), width=2*blob[2]*flat, height=2*blob[2], angle=pa*u.rad)
        region.plot(ax=ax, color='red', lw=2.0)
    plt.show()


def fit_ellipse_with_coordinates(img_filename, icrs_coord, axis_ratio, pa, background_flux):
    #####################
    # Set some stuff up #
    #####################

    # Define some constants
    major_ax = 0.01
    rad_increment = 0.005
    threshold = 0.00001

    # Open image as 2D numpy array
    img = fits.open(img_filename)[0].data

    # Define WCS system to convert to pixel apertures
    wcs = WCS(fits.open(img_filename)[0].header, naxis=[1,2])

    ##############################
    # Find correct aperture size #
    ##############################

    # Make aperture bigger, iterating until a threshold is reached
    flux_max = 10  # Arbitrary starting value
    flux_edge_avg = flux_max
    while (flux_edge_avg-background_flux) > (threshold*(flux_max-background_flux)):
        # Increment aperture size
        major_ax += rad_increment

        # Define the main aperture
        main_ap_sky = EllipseSkyRegion(
                center=SkyCoord(icrs_coord[0], icrs_coord[1], unit='deg', frame='icrs'), 
                height=major_ax*u.deg, 
                width=major_ax*axis_ratio*u.deg, 
                angle=pa*u.deg
                )
        main_ap_pix = main_ap_sky.to_pixel(wcs)

        # Get cutouts
        aperture_full_mask = main_ap_pix.to_mask()
        aperture_edge_mask = aperture_full_mask - erosion(aperture_full_mask, disk(2))
        aperture_full_cutout = aperture_full_mask.multiply(img)
        aperture_edge_cutout = aperture_full_cutout * aperture_edge_mask 

        # Find max flux in aperture
        flux_max = np.amax(aperture_full_cutout)

        # Find the average flux around the outside edge of the aperture
        flux_edge_avg = np.median(aperture_edge_cutout[np.nonzero(aperture_edge_cutout)])

    plt.imshow(img)
    main_ap_pix.plot()
    plt.show()

# fit_with_coordinates(
#         "/home/ben/Desktop/research/research_boizelle_working/FITS/NGC6861DES_i.fits",
#         [301.8310264, -48.3699881],
#         gnd.get_ellipse_parameters("NGC 6861")["axis_ratio"],
#         gnd.get_ellipse_parameters("NGC 6861")["position_angle"]
#         )

b = find_background_coord("/home/ben/Desktop/research/research_boizelle_working/FITS/NGC1380DES_r.fits")
