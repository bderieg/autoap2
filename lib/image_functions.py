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
import threading
from matplotlib.patches import Ellipse
import subprocess as sp
from multiprocessing import Process


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

    ##########################
    # Perform blob detection #
    ##########################

    # Run Laplacian of Gaussian method for detection
    blob_log = feature.blob_log(
            img, 
            min_sigma=0.5, 
            max_sigma=0.1*len(img), 
            num_sigma=40,
            threshold_rel=0.005,
            overlap=0.0
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

        # Recover radius (blow up to a reasonable size)
        r *= 2

        # Define region for use in moment calculation
        source_region = img[int(y-r):int(y+r), int(x-r):int(x+r)].copy()
        
        # Filter some apertures
        ## Blank aperture
        if source_region.size == 0:
            continue
        ## Not masking anything significant
        if np.nanmax(source_region) < bth:
            continue

        # Arbitrarily threshold image to make moment calculation more effective
        thresh_bin_mask = source_region < (0.8*(np.max(source_region)-np.min(source_region))-np.min(source_region))
        source_region[thresh_bin_mask] = 0

        # Calculate raw moment and centroid
        raw_moment = measure.moments(source_region)

        # Calculate central moment and covariance matrix
        central_moment = measure.moments_central(source_region)
        central_moment_norm = central_moment / raw_moment[0,0]
        cov_mat = np.asarray([[central_moment_norm[2,0],central_moment_norm[1,1]],[central_moment_norm[1,1],central_moment_norm[0,2]]])

        # Get eigenvalues of covariance matrix
        eigenvalues = np.linalg.eigvals(cov_mat)

        # Calculate position angle
        pa = 0.5 * np.arctan2(2 * central_moment_norm[1,1], central_moment_norm[2,0] - central_moment_norm[0,2])

        # Calculate flattening
        flat = np.sqrt(np.min(eigenvalues) / np.max(eigenvalues))

        # Throw away aperture in case of weirdness
        if flat == 0 or\
                flat == 1 or\
                flat < 0.05 or\
                r <= 1.0 or\
                r >= 0.5*len(img):
            continue

        # Add values to lists
        centroids_x.append(x - (len(source_region[0])/2 - raw_moment[1,0]/raw_moment[0,0]))
        centroids_y.append(y - (len(source_region)/2 - raw_moment[0,1]/raw_moment[0,0]))
        flattenings.append(flat)
        position_angles.append(pa)
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

    #############################
    # Store apertures as a list #
    #############################

    sub_aps = []
    for cx, cy, r, pa, flat in zip(centroids_x, centroids_y, radii, position_angles, flattenings):
        sub_aps.append(EllipsePixelRegion(center=PixCoord(x=cx,y=cy), width=2*r, height=2*r/flat, angle=pa*u.rad, visual={'edgecolor':'red'}))

    ###############################
    # Confirm apertures with user #
    ###############################

    # Display aperture with matplotlib
    fig, ax = plt.subplots()
    cm = plt.get_cmap('cividis').copy()
    cm.set_under('black')
    cm.set_over('black')
    ax.imshow(img, cmap=cm, norm=colors.SymLogNorm(linthresh=0.01, linscale=0.5, vmin=1))
    ## Plot subtraction apertures over image
    for ap in sub_aps:
        ap.plot(ax=ax, lw=2.0)
    plt.show(block=False)  # To keep window open with terminal access

    # Is the user satisfied?
    print(' ')
    curinp = input('Are you satisfied with these subtraction apertures? (y/n) : ')

    # If the user is not satisfied . . . 
    if curinp == "n" or curinp == "N":
        # Temporarily save apertures
        Regions(sub_aps).write('./_temp_subaps_autoap2.reg', format='ds9', overwrite=True)

        # Temporarily save image as fits
        fits.writeto('./_temp_main_ap_autoap2.fits', img, header=header, overwrite=True)

        # Open DS9 for user to edit aperture
        def open_ds9():
            sp.run(["ds9", './_temp_main_ap_autoap2.fits', "-region", './_temp_subaps_autoap2.reg'])
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

    # Close pyplot if still open
    plt.close()

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
    except FileNotFoundError:
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

    #######################################
    # Rotate ellipse to align with galaxy #
    #######################################

    # Shrink the ellipse a little for better fitting
    major_ax /= 4

    # Calculate edge flux variance for each pa
    pas = []
    variances = []
    for curpa in range(0, 360, 2):
        # Define the main aperture
        main_ap_sky = EllipseSkyRegion(
                center=SkyCoord(icrs_coord[0], icrs_coord[1], unit='deg', frame='icrs'), 
                height=major_ax*axis_ratio*u.deg, 
                width=major_ax*u.deg, 
                angle=curpa*u.deg,
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

        # Calculate variance
        flux_edge_var = np.var(aperture_edge_cutout[np.nonzero(aperture_edge_cutout)])
        pas.append(curpa)
        variances.append(flux_edge_var)

    # Reset aperture size
    major_ax *= 4

    # Get minimum variance and set pa
    minvarind = np.where(variances == np.min(variances))[0][0]
    pa = pas[minvarind]
    main_ap_sky = EllipseSkyRegion(
            center=SkyCoord(icrs_coord[0], icrs_coord[1], unit='deg', frame='icrs'), 
            height=major_ax*axis_ratio*u.deg, 
            width=major_ax*u.deg, 
            angle=pa*u.deg,
            visual={'edgecolor':'green'}
            )

    ########################################
    # Calculate surface brightness profile #
    ########################################

    # Define some lists
    rads = []
    sbs = []

    # Continually shrink the ellipse and calculate surface brightness
    for temp_major_ax in np.logspace(np.log10(major_ax), np.log10(major_ax/sb_res), sb_res, base=10):
        # Define the main aperture
        main_ap_sky = EllipseSkyRegion(
                center=SkyCoord(icrs_coord[0], icrs_coord[1], unit='deg', frame='icrs'), 
                height=temp_major_ax*axis_ratio*u.deg, 
                width=temp_major_ax*u.deg, 
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

        # Get edge surface brightness
        flux_edge_avg = np.median(aperture_edge_cutout[np.nonzero(aperture_edge_cutout)])

        # Update lists
        rads.append(temp_major_ax)
        sbs.append(flux_edge_avg)

    #########################
    # Restore aperture size #
    #########################

    main_ap_sky = EllipseSkyRegion(
            center=SkyCoord(icrs_coord[0], icrs_coord[1], unit='deg', frame='icrs'), 
            height=major_ax*axis_ratio*u.deg, 
            width=major_ax*u.deg, 
            angle=pa*u.deg,
            visual={'edgecolor':'green'}
            )

    ###################################
    # Plot surface brightness profile #
    ###################################

    plt.loglog([i*3600 for i in rads], sbs, 'k.')
    plt.title('Surface Brightness Profile')
    plt.xlabel('Major Axis (arcsec)')
    plt.ylabel('Average Surface Brightness (image unit)')
    plt.show(block=False)  # Open matplotlib but keep terminal access

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
    plt.show(block=False)  # Open matplotlib but keep terminal access

    # Is the user satisfied?
    print(' ')
    curinp = input('Are you satisfied with this aperture? (y/n) : ')

    # If the user is not satisfied . . . 
    if curinp == "n" or curinp == "N":
        # Temporarily save aperture
        main_ap_sky.write('./_temp_mainap_autoap2.reg', format='ds9', overwrite=True)

        # Open DS9 for user to edit aperture
        def open_ds9(img_path):
            sp.run(["ds9", img_path, "-region", './_temp_mainap_autoap2.reg'])
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

    # Close pyplot if still open
    plt.close('all')

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
    cm.set_under('black')
    cm.set_over('black')
    ax.imshow(img, cmap=cm, norm=colors.SymLogNorm(linthresh=0.01, linscale=0.5, vmin=1))
    (main_ap_sky.to_pixel(wcs)).plot(ax=ax, lw=2.0)
    cid = fig.canvas.mpl_connect('button_press_event', onclick)
    plt.show()

    # Define background aperture based on final click location
    bg_ap_pix = RectanglePixelRegion(PixCoord(bg_x,bg_y), width=10, height=10, visual={'edgecolor':'blue'})

    return main_ap_sky, bg_ap_pix
