# autoap2
A set of Python scripts to: draw apertures quickly (semi-automatically with user verification), do photometry on these, pull select data from NED, and perform other SED analyses.

# Overview
In the main folder, there is a collection of Python scripts meant to be run individually. The 'param_files/' folder has a collection of configuration files which can be altered by the user (should be pretty self-explanatory). The 'lib/' folder has a collection of Python scripts which are used by other functions and should not need to be altered by the user. Here I'll explain what each script in the main directory does:

## create_apertures.py
When run, this will prompt the user to select a data folder (with structure according to 'structure' section) and go through the process of drawing regions around the detected flux. A main aperture will be drawn automatically (though the user can edit it), then the user will be prompted to select source-free region to use as the background, as well as contaminated regions near the target. Assuming this was done with the correct folder structure, these will be automatically used by the 'do_ap_phot.py' routine.

## do_ap_phot.py
Prompts the user to select a data folder which presumably contains all the requisite DS9 region and fits files. The user will be prompted to select which targets to do photometry on, and resulting photometry data will be dumped to a JSON file in the user-specified folder. Note that this will **always overwrite all previously-gathered data** for this target (except for NED data)!

## manual_entry.py
After photometry has been done, the user may wish to manually change measurements by running this script, which will prompt for a target, and then to change the flux value and upper/lower uncertainties.

## add_image_flags.py
Prompts the user to select a data folder where a JSON file as will be searched for (as created by the 'do_ap_phot.py' routine). For measurements which already exist in this data file, flags can be added as specified in the routine. Among other things, this is most useful for specifying how photometry is to be done on specific bands, such as taking an upper limit, etc.

### Flag Options
The following options are available for flags on each photometric band:
- u

    Specifies that when photometry is run on this band, it should be calculated as an upper limit (for use when there is no apparent emission). This takes the RMS average inside the aperture and multiplies it by 4.5.
    
- a

    This is the same as 'u' but for use with ALMA bands.
    
- r

    Specifies that there is suspected to be some significant unresolved contamination in the image (as determined from other images). If this flag is used, there should also be an extra aperture file with the same name as the main one but with "Lower" appended (e.g., 'SPIRE250Lower.fits'). This aperture should be a lower-end estimate for where only the flux from the target itself is contained (i.e., excluding contamination), while the main aperture file should include everything. The photometry is then measured with the main aperture, but the lower limit is adjusted lower for this lower bound.
    
- n

    **Not for explicit use by the user**, but specifies internally that a data point was gathered from NED.

## sed_pop_ned.py
When run, this will gather photometric data from NED to populate the SED for the user-specified target. Points are selected based on some hard-coded logic . . . if different logic is needed, the code itself needs to be altered. Running this on a target will always overwrite previously-gathered NED data (but data that was gathered otherwise).

## stitch_images.py
This will prompt the user for two .fits files, as well as an output file. It attempts to stitch the two images (preserving flux) by matching the WCS coordinates. The output .fits file header will contain the correct scaling/WCS, and everything else is copied from the first image. See the [reproject](https://reproject.readthedocs.io/) documentation for more details about how the stitching works.

## plot_seds.py
Will prompt the user for a data folder and a target from such to plot the SED. Data taken from the JSON file found in the data folder (assuming it has been created with either 'sed_pop_ned.py' or 'do_ap_phot.py'). Plotting options are found and can be edited in the 'param_files/' folder.

# Data structure
To be recognizable, the 'data' should be organized like the following example (this is the folder which should be selected when prompted in the above routines). This structure should be set up manually, and the apertures will be automatically created and placed in the correct location.
- data/
  - NGC1380/
    - fits/
      - W1.fits
      - SPIRE250.fits
      - etc.
    - apertures/
  - NGC3557/
    - fits/
      - W3.fits
      - PACS100.fits
      - etc.
    - apertures/
  - etc.

The following band names are currently supported (with predefined frequencies in the 'param_files/' folder, except for ALMA, which reads the frequency from the .fits header):
- ALMAExtended
- ALMAExtendedLowerNat
- ALMAExtendedUpperNat
- ALMAExtendedLowerBri
- ALMAExtendedUpperBri
- ALMANuclearNat
- ALMANuclearBri
- IRAC1
- IRAC2
- W1
- W2
- W3
- W4
- PACS70
- PACS100
- PACS160
- SPIRE250
- SPIRE350
- SPIRE500
- SDSS_u
- SDSS_g
- SDSS_r
- SDSS_i
- SDSS_z
- GALEX_NUV
- GALEX_FUV
- DES_u
- DES_Y
- DES_g
- DES_r
- DES_i
- DES_z
- 2MASS_J
- 2MASS_H
- 2MASS_K

# Dependencies
The following Python packages are required for one or more of the scripts:
- [astropy](https://www.astropy.org/)
- [matplotlib](https://matplotlib.org/)
- [numpy](https://numpy.org/)
- [pandas](https://pypi.org/project/pandas/)
- [regions](https://astropy-regions.readthedocs.io/)
- [scikit-image](https://scikit-image.org/)
- [astroquery](https://github.com/astropy/astroquery/)
- [termcolor](https://pypi.org/project/termcolor/)
- [reproject](https://reproject.readthedocs.io/)
- [DS9](https://ds9.si.edu/)
