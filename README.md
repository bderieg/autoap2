(in progress)

# autoap2
A set of Python scripts to: draw apertures quickly (semi-automatically with user verification), do photometry on these, pull select data from NED, and perform other SED analyses.

# Overview
In the main folder, there is a collection of Python scripts meant to be run individually. The 'param_files/' folder has a collection of configuration files which can be altered by the user (should be pretty self-explanatory). The 'lib/' folder has a collection of Python scripts which are used by other functions and should not need to be altered by the user. Here I'll explain what each script in the main directory does:

## create_apertures.py
When run, this will prompt the user to select a data folder (with structure according to 'structure' section) and go through the process of drawing regions around the detected flux. A main aperture will be drawn automatically (though the user can edit it), then the user will be prompted to select source-free region to use as the background, as well as contaminated regions near the target. Assuming this was done with the correct folder structure, these will be automatically used by the 'do_ap_phot.py' routine.

## do_ap_phot.py

## add_image_flags.py

## sed_pop_ned.py

## stitch_images.py

## plot_seds.py

# Data structure
blah

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
