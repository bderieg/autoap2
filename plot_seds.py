import os

import pandas as pd
import numpy as np
import json
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter, LogLocator
from decimal import Decimal

#################
# Set constants #
#################

workdir = ""

#######################
#######################
## Start main script ##
#######################
#######################

#################################
# Set constants with user input #
#################################

# Prompt to change working directory
print(' ')
workdir = input("Enter working directory (relative or absolute) : ")
if workdir[-1] != '/':
    workdir += '/'

# Prompt to specify targets
print(' ')
curinp = input("Plot SEDs for ALL subfolders in this directory? (y/n) : ")
subdirs = [f.name for f in os.scandir(workdir) if f.is_dir()]
if curinp=='y' or curinp=='Y':
    target_names = subdirs
else:
    print(' ')
    target_names = input("Enter (as a comma-separated list) the desired targets : ")
    target_names = target_names.split(',')

######################
# Open SED data file #
######################

# Set SED data file path
if workdir[-1] != '/':
    workdir += '/'
sed_data_loc = workdir + "sed_data.json"

# Import file
try:
    seds = pd.read_json(sed_data_loc)
except FileNotFoundError:
    print(' ')
    print('SED data file not found . . . exiting')
    exit()

########################
# Loop through targets #
########################

sed_ps = json.load(open('./param_files/telescope_sed_pointstyles.json'))
sed_colors = json.load(open('./param_files/telescope_sed_colors.json'))

for target in target_names:
    if target in subdirs:
        # Import SED data
        sed_data = seds[target]['sed_data']
        sed_telenames = seds[target]['sed_telescopenames']
        sed_filternames = seds[target]['sed_filternames']
        sed_data_arr = np.transpose(np.asarray([[float(key), sed_data[key]] for key in sed_data]))
        sed_telenames_arr = np.transpose(np.asarray([sed_telenames[key] for key in sed_telenames]))
        sed_filternames_arr = np.transpose(np.asarray([sed_filternames[key] for key in sed_filternames]))
        # Make unique legend list
        unique_legend_points = ["_"] * len(sed_telenames_arr)
        for ind in np.unique(np.asarray(sed_telenames_arr), return_index=True)[1]:
            unique_legend_points[ind] = sed_telenames_arr[ind] 
        # Plot
        fig, ax = plt.subplots()
        for xval,yval,filtername,telename in zip(sed_data_arr[0], sed_data_arr[1], sed_filternames_arr, unique_legend_points):
            ax.scatter(xval, yval, marker=sed_ps[filtername], color=sed_colors[filtername], label=telename)
        # Set other plot parameters
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.legend()
        ax.set_title(target)
        ax.set_xlabel('log$_{10}$ Rest Frequency (Hz)')
        ax.set_ylabel('log$_{10}$ Flux Density (Jy)')
        ax.xaxis.set_major_formatter(FuncFormatter(lambda x,pos : np.log10(x)))
        ax.yaxis.set_major_formatter(FuncFormatter(lambda x,pos : np.log10(x)))
        # Set twin wavelength axis
        axwl = ax.twiny()
        axwl.set_xlabel('Rest Wavelength (\u03BCm)')
        axwl.set_xlim([3e8/i for i in ax.get_xlim()])
        axwl.set_xscale('log')
        axwl.xaxis.set_major_locator(LogLocator(base=10, numticks=10))
        axwl.xaxis.set_major_formatter(FuncFormatter(lambda x,pos : int(1e6*x)))
        # Show plot
        plt.show()
    else:  # If a given image doesn't exist
        logging.warning("A target was specified, but no SED data was found, so nothing was done.")
