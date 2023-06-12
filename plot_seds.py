import sys
import os
sys.path.insert(0, './lib')

import pandas as pd
import numpy as np
import json
import pickle
import base64
import astropy.units as u
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter, LogLocator
from decimal import Decimal

import logging

import physical_functions as pf

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
fit_data_loc = workdir + "fit_data.json"
galaxy_properties_loc = workdir + "galaxy_properties.csv"

# Import SED file
try:
    seds = pd.read_json(sed_data_loc)
except FileNotFoundError:
    print(' ')
    print('SED data file not found . . . exiting')
    exit()

# Import fit file
try:
    fits = json.load(open(fit_data_loc))
except FileNotFoundError:
    fits = {}

# Import galaxy properties file
try:
    galaxy_properties = pd.read_csv(galaxy_properties_loc, index_col='Galaxy')
    new_index = [old_index.replace(" ","") for old_index in galaxy_properties.index]
    galaxy_properties = galaxy_properties.rename(index=dict(zip(galaxy_properties.index, new_index)))
except FileNotFoundError:
    print(' ')
    logging.warning('galaxy property data file not found . . . fits will not be plotted')

########################
# Loop through targets #
########################

sed_ps = json.load(open('./param_files/telescope_sed_pointstyles.json'))
sed_colors = json.load(open('./param_files/telescope_sed_colors.json'))

for target in target_names:
    if target in subdirs:
        # Import SED data
        sed_freq = seds[target]['sed_freq']
        sed_flux = seds[target]['sed_flux']
        sed_unc_upper = seds[target]['sed_unc_upper']
        sed_unc_lower = seds[target]['sed_unc_lower']
        sed_telenames = seds[target]['sed_telescopenames']
        sed_data_arr = np.transpose(np.asarray([[sed_freq[key], sed_flux[key]] for key in sed_flux]))
        sed_unc_lower_arr = np.transpose(np.asarray([sed_unc_lower[key] for key in sed_flux]))
        sed_unc_upper_arr = np.transpose(np.asarray([sed_unc_upper[key] for key in sed_flux]))
        sed_telenames_arr = np.transpose(np.asarray([sed_telenames[key] for key in sed_flux]))
        # Import fit data
        fitparams = {}
        try:
            for key in fits[target]:
                fitparams[key] = {}
                if "mb" in key:
                    fitparams[key]['mass'] = fits[target][key]['mass']
                    fitparams[key]['temperature'] = fits[target][key]['temperature']
                    fitparams[key]['beta'] = fits[target][key]['beta']
                    fitparams[key]['distance'] = fits[target][key]['distance']
                    fitparams[key]['posterior_spread_obj'] = fits[target][key]['posterior_spread_obj']
                elif "radio" in key:
                    fitparams[key]['radio_slope'] = fits[target][key]['slope']
                    fitparams[key]['radio_coef'] = fits[target][key]['coef']
                elif "stellar" in key:
                    fitparams[key]['stellar_slope'] = fits[target][key]['slope']
                    fitparams[key]['stellar_coef'] = fits[target][key]['coef']
        except KeyError:
            pass
        # Make unique legend list
        unique_legend_points = ["_"] * len(sed_telenames_arr)
        for ind in np.unique(np.asarray(sed_telenames_arr), return_index=True)[1]:
            unique_legend_points[ind] = sed_telenames_arr[ind] 
        # Plot data
        fig, ax = plt.subplots()
        for xval,yval,telename,legendname,unc_upper,unc_lower\
                in zip(sed_data_arr[0], sed_data_arr[1],\
                sed_telenames_arr,\
                unique_legend_points,\
                sed_unc_upper_arr, sed_unc_lower_arr):
            ax.scatter(xval, yval, marker=sed_ps[telename], color=sed_colors[telename], label=legendname, zorder=1)
            if unc_upper < 0.0:
                ax.errorbar(xval, yval, yerr=[[0.3*yval],[0]], uplims=True, fmt='none', ecolor='black', zorder=0)
            else:
                ax.errorbar(xval, yval, yerr=[[unc_lower],[unc_upper]], fmt='none', ecolor='black', capsize=0.0, zorder=0)
        # Plot fits
        try:
            # Define a frequency array
            basis = np.logspace(np.log10(min(sed_data_arr[0,:])), np.log10(max(sed_data_arr[0,:])), num=500)
            # Define an empty total curve to plot over everything
            total_curve = np.zeros_like(basis)
            # Plot all existing fits
            for key in fitparams:
                if "mb" in key:
                    # Plot mpfit curve
                    fitted_curve = pf.mb_model(
                            [
                                fitparams[key]['mass'],
                                fitparams[key]['temperature'],
                                fitparams[key]['beta'],
                                fitparams[key]['distance']
                            ], 
                            x=1e-9*basis
                            )
                    ax.plot(basis, fitted_curve, c='maroon', ls='dashdot', label='modified blackbody')
                    # Plot emcee posterior distribution
                    pso = pickle.loads(
                            base64.b64decode(fitparams[key]['posterior_spread_obj'].encode('utf-8'))
                        )
                    ax.fill_between(
                                pso[0], 
                                pso[1]-pso[2], 
                                pso[1]+pso[2], 
                                color='grey', alpha=0.5, label='1$\sigma$ Posterior Spread'
                            )
                elif "radio" in key:
                    # Plot fitted curve
                    fitted_curve = pf.pl_model(
                            1e-9*basis,
                            fitparams[key]['radio_slope'],
                            fitparams[key]['radio_coef']
                            )
                    ax.plot(basis, fitted_curve, c='blue', ls='dashed', label='power law')
                elif "stellar" in key:
                    # Plot fitted curve
                    fitted_curve = pf.pl_model(
                            1e-9*basis,
                            fitparams[key]['stellar_slope'],
                            fitparams[key]['stellar_coef']
                            )
                    ax.plot(basis, fitted_curve, c='blue', ls='dashed')
                # Update the total curve
                total_curve += fitted_curve
            # Overplot the sum of all fits
            ax.plot(basis, total_curve, c='black', ls='solid', label='total emission')
        # If no fits exist
        except KeyError:
            pass
        # Set other plot parameters
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.legend(loc='upper center', bbox_to_anchor=(0.5,-0.2), ncol=4)
        ax.set_title(target)
        ax.set_xlabel('log$_{10}$ Rest Frequency (Hz)')
        ax.set_ylabel('log$_{10}$ Flux Density (Jy)')
        ax.set_ylim([0.3*min([aa if aa > 0 else 1e10 for aa in sed_data_arr[1]]), 3.0*max(sed_data_arr[1])])
        ax.xaxis.set_major_formatter(FuncFormatter(lambda x,pos : int(np.log10(x))))
        ax.yaxis.set_major_formatter(FuncFormatter(lambda x,pos : int(np.log10(x))))
        # Set twin wavelength axis
        axwl = ax.twiny()
        axwl.set_xlabel('log$_{10}$ Rest Wavelength (m)')
        axwl.set_xlim([3e8/i for i in ax.get_xlim()])
        axwl.set_xscale('log')
        axwl.xaxis.set_major_locator(LogLocator(base=10, numticks=10))
        axwl.xaxis.set_major_formatter(FuncFormatter(lambda x,pos : int(np.log10(x))))
        # Show plot
        plt.tight_layout()
        plt.show()
    else:  # If a given image doesn't exist
        logging.warning("A target was specified, but no SED data was found, so nothing was done.")
