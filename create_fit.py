import sys
import os
sys.path.insert(0, './lib')

import lmfit
import numpy as np
import json

import physical_functions as pf


#####################
#####################
# Start main script #
#####################
#####################

############################
# Prompt user for SED data #
############################

# Prompt to change working directory
print(' ')
workdir = input("Enter working directory: ")
if workdir[-1] != '/':
    workdir += '/'
sed_data_loc = workdir + "sed_data.json"
fit_data_loc = workdir + "fit_data.json"

# Prompt to specify targets for photometry
print(' ')
target = input("Enter (as a comma-separated list) the desired targets : ")

# Import data

sed_data = json.load(open(sed_data_loc))
target_data = sed_data[target]
freq = ([ target_data['sed_freq'][key] for key in target_data['sed_freq'] ])[2:9]
flux = ([ target_data['sed_flux'][key] for key in target_data['sed_flux'] ])[2:9]

############
# Fit data #
############

# Initialize fit
bparams = lmfit.Parameters()
## Initialize fixed variables
bparams.add('distance', value=5.67765e23, vary=False)
## Initialize fitted variables
bparams.add('mass', value=1e36)
bparams.add('temperature', value=25)
bparams.add('beta', value=1.8)

# Perform fit
minner = lmfit.Minimizer(pf.mb, bparams, fcn_args=(freq,flux))
result = minner.minimize()

#####################
# Write fit to file #
#####################

# Make data structure
target_fit = {
        "NGC1380" : {
            "mb" : {
                "mass" : result.params['mass'].value,
                "mass_unc" : result.params['mass'].stderr,
                "temperature" : result.params['temperature'].value,
                "temperature_unc" : result.params['temperature'].stderr,
                "beta" : result.params['beta'].value,
                "beta_unc" : result.params['beta'].stderr,
                "distance" : result.params['distance'].value
                }
            }
        }

fit_outfile = open(fit_data_loc, 'w')
json.dump(target_fit, fit_outfile, indent=5)
