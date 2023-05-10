import sys
import os
sys.path.insert(0, './lib')

import lmfit
import numpy as np
import json
import astropy.units as u
import pandas as pd

import physical_functions as pf


#####################
#####################
# Start main script #
#####################
#####################

################################
# Prompt user for fitting data #
################################

# Prompt to change working directory
print(' ')
workdir = input("Enter working directory (relative or absolute) : ")
if workdir[-1] != '/':
    workdir += '/'
sed_data_loc = workdir + "sed_data.json"
fit_data_loc = workdir + "fit_data.json"
galaxy_properties_loc = workdir + "galaxy_properties.csv"

# Import galaxy properties file
try:
    galaxy_properties = pd.read_csv(galaxy_properties_loc, index_col='Galaxy')
    new_index = [old_index.replace(" ","") for old_index in galaxy_properties.index]
    galaxy_properties = galaxy_properties.rename(index=dict(zip(galaxy_properties.index, new_index)))
except FileNotFoundError:
    print(' ')
    print('galaxy property data file not found . . . exiting')
    exit()

# Prompt to specify targets for photometry
print(' ')
target = input("Enter (as a comma-separated list) the desired targets : ")

# Prompt to specify type of fit
print(' ')
fitfunc = input("What function would you like to fit (\'mb\' (modified blackbody) or \'pow\' (power law)) : ")

# Prompt to specify initial values
if fitfunc.lower() == "mb":
    print(' ')
    print('Enter initial guesses for the following quantities : ')
    print (' ')
    init_mass = float( ( u.solMass*input("\tmass (in solar masses) : ") ).to(u.kg) )
    print(' ')
    tempinp = input('\tfix mass? (y/n) : ')
    vary_mass = False if tempinp.lower() == "y" else True
    print (' ')
    init_temp = float( input("\ttemperature (in kelvin) : ") )
    print(' ')
    tempinp = input('\tfix temperature? (y/n) : ')
    vary_temp = False if tempinp.lower() == "y" else True
    print (' ')
    init_beta = float( input("\tbeta : ") )
    print(' ')
    tempinp = input('\tfix beta? (y/n) : ')
    vary_beta = False if tempinp.lower() == "y" else True
elif fitfunc.lower() == "pow":
    print(' ')
    print('Enter initial guesses for the following quantities : ')
    print (' ')
    init_slope = float( input("\tslope : ") )
    print(' ')
    init_intercept = float( input('\tcoefficient : ') )
else:
    print(' ')
    print('not a valid fit function . . . exiting')
    exit()

# Prompt to subtract off other fits
print(' ')
tempinp = input('Subtract existing fits off from the data points? (Y/n) : ')
subfits = False if tempinp.lower() == "n" else True

# Prompt to delete other fits
print(' ')
tempinp = input('Delete other fits of this type? (y/N) : ')
overwrite = True if tempinp.lower() == "y" else False

# Prompt to specify points to fit
print(' ')
points = input('Specify the point numbers to fit (as a comma-separated list) : ')
points = points.split(',')
points = [ int(p) for p in points ]

# Import data, sort, and select
sed_data = json.load(open(sed_data_loc))
target_data = sed_data[target]
freq = np.array([ target_data['sed_freq'][key] for key in target_data['sed_freq'] ])
flux = np.array([ target_data['sed_flux'][key] for key in target_data['sed_flux'] ])
unc_upper = np.array([ target_data['sed_unc_upper'][key] for key in target_data['sed_unc_upper'] ])
unc_lower = np.array([ target_data['sed_unc_lower'][key] for key in target_data['sed_unc_lower'] ])
tele_names = np.array([ target_data['sed_telescopenames'][key] for key in target_data['sed_telescopenames'] ])

sort_ind = np.argsort(freq)
freq = freq[sort_ind]
flux = flux[sort_ind]
unc_upper = unc_upper[sort_ind]
unc_lower = unc_lower[sort_ind]
tele_names = tele_names[sort_ind]

freq = [freq[i] for i in points]
flux = [flux[i] for i in points]
unc_upper = [unc_upper[i] for i in points]
unc_lower = [unc_lower[i] for i in points]
tele_names = [tele_names[i] for i in points]

print("Fitting the following points:")
for i in tele_names:
    print("\t"+i)

# Import existing fit data structure
target_fit = {}
try:
    target_fit = json.load(open(fit_data_loc))
    if target not in target_fit:
        target_fit[target] = {}
except FileNotFoundError:
    target_fit[target] = {}

if overwrite:
    del_keys = []
    for key in target_fit[target]:
        if fitfunc.lower() in key:
            del_keys.append(key)
    for key in del_keys:
        del target_fit[target][key]

# Subtract off existing fits if necessary
flux_sub = np.zeros_like(flux)
for key in target_fit[target]:
    if "pow" in key:
        flux_sub += pf.powerlaw_basic(target_fit[target][key], freq)
    elif "mb" in key:
        target_fit[target][key]['distance'] = galaxy_properties['D_L (Mpc)'][target]*u.Mpc.to(u.m)
        flux_sub += pf.mb_basic(target_fit[target][key], freq)
        del target_fit[target][key]['distance']
for itr, band in zip(range(len(tele_names)), tele_names):
    if "ALMA Extended" in band:
        flux_sub[itr] = 0.0
flux -= flux_sub

############
# Fit data #
############

# Initialize fit
bparams = lmfit.Parameters()
## For modified blackbody
if fitfunc.lower() == "mb":
    # Initialize fixed variables
    bparams.add('distance', value=galaxy_properties['D_L (Mpc)'][target]*u.Mpc.to(u.m), vary=False)
    # Initialize fitted variables
    bparams.add('mass', value=init_mass, vary=vary_mass)
    bparams.add('temperature', value=init_temp, vary=vary_temp)
    bparams.add('beta', value=init_beta, vary=vary_beta)
elif fitfunc.lower() == "pow":
    # Initialize fixed variables
    ## (none)
    # Initialize fitted variables
    bparams.add('b', value=init_intercept, vary=True)
    bparams.add('alpha', value=init_slope, vary=True)

# Perform fit
minner = lmfit.Minimizer(pf.mb, bparams, fcn_args=(freq,flux,unc_upper,unc_lower))
if fitfunc.lower() == "pow":
    minner = lmfit.Minimizer(pf.powerlaw, bparams, fcn_args=(freq,flux,unc_upper,unc_lower))
result = minner.minimize()

#####################
# Write fit to file #
#####################

# Populate with new values
if fitfunc.lower() == "mb":
    itr = 0
    while "mb "+str(itr) in target_fit[target]:
        itr += 1
    try:
        target_fit[target]["mb "+str(itr)] = {
                        "mass" : result.params['mass'].value*u.kg.to(u.solMass),
                        "mass_unc" : result.params['mass'].stderr*u.kg.to(u.solMass),
                        "temperature" : result.params['temperature'].value,
                        "temperature_unc" : result.params['temperature'].stderr,
                        "beta" : result.params['beta'].value,
                        "beta_unc" : result.params['beta'].stderr
                    }
    except TypeError:
        target_fit[target]["mb "+str(itr)] = {
                        "mass" : result.params['mass'].value*u.kg.to(u.solMass),
                        "mass_unc" : np.nan,
                        "temperature" : result.params['temperature'].value,
                        "temperature_unc" : result.params['temperature'].stderr,
                        "beta" : result.params['beta'].value,
                        "beta_unc" : result.params['beta'].stderr
                    }
elif fitfunc.lower() == "pow":
    itr = 0
    while "pow "+str(itr) in target_fit[target]:
        itr += 1
    target_fit[target]["pow "+str(itr)] = {
                    "alpha" : result.params['alpha'].value,
                    "alpha_unc" : result.params['alpha'].stderr,
                    "b" : result.params['b'].value,
                    "b_unc" : result.params['b'].stderr
                }

# Write
fit_outfile = open(fit_data_loc, 'w')
json.dump(target_fit, fit_outfile, indent=5)
