import sys
import os
sys.path.insert(0, './lib')

import lmfit
import numpy as np
import json
import astropy.units as u
import pandas as pd

import physical_functions as pf
from read_params import read_params as rp
import mcmc_functions as mcmcf
import mpfit

# Define some constants

fit_params = {
        'mb_mass' : 1e6,
        'mb_mass_uplim' : 1e20,
        'mb_mass_lolim' : 1.0,
        'mb_mass_hold' : False,
        'mb_temp' : 25,
        'mb_temp_uplim' : 25.0,
        'mb_temp_lolim' : 0.0,
        'mb_temp_hold' : False,
        'mb_beta' : 2.0,
        'mb_beta_uplim' : 16.0,
        'mb_beta_lolim' : 0.0,
        'mb_beta_hold' : False,
        'verbose' : 0
        }


###############################
###############################
# Import data with user input #
###############################
###############################

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
target = input("Enter the desired target : ")

# Import fit parameters file
try:
    fit_params |= rp(workdir+target+'/sed_fit.param')
except FileNotFoundError:
    print(' ')
    print('fit parameters file not found . . . exiting')
    exit()
for key in fit_params:
    if "hold" in key:
        if fit_params[key] == False:
            fit_params[key] = 0
        elif fit_params[key] == True:
            fit_params[key] = 1
assert 'points' in fit_params, "\'points\' keyword not specified in parameter file"

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

points = fit_params['points']
freq = [1e-9*freq[i] for i in points]  # Also convert to GHz
flux = [flux[i] for i in points]
unc_upper = [unc_upper[i] for i in points]
unc_lower = [unc_lower[i] for i in points]
tele_names = [tele_names[i] for i in points]

print(' ')
print("Fitting the following points:")
for i in tele_names:
    print("\t"+i)

############
############
# Fit data #
############
############

if fit_params['verbose'] >= 1:
    quiet = False

#######################
# Fit radio power law #
#######################

##########################
# Fit modified blackbody #
##########################

# Create fitted object
pv = [
        fit_params['mb_mass'], 
        fit_params['mb_temp'], 
        fit_params['mb_beta'], 
        galaxy_properties.loc[target,'D_L (Mpc)']
    ]
p_fixed = [
        fit_params['mb_mass_hold'], 
        fit_params['mb_temp_hold'], 
        fit_params['mb_beta_hold'], 
        True
    ]
p_ltd = [
        [1,1],
        [1,1],
        [1,1],
        [0,0]
    ]
p_lims = [
        [fit_params['mb_mass_lolim'],fit_params['mb_mass_uplim']], 
        [fit_params['mb_temp_lolim'],fit_params['mb_temp_uplim']], 
        [fit_params['mb_beta_lolim'],fit_params['mb_beta_uplim']], 
        [0.,0.]
    ]
parinfo = [
                {'value':p, 'fixed':pf, 'limited':pd, 'limits':pl} 
            for p,pf,pd,pl in zip(pv,p_fixed,p_ltd,p_lims)
        ]
quiet = True
mbfit = mpfit.mpfit(
                    pf.mb_fit, 
                    functkw={
                        'x':np.array(freq), 
                        'y':np.array(flux), 
                        'err':np.array([(hi+lo)/2 for hi,lo in zip(unc_upper,unc_lower)])
                        }, 
                    parinfo=parinfo,
                    quiet=quiet
                )

#########################
# Fit stellar power law #
#########################


##################
##################
# Emcee analysis #
##################
##################

emcee_params = mcmcf.mcmc_full_run(
            np.array(freq), 
            np.array(flux), 
            np.array([(hi+lo)/2 for hi,lo in zip(unc_upper,unc_lower)]), 
            {
                "mb_mass" : mbfit.params[0],
                "mb_temp" : mbfit.params[1],
                "mb_beta" : mbfit.params[2],
                "distance" : mbfit.params[3]
            },
            {
                "mb_mass_lolim" : fit_params["mb_mass_lolim"],
                "mb_mass_uplim" : fit_params["mb_mass_uplim"],
                "mb_temp_lolim" : fit_params["mb_temp_lolim"],
                "mb_temp_uplim" : fit_params["mb_temp_uplim"],
                "mb_beta_lolim" : fit_params["mb_beta_lolim"],
                "mb_beta_uplim" : fit_params["mb_beta_uplim"]
            },
            200, 
            1000
        )

#####################
#####################
# Write fit to file #
#####################
#####################

# Open existing data file
target_fit = {}
try:
    target_fit = json.load(open(fit_data_loc))
except FileNotFoundError:
    print(' ')
    print("fit data file not found . . . creating a new one here")

# Populate with new values
target_fit[target] = {}
target_fit[target]["mb"] = {
                "mass" : mbfit.params[0],
                "mass_emcee" : emcee_params['mass'],
                "mass_unc" : mbfit.perror[0],
                "mass_unc_emcee" : emcee_params['mass_spread'],
                "temperature" : mbfit.params[1],
                "temperature_emcee" : emcee_params['temp'],
                "temperature_unc" : mbfit.perror[1],
                "temperature_unc_emcee" : emcee_params['temp_spread'],
                "beta" : mbfit.params[2],
                "beta_emcee" : emcee_params['beta'],
                "beta_unc" : mbfit.perror[2],
                "beta_unc_emcee" : emcee_params['beta_spread'],
                "distance" : mbfit.params[3]
            }

# Write
fit_outfile = open(fit_data_loc, 'w')
json.dump(target_fit, fit_outfile, indent=5)
