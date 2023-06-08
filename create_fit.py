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
import mpfit


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
    fit_params = rp(workdir+target+'/sed_fit.param')
except FileNotFoundError:
    print(' ')
    print('fit parameters file not found . . . exiting')
    exit()

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

print("Fitting the following points:")
for i in tele_names:
    print("\t"+i)

############
############
# Fit data #
############
############

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
parinfo = [
                {'value':p, 'fixed':0, 'limited':[0,0], 'limits':[0.,0.]} 
            for p,itr in zip(pv,range(len(pv)))
        ]
mbfit = mpfit.mpfit(
                    pf.mb_fit, 
                    functkw={
                        'x':np.array(freq), 
                        'y':np.array(flux), 
                        'err':np.array([(hi+lo)/2 for hi,lo in zip(unc_upper,unc_lower)])
                        }, 
                    parinfo=parinfo
                )

#########################
# Fit stellar power law #
#########################

exit()

#####################
#####################
# Write fit to file #
#####################
#####################

# Populate with new values
if fitfunc.lower() == "mb":
    itr = 0
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
