import os

import json

import logging

#######################
#######################
## Start main script ##
#######################
#######################

# Prompt to change working directory
print(' ')
workdir = input("Enter working directory: ")
if workdir[-1] != '/':
    workdir += '/'
sed_data_loc = workdir + "sed_data.json"

# Prompt to specify target
print(' ')
target = input("Enter the desired target : ")

# Check if target is valid
subdirs = [f.name for f in os.scandir(workdir) if f.is_dir()]
if target in subdirs:

    # Get filter list
    fltr_list = [(f.name).replace('.fits','') for f in os.scandir(workdir + target + '/fits/')]

    # Import necessary static data
    tele_wl = json.load(open('./param_files/telescope_wavelengths.json'))
    tele_filters = json.load(open('./param_files/telescope_filter_names.json'))

    # Set up data structure (retrieve if it exists)
    sed_data = {}
    try:
        sed_data = json.load(open(sed_data_loc))
        if target not in sed_data:
            sed_data[target] = {
                        "sed_flux" : {},
                        "sed_freq" : {},
                        "sed_unc_lower" : {},
                        "sed_unc_upper" : {},
                        "sed_telescopenames" : {},
                        "sed_flags" : {}
                }
    except FileNotFoundError:
        logging.warning("No SED data file found in the specified location . . . exiting")
        exit()

    # List filters and prompt for a specific filter
    print(' ')
    print('Filters for ' + target + ':')
    for key in sed_data[target]["sed_flux"]:
        print('\t' + key)
    if len(sed_data[target]["sed_flux"]) == 0:
        print('\tnone found')
    print(' ')
    fltr = input('Choose a filter to manually adjust (include suffix code; or enter a new one) : ')

    # Check if filter exists
    if fltr in sed_data[target]["sed_flux"]:

        # Enter new flux value
        print(' ')
        sed_data[target]["sed_flux"][fltr] = float(input(
            'Enter a new flux value (previously ' 
            + str(sed_data[target]["sed_flux"][fltr]) + ') : '
            ))

        # Enter new upper uncertainty
        print(' ')
        sed_data[target]["sed_unc_upper"][fltr] = float(input(
            'Enter a new upper uncertainty (previously ' 
            + str(sed_data[target]["sed_unc_upper"][fltr]) + ') : '
            ))

        # Enter new lower uncertainty
        print(' ')
        sed_data[target]["sed_unc_lower"][fltr] = float(input(
            'Enter a new lower uncertainty (previously ' 
            + str(sed_data[target]["sed_unc_lower"][fltr]) + ') : '
            ))

        # Write data
        sed_outfile = open(sed_data_loc, 'w')
        json.dump(sed_data, sed_outfile, indent=5)

    else:  # If the filter does not exist
        logging.warning("Specified filter does not exist . . . exiting")
        exit()

else:  # If the target is not valid
    logging.warning("A target was specified, but no corresponding folder was found")
