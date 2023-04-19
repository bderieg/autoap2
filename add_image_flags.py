import os
import logging

import json

##################
# Types of flags #
##################

flag_names = [
        "u (to be calculated as upper limit)",
        "r (unresolved contamination present)",
        "n (from NED)"
        ]
flag_names_short = [
        "u",
        "r",
        "n"
        ]

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

# Prompt to specify targets for photometry
print(' ')
target = input("Enter the desired target : ")

# Check if target is valid
subdirs = [f.name for f in os.scandir(workdir) if f.is_dir()]
if target in subdirs:

    # Get filter list
    fltr_list = [(f.name).replace('.fits','') for f in os.scandir(workdir + target + '/fits/')]

    # List filters and prompt for a specific filter
    print(' ')
    print('Filters for ' + target + ':')
    for f in fltr_list:
        print('\t' + f)
    print(' ')
    fltr = input('Choose a filter to add a flag to : ')

    # Check if if filter is valid
    if fltr in fltr_list:
        print(' ')
        print('Here is a list of possible flags: ')
        for fl in flag_names:
            print('\t' + fl)
        print(' ')
        flag = input('Choose one of those flags to add : ')
        if flag not in flag_names_short:
            logging.warning("That's not a valid flag . . . exiting")
            exit()

        # Get SED data or create new file
        try:
            all_sed_data = json.load(open(sed_data_loc))
            sed_outfile = open(sed_data_loc, 'w')
            if target not in all_sed_data:
                all_sed_data[target] = {}
            if "sed_flags" not in all_sed_data[target]:
                all_sed_data[target]["sed_flags"] = {}
            if fltr in all_sed_data[target]["sed_flags"]:
                flagset = set(all_sed_data[target]["sed_flags"][fltr])
                flagset.add(flag)
                outstring = ''
                for i in flagset:
                    outstring += i
                all_sed_data[target]["sed_flags"][fltr] = outstring
            else:
                all_sed_data[target]["sed_flags"][fltr] = flag
            json.dump(all_sed_data, sed_outfile, indent=5)
        except FileNotFoundError:
            print("\nSED data file not found ... creating a new one here")
            sed_outfile = open(sed_data_loc, 'w')
            json.dump({target:{"sed_data":{fltr:flag}}}, sed_outfile, indent=5)
        
    else:  # If the filter is not valid
        logging.warning("Specified filter does not exist . . . exiting")
        exit()

else:  # If the target is not valid
    logging.warning("A target was specified, but no corresponding folder was found, so no flags could be added")
