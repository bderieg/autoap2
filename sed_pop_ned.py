import sys
import os
sys.path.insert(0, './lib')
sys.path.insert(0, './param_files')

import json
from termcolor import colored

import get_ned_data as gnd

# Prompt to change working directory
print(' ')
workdir = input("Enter working directory: ")
if workdir[-1] != '/':
    workdir += '/'
sed_data_loc = workdir + "sed_data.json"

# Prompt to specify targets for photometry
print(' ')
allq = input("Gather data from NED for all targets in this folder? (y/N): ")
subdirs = [f.name for f in os.scandir(workdir) if f.is_dir()]
if allq == 'y' or allq == 'Y':
    target_names = subdirs
else:
    print(' ')
    target_names = input("Enter (as a comma-separated list) the desired target(s) : ")
    target_names = target_names.split(',')

# Get SED data or create new file
for target in target_names:
    try:
        
        print(' ')
        print(colored('Attempting to gather data for ' + target,'green'))
        print(' ')

        # Get SED data from NED
        new_sed_data = gnd.get_sed_data(target)

        # Get previously-existing local SED data
        sed_data = json.load(open(sed_data_loc))
        sed_outfile = open(sed_data_loc, 'w')

        # Add target if it's not there
        if target not in sed_data:
            sed_data[target] = {
                        "sed_flux" : {},
                        "sed_freq" : {},
                        "sed_unc_upper" : {},
                        "sed_unc_lower" : {},
                        "sed_telescopenames" : {},
                        "sed_flags" : {}
                    }

        # Remove previous NED measurements from local data
        if "sed_flux" in sed_data[target]:
            bands_to_remove = []
            for band in sed_data[target]["sed_flux"]:
                if band in sed_data[target]["sed_flags"]:
                    if "n" in sed_data[target]["sed_flags"][band]:
                        bands_to_remove.append(band)
            for band in bands_to_remove:
                sed_data[target]["sed_flux"].pop(band, None)
                sed_data[target]["sed_freq"].pop(band, None)
                sed_data[target]["sed_unc_lower"].pop(band, None)
                sed_data[target]["sed_unc_upper"].pop(band, None)
                sed_data[target]["sed_telescopenames"].pop(band, None)
                sed_data[target]["sed_flags"].pop(band, None)
        else:
            bands_to_remove = []
            for band in sed_data[target]["sed_flags"]:
                if "n" in sed_data[target]["sed_flags"][band]:
                    bands_to_remove.append(band)
            for band in bands_to_remove:
                sed_data[target]["sed_flags"].pop(band, None)
            sed_data[target]["sed_flux"] = {}
            sed_data[target]["sed_freq"] = {}
            sed_data[target]["sed_unc_lower"] = {}
            sed_data[target]["sed_unc_upper"] = {}
            sed_data[target]["sed_telescopenames"] = {}

        # Merge the NED data with the local data
        for band in new_sed_data["sed_flux"]:
            fltr_itr = 0
            while (band+str(fltr_itr)) in sed_data[target]["sed_flux"]:
                fltr_itr += 1
            sed_data[target]["sed_flux"][band+str(fltr_itr)] = new_sed_data["sed_flux"][band]
            sed_data[target]["sed_freq"][band+str(fltr_itr)] = new_sed_data["sed_freq"][band]
            sed_data[target]["sed_unc_lower"][band+str(fltr_itr)] = new_sed_data["sed_unc_lower"][band]
            sed_data[target]["sed_unc_upper"][band+str(fltr_itr)] = new_sed_data["sed_unc_upper"][band]
            sed_data[target]["sed_telescopenames"][band+str(fltr_itr)] = new_sed_data["sed_telescopenames"][band]
            sed_data[target]["sed_flags"][band+str(fltr_itr)] = new_sed_data["sed_flags"][band]

        # Write
        json.dump(sed_data, sed_outfile, indent=5)
        sed_outfile.close()

    except FileNotFoundError:
        print("\nSED data file not found ... creating a new one here")
        sed_outfile = open(sed_data_loc, 'w')
        json.dump({target:gnd.get_sed_data(target)}, sed_outfile, indent=5)
        sed_outfile.close()

    # except:
    #     print(colored("Something bad happened (╯°□°）╯︵ ┻━┻","red"))
    #     exit()
