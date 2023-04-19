import sys
import os
sys.path.insert(0, './lib')
sys.path.insert(0, './param_files')

import json

import get_ned_data as gnd

# Prompt to change working directory
print(' ')
workdir = input("Enter working directory: ")
if workdir[-1] != '/':
    workdir += '/'
sed_data_loc = workdir + "sed_data.json"

# Prompt to specify targets for photometry
print(' ')
target = input("Enter the desired target : ")

# Get SED data or create new file
try:
    new_sed_data = gnd.get_sed_data(target)
    all_sed_data = json.load(open(sed_data_loc))
    sed_outfile = open(sed_data_loc, 'w')
    if target not in all_sed_data:
        all_sed_data[target] = new_sed_data
    else:
        # Remove any previously-found NED points
        bandremoval = []
        freqremoval = []
        for band in all_sed_data[target]["sed_flags"]:
            if "n" in all_sed_data[target]["sed_flags"][band]:
                bandremoval.append(band)
        for key in bandremoval:
            del all_sed_data[target]["sed_flags"][key]
        for freq in all_sed_data[target]["sed_filternames"]:
            if all_sed_data[target]["sed_filternames"][freq] in bandremoval:
                freqremoval.append(freq)
                all_sed_data[target]["sed_data"].pop(freq, None)
                all_sed_data[target]["sed_unc_upper"].pop(freq, None)
                all_sed_data[target]["sed_unc_lower"].pop(freq, None)
                all_sed_data[target]["sed_telescopenames"].pop(freq, None)
        for key in freqremoval:
            all_sed_data[target]["sed_filternames"].pop(key, None)
        # Put in new NED points
        for freq in new_sed_data["sed_data"]:
            tempfreq = float(freq)
            while tempfreq in all_sed_data[target]["sed_data"]:
                tempfreq += 1
            all_sed_data[target]["sed_data"][str(tempfreq)] = new_sed_data["sed_data"][freq]
            all_sed_data[target]["sed_unc_upper"][str(tempfreq)] = new_sed_data["sed_unc_upper"][freq]
            all_sed_data[target]["sed_unc_lower"][str(tempfreq)] = new_sed_data["sed_unc_lower"][freq]
            all_sed_data[target]["sed_telescopenames"][str(tempfreq)] = new_sed_data["sed_telescopenames"][freq]
            all_sed_data[target]["sed_filternames"][str(tempfreq)] = new_sed_data["sed_filternames"][freq]
        if "sed_flags" not in all_sed_data[target]:
            all_sed_data[target]["sed_flags"] = {}
        all_sed_data[target]["sed_flags"] |= new_sed_data["sed_flags"]
    json.dump(all_sed_data, sed_outfile, indent=5)
except FileNotFoundError:
    print("\nSED data file not found ... creating a new one here")
    sed_outfile = open(sed_data_loc, 'w')
    json.dump({target:gnd.get_sed_data(target)}, sed_outfile, indent=5)
