import json
import os
from termcolor import colored

# Prompt to change working directory
print(' ')
workdir = input("Enter working directory (relative or absolute) : ")
if workdir[-1] != '/':
    workdir += '/'
sed_data_loc = workdir + "sed_data.json"

# Prompt to specify targets for photometry
print(' ')
curinp = input("Combine ALMA measurements from all targets in this directory? (y/n) : ")
subdirs = [f.name for f in os.scandir(workdir) if f.is_dir()]
if curinp=='y' or curinp=='Y':
    target_names = subdirs
else:
    print(' ')
    target_names = input("Enter (as a comma-separated list) the desired targets : ")
    target_names = target_names.split(',')

sed_data = json.load(open(sed_data_loc))

for target in target_names:

    print(' ')
    print(colored('Combining extended ALMA measurements for : '+target, 'green'))
    print(' ')

    # Find extended keys
    extended_keys_nat = []
    extended_keys_bri = []
    extended_keys_other = []
    for key in sed_data[target]['sed_flux']:
        if "Extended" in key and "Bri" in key and ("Upper" in key or "Lower" in key):
            extended_keys_bri.append(key)
        elif "Extended" in key and "Nat" in key and ("Upper" in key or "Lower" in key):
            extended_keys_nat.append(key)
        elif "Extended" in key and ("Upper" in key or "Lower" in key):
            extended_keys_other.append(key)

    # Get combined flux and uncertainty (natural weighting)
    lower_flux_nat = 0.0
    lower_unc_nat = 0.0
    upper_unc_nat = 0.0
    for key in extended_keys_nat:
        if "Lower" in key:
            lower_flux_nat = sed_data[target]['sed_flux'][key]
            lower_unc_nat = sed_data[target]['sed_unc_lower'][key]
    for key in extended_keys_nat:
        if "Upper" in key:
            upper_unc_nat = (sed_data[target]['sed_flux'][key] - lower_flux_nat) + sed_data[target]['sed_unc_upper'][key]

    # Get combined flux and uncertainty (Briggs weighting)
    lower_flux_bri = 0.0
    lower_unc_bri = 0.0
    upper_unc_bri = 0.0
    for key in extended_keys_bri:
        if "Lower" in key:
            lower_flux_bri = sed_data[target]['sed_flux'][key]
            lower_unc_bri = sed_data[target]['sed_unc_lower'][key]
    for key in extended_keys_bri:
        if "Upper" in key:
            upper_unc_bri = (sed_data[target]['sed_flux'][key] - lower_flux_bri) + sed_data[target]['sed_unc_upper'][key]

    # Get combined flux and uncertainty (unknown weighting)
    lower_flux_other = 0.0
    lower_unc_other = 0.0
    upper_unc_other = 0.0
    for key in extended_keys_other:
        if "Lower" in key:
            lower_flux_other = sed_data[target]['sed_flux'][key]
            lower_unc_other = sed_data[target]['sed_unc_lower'][key]
    for key in extended_keys_other:
        if "Upper" in key:
            upper_unc_other = (sed_data[target]['sed_flux'][key] - lower_flux_other) + sed_data[target]['sed_unc_upper'][key]

    # Make new extended keys
    try:
        sed_data[target]['sed_freq']['ALMAExtendedCombinedNat'] = sed_data[target]['sed_freq'][extended_keys_nat[0]]
        sed_data[target]['sed_flux']['ALMAExtendedCombinedNat'] = lower_flux_nat
        sed_data[target]['sed_unc_upper']['ALMAExtendedCombinedNat'] = upper_unc_nat
        sed_data[target]['sed_unc_lower']['ALMAExtendedCombinedNat'] = lower_unc_nat
        sed_data[target]['sed_telescopenames']['ALMAExtendedCombinedNat'] = "ALMA Extended"
    except IndexError:
        pass

    try:
        sed_data[target]['sed_freq']['ALMAExtendedCombinedBri'] = sed_data[target]['sed_freq'][extended_keys_bri[0]]
        sed_data[target]['sed_flux']['ALMAExtendedCombinedBri'] = lower_flux_bri
        sed_data[target]['sed_unc_upper']['ALMAExtendedCombinedBri'] = upper_unc_bri
        sed_data[target]['sed_unc_lower']['ALMAExtendedCombinedBri'] = lower_unc_bri
        sed_data[target]['sed_telescopenames']['ALMAExtendedCombinedBri'] = "ALMA Extended"
    except IndexError:
        pass

    try:
        sed_data[target]['sed_freq']['ALMAExtendedCombined'] = sed_data[target]['sed_freq'][extended_keys_other[0]]
        sed_data[target]['sed_unc_upper']['ALMAExtendedCombined'] = upper_unc_other
        sed_data[target]['sed_flux']['ALMAExtendedCombined'] = lower_flux_other
        sed_data[target]['sed_unc_lower']['ALMAExtendedCombined'] = lower_unc_other
        sed_data[target]['sed_telescopenames']['ALMAExtendedCombined'] = "ALMA Extended"
    except IndexError:
        pass

    # Delete the old extended keys
    for key1 in sed_data[target]:
        for key2 in extended_keys_nat:
            sed_data[target][key1].pop(key2,None)
        for key2 in extended_keys_bri:
            sed_data[target][key1].pop(key2,None)
        for key2 in extended_keys_other:
            sed_data[target][key1].pop(key2,None)


# Write out new data
sed_outfile = open(sed_data_loc, 'w')
json.dump(sed_data, sed_outfile, indent=5)
