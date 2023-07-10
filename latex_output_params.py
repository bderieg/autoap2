import json
import numpy as np
import re

# Set some constants
prefixes = [
        "ESO",
        "NGC",
        "FCC",
        "GAMA",
        "CGCG",
        "Hydra",
        "IC",
        "PGC",
        "PKS",
        "UGC",
        "MCG",
        "WISEA"
        ]
targets_trunc = {
        "WISEAJ090330.81-010612.3" : "WISEA J090",
        "WISEAJ130100.81+270131.4" : "WISEA J130",
        "WISEAJ133333.03+261619.4" : "WISEA J133",
        "WISEAJ141611.82+015205.1" : "WISEA J141"
        }


# Define rounding function
def round_to_n(x,n=3):
    if x == 0.0 or x is None:
        return 0
    r = round(x, -int( np.floor(np.log10(abs(x))) ) + (n-1))
    if r > 10**(n-1):
        r = int(r)
    return r

# Prompt to change working directory
print(' ')
workdir = input("Enter working directory (relative or absolute) : ")
if workdir[-1] != '/':
    workdir += '/'
sed_data_loc = workdir + "sed_data.json"
fit_data_loc = workdir + "fit_data.json"
kin_data_loc = workdir + "all_parameters_link.json"

# Prompt for export directory
print(' ')
exportdir = input("Enter desired table output directory (relative or absolute) : ")
if exportdir[-1] != '/':
    exportdir += '/'
outputfile_loc = exportdir + "fit_parameters.tex"

# Open SED data file
sed_data = json.load(open(sed_data_loc))
fit_data = json.load(open(fit_data_loc))
kin_data = json.load(open(kin_data_loc))

# Write LaTeX file
with open(outputfile_loc,"w") as f:
    for target in sorted(sed_data):

        # Set up initial target string
        target_ms = ""

        # Target name
        pred = False
        if target in targets_trunc:
            target_ms += "\\text{"+targets_trunc[target]+"}\\footnote{"+target+"} &"
        else:
            for p in prefixes:
                if p in target:
                    pred = True
                    target_ms += "\\text{"+re.sub(r'('+p+')', r'\1 ', target)+"} &"
                    break
            if pred is False:
                target_ms += "\\text{"+target+"} &"

        # Intensity
        if target in kin_data:
            target_ms += " $" +\
                    str(round_to_n(kin_data[target]['intensity (Jy km/s)'])).rstrip('.') +\
                    " (" +\
                    str(round_to_n(kin_data[target]['intensity uncertainty (Jy km/s)'])).rstrip('.') +\
                    ")$ &"
        else:
            target_ms += " \\nodata &"

        # Gas mass
        if target in kin_data:
            target_ms += " $" +\
                    str(round_to_n(kin_data[target]['gas mass (M_sol)']/1e6)).rstrip('.') +\
                    " (" +\
                    str(round_to_n(kin_data[target]['gas mass uncertainty (M_sol)']/1e6)).rstrip('.') +\
                    ")$ &"
        else:
            target_ms += " \\nodata &"

        # Dust mass
        if target in fit_data:
            target_ms += " $" +\
                    str(round_to_n(fit_data[target]['mb']['mass']/1e6)).rstrip('.') +\
                    " (" +\
                    str(round_to_n(fit_data[target]['mb']['mass_unc']/1e6)).rstrip('.') +\
                    ")$ &"
        else:
            target_ms += " \\nodata &"

        # Mass ratio
        if target in kin_data and target in fit_data:
            target_ms += " $" +\
                    str(round_to_n(kin_data[target]['gas mass (M_sol)']/fit_data[target]['mb']['mass'])).rstrip('.') +\
                    " (" +\
                    str(round_to_n(kin_data[target]['gas mass (M_sol)']/fit_data[target]['mb']['mass'] * np.sqrt((kin_data[target]['gas mass uncertainty (M_sol)']/kin_data[target]['gas mass (M_sol)'])**2+(fit_data[target]['mb']['mass_unc']/fit_data[target]['mb']['mass'])**2))).rstrip('.') +\
                    ")$ &"
        else:
            target_ms += " \\nodata &"

        # Dust temperature
        if target in fit_data:
            target_ms += " $" +\
                    str(round_to_n(fit_data[target]['mb']['temperature'])).rstrip('.') +\
                    " (" +\
                    str(round_to_n(fit_data[target]['mb']['temperature_unc'])).rstrip('.') +\
                    ")$ &"
        else:
            target_ms += " \\nodata &"

        # β
        if target in fit_data:
            if fit_data[target]['mb']['beta_unc'] == 0.0:
                target_ms += " $" +\
                        str(round_to_n(fit_data[target]['mb']['beta'])).rstrip('.') +\
                        "$ &"
            else:
                target_ms += " $" +\
                        str(round_to_n(fit_data[target]['mb']['beta'])).rstrip('.') +\
                        " (" +\
                        str(round_to_n(fit_data[target]['mb']['beta_unc'])).rstrip('.') +\
                        ")$ &"
        else:
            target_ms += " \\nodata &"

        # S0
        if target in fit_data:
            if 'stellar' in fit_data[target]:
                target_ms += " $" +\
                        str(round_to_n(fit_data[target]['stellar']['coef']*1e29)).rstrip('.') +\
                        "$"
            else:
                target_ms += " \\nodata"
        else:
            target_ms += " \\nodata"

        target_ms += " \\\\"
        print(target_ms, file=f)
