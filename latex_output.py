import json
import numpy as np
import re

# Set some constants
band_contains = [
        "ALMAExtended",
        "ALMAExtended",
        "SPIRE500",
        "SPIRE350",
        "SPIRE250",
        "PACS160",
        "PACS100",
        "PACS70"
        ]
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

# Define rounding function
def round_to_n(x,n=3):
    if x == 0.0 or x is None:
        return 0
    r = round(x, -int( np.floor(np.log10(abs(x))) ) + (n-1))
    if r > 10**(n-1):
        r = int(r)
    return r

def is_close(a,b):
    if (a < 1.05*b) and (a > 0.95*b):
        return True
    else:
        return False

# Prompt to change working directory
print(' ')
workdir = input("Enter working directory (relative or absolute) : ")
if workdir[-1] != '/':
    workdir += '/'
sed_data_loc = workdir + "sed_data.json"

# Prompt for export directory
print(' ')
exportdir = input("Enter desired table output directory (relative or absolute) : ")
if exportdir[-1] != '/':
    exportdir += '/'
outputfile_loc = exportdir + "new_measurements.tex"

# Open SED data file
sed_data = json.load(open(sed_data_loc))

# Write LaTeX file
with open("/home/ben/Desktop/new_measurements.tex","w") as f:
    for target in sorted(sed_data):
        target_ms = ""
        pred = False
        for p in prefixes:
            if p in target:
                pred = True
                target_ms += "\\text{"+re.sub(r'('+p+')', r'\1 ', target)+"} &"
                break
        if pred is False:
            target_ms += "\\text{"+target+"} &"
        for band in band_contains:
            contains = False
            for key in sed_data[target]["sed_flux"]:
                almaitr = 0
                if band in key:
                    flux = str(round_to_n(sed_data[target]["sed_flux"][key]))
                    if band == "ALMAExtended":
                        freq = sed_data[target]["sed_freq"][key]
                        if not is_close(2.5e9,freq) and almaitr == 0:
                            continue
                        if not is_close(3.6e9,freq) and almaitr == 1:
                            continue
                        flux = str(round_to_n(1e3*float(flux)))
                        almaitr += 1
                    unc = ""
                    unc_upper = str(round_to_n(sed_data[target]["sed_unc_upper"][key],n=2))
                    unc_lower = str(round_to_n(sed_data[target]["sed_unc_lower"][key],n=2))
                    if unc_upper == unc_lower:
                        unc = unc_upper
                    else:
                        unc = "^{+"+unc_upper+"}_{-"+unc_lower+"}"
                    if float(unc_upper) < 0.0:
                        target_ms += " $<"+flux+"$ &"
                    else:
                        if "+" in unc:
                            target_ms += " $"+flux+" ("+unc+")$ &"
                        elif not is_close(float(unc),0.1*float(flux)) and "SPIRE" in key:
                            target_ms += " $"+flux+" ("+unc+")$ &"
                        elif not is_close(float(unc),0.07*float(flux)) and "PACS" in key:
                            target_ms += " $"+flux+" ("+unc+")$ &"
                        else:
                            target_ms += " $"+flux+"$ &"
                    contains = True
                    break
            if contains is False:
                target_ms += " &"
        target_ms = target_ms[:-2]
        target_ms += " \\\\"
        print(target_ms, file=f)
