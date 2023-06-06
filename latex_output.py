import json
import numpy as np
import re

# Set some constants
band_contains = [
        "ALMANuclear",
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

def is_close(a,b):
    if (a < 1.08*b) and (a > 0.92*b):
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
outputfile_loc = exportdir + "new_measurements_table.tex"

# Open SED data file
sed_data = json.load(open(sed_data_loc))

# Write LaTeX file
with open(outputfile_loc,"w") as f:
    for target in sorted(sed_data):
        target_ms = ""
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

        almaband = 0
        almakeylist = [key for key in sed_data[target]["sed_freq"] if "ALMA" in key]
        almakey = ""
        if len(almakeylist) > 0:
            almakey = almakeylist[0]
            almafreq = sed_data[target]["sed_freq"][almakey]
            if almafreq > 84e9 and almafreq < 116e9:
                almaband = 3
            elif almafreq > 125e9 and almafreq < 163e9:
                almaband = 4
            elif almafreq > 163e9 and almafreq < 211e9:
                almaband = 5
            elif almafreq > 211e9 and almafreq < 275e9:
                almaband = 6
            elif almafreq > 275e9 and almafreq < 373e9:
                almaband = 7
            elif almafreq > 385e9 and almafreq < 500e9:
                almaband = 8
            elif almafreq > 602e9 and almafreq < 720e9:
                almaband = 9
            elif almafreq > 787e9 and almafreq < 950e9:
                almaband = 10
        else:
            almaband = "\\nodata"

        target_ms += " " + str(almaband) + " &"

        for band in band_contains:
            contains = False
            for key in sed_data[target]["sed_flux"]:
                if band in key:

                    flux = "{:f}".format(round_to_n(sed_data[target]["sed_flux"][key])).rstrip('0').rstrip('.')
                    unc = ""
                    unc_upper = "{:f}".format(round_to_n(sed_data[target]["sed_unc_upper"][key],n=2)).rstrip('0').rstrip('.')
                    unc_lower = "{:f}".format(round_to_n(sed_data[target]["sed_unc_lower"][key],n=2)).rstrip('0').rstrip('.')
                    if unc_upper == unc_lower:
                        unc = unc_upper
                    else:
                        unc = "^{+"+unc_upper+"}_{-"+unc_lower+"}"
                    if float(unc_upper) < 0.0:
                        target_ms += " $<"+flux+"$ &"
                    else:
                        if "+" in unc:
                            target_ms += " $"+flux+" ("+unc+")$ &"
                        elif not is_close(float(unc),0.1*float(flux)) and ("SPIRE" in key or "ALMA" in key):
                            target_ms += " $"+flux+" ("+unc+")$ &"
                        elif not is_close(float(unc),0.07*float(flux)) and "PACS" in key:
                            target_ms += " $"+flux+" ("+unc+")$ &"
                        else:
                            target_ms += " $"+flux+"$ &"
                    contains = True
                    break
            if contains is False:
                target_ms += " \\nodata &"
        target_ms = target_ms[:-2]
        target_ms += " \\\\"
        print(target_ms, file=f)
