import pandas as pd
import numpy as np
import json

properties = json.load(open('/home/ben/Desktop/research/research_boizelle_working/kinemetry_working/output/for_paper/all_parameters.json'))

outputfile_loc = "/home/ben/Desktop/co_props.tex"

# Define rounding function
def round_to_n(x,n=3):
    if x == 0.0 or x is None:
        return 0
    r = round(x, -int( np.floor(np.log10(abs(x))) ) + (n-1))
    if r > 10**(n-1):
        r = int(r)
    return r

itr = 0
target_ms = ""
prev_target_name = ""
for key in properties:
    itr += 1

    target_name = key.split('_')[0].replace("GC0", "GC")
    target_name = target_name.replace("GC", "GC ")
    project_code = key.split('_')[1]
    line_name = key.split('_')[2]
    if line_name == "CO10":
        line_name = "\\coone"
    elif line_name == "CO21":
        line_name = "\\cotwo"
    elif line_name == "CO32":
        line_name = "\\cothree"
    
    if prev_target_name == target_name:
        target_ms += "--- & " +\
                    line_name + " & " +\
                    str(int(round(properties[key]['disk radius (pc)'],0))) + " & " +\
                    str(int(round(2*properties[key]['max k1 (km/s)'],0))) + " & " +\
                    str(round(properties[key]['average pa (deg)'],1)) + "/" + str(round(properties[key]['range pa (deg)'],1)) + " & " +\
                    str(round(properties[key]['average q'],2)) + "/" + str(round(properties[key]['range q'],2)) + " & " +\
                    str(round_to_n(properties[key]['intensity (Jy km/s)'],4)) + " (" + str(round_to_n(properties[key]['intensity uncertainty (Jy km/s)'],3)) + ") & " +\
                    str(round_to_n(1e-6*properties[key]['luminosity (K km s^-1 pc^2)'],3)) + " (" + str(round_to_n(1e-6*properties[key]['luminosity uncertainty (K km s^-1 pc^2)'],2)) + ") & " +\
                    str(round_to_n(1e-7*properties[key]['gas mass (M_sol)'],3)) + " (" + str(round_to_n(1e-7*properties[key]['gas mass uncertainty (M_sol)'],3)) + ")"
    else:
        target_ms += target_name + " & " +\
                    line_name + " & " +\
                    str(int(round(properties[key]['disk radius (pc)'],0))) + " & " +\
                    str(int(round(2*properties[key]['max k1 (km/s)'],0))) + " & " +\
                    str(round(properties[key]['average pa (deg)'],1)) + "/" + str(round(properties[key]['range pa (deg)'],1)) + " & " +\
                    str(round(properties[key]['average q'],2)) + "/" + str(round(properties[key]['range q'],2)) + " & " +\
                    str(round_to_n(properties[key]['intensity (Jy km/s)'],4)) + " (" + str(round_to_n(properties[key]['intensity uncertainty (Jy km/s)'],3)) + ") & " +\
                    str(round_to_n(1e-6*properties[key]['luminosity (K km s^-1 pc^2)'],3)) + " (" + str(round_to_n(1e-6*properties[key]['luminosity uncertainty (K km s^-1 pc^2)'],2)) + ") & " +\
                    str(round_to_n(1e-7*properties[key]['gas mass (M_sol)'],3)) + " (" + str(round_to_n(1e-7*properties[key]['gas mass uncertainty (M_sol)'],3)) + ")"

    if itr < len(properties):
        target_ms += " \\\\\n"

    prev_target_name = target_name

print(target_ms, file=open(outputfile_loc,mode='w'))
