from astroquery.ipac.ned import Ned as ned
from statistics import mean
from termcolor import colored
import numpy as np

def get_ellipse_parameters(target_name):
    diameters_table = ned.get_table(target_name, table="diameters")
    a_to_b_ratios = []
    pa = []
    for itr in range(len(diameters_table)):
        if isinstance(diameters_table[itr]["Axis Ratio"], float):
            if diameters_table[itr]["Axis Ratio"] <= 1.0:
                a_to_b_ratios.append(diameters_table[itr]["Axis Ratio"])
            else:
                a_to_b_ratios.append(1/diameters_table[itr]["Axis Ratio"])
    avg_ratio = mean(a_to_b_ratios)
    return {"axis_ratio" : avg_ratio}


def get_coords(target_name):
    main_table = ned.query_object(target_name)
    return float(main_table["RA"]), float(main_table["DEC"])


def get_sed_data(target_name):
    sed_table = ned.get_table(target_name, table="photometry")

    sed = {
            "sed_flux" : {},
            "sed_freq" : {},
            "sed_unc_upper" : {},
            "sed_unc_lower" : {},
            "sed_telescopenames" : {},
            "sed_flags" : {}
        }

    band_itr = 0
    for row in sed_table:
        if row["Flux Density"] > 0.0:
            # IRAS
            if "quality flag" in row["Qualifiers"] and\
                    (("100" in row["Observed Passband"]) or ("60" in row["Observed Passband"])) and\
                    "IRAS" in row["Observed Passband"]:
                sed["sed_flux"][row["Observed Passband"]+" "+str(band_itr)] = row["Flux Density"]
                sed["sed_freq"][row["Observed Passband"]+" "+str(band_itr)] = row["Frequency"]
                sed["sed_unc_upper"][row["Observed Passband"]+" "+str(band_itr)] = row["Upper limit of uncertainty"]
                sed["sed_unc_lower"][row["Observed Passband"]+" "+str(band_itr)] = row["Lower limit of uncertainty"]
                sed["sed_telescopenames"][row["Observed Passband"]+" "+str(band_itr)] = "IRAS"
                sed["sed_flags"][row["Observed Passband"]+" "+str(band_itr)] = "n"
                band_itr += 1
            # ISO
            if "ISO" in row["Observed Passband"] and\
                    "[" not in row["Observed Passband"] and\
                    row["Frequency"] < 1e13:
                sed["sed_flux"][row["Observed Passband"]+" "+str(band_itr)] = row["Flux Density"]
                sed["sed_freq"][row["Observed Passband"]+" "+str(band_itr)] = row["Frequency"]
                sed["sed_unc_upper"][row["Observed Passband"]+" "+str(band_itr)] = row["Upper limit of uncertainty"]
                sed["sed_unc_lower"][row["Observed Passband"]+" "+str(band_itr)] = row["Lower limit of uncertainty"]
                sed["sed_telescopenames"][row["Observed Passband"]+" "+str(band_itr)] = "ISO"
                sed["sed_flags"][row["Observed Passband"]+" "+str(band_itr)] = "n"
                band_itr += 1
            # Spitzer MIPS
            if "MIPS" in row["Observed Passband"] and\
                    "Total" in row["Spatial Mode"] and\
                    row["Frequency"] < 1e13:
                sed["sed_flux"][row["Observed Passband"]+" "+str(band_itr)] = row["Flux Density"]
                sed["sed_freq"][row["Observed Passband"]+" "+str(band_itr)] = row["Frequency"]
                sed["sed_unc_upper"][row["Observed Passband"]+" "+str(band_itr)] = row["Upper limit of uncertainty"]
                sed["sed_unc_lower"][row["Observed Passband"]+" "+str(band_itr)] = row["Lower limit of uncertainty"]
                sed["sed_telescopenames"][row["Observed Passband"]+" "+str(band_itr)] = "Spitzer"
                sed["sed_flags"][row["Observed Passband"]+" "+str(band_itr)] = "n"
                band_itr += 1
            # Parkes
            if row["Refcode"] == "1990PKS90.C...0000W":
                print(colored("yee","yellow"))
                sed["sed_flux"][row["Observed Passband"]+" "+str(band_itr)] = row["Flux Density"]
                sed["sed_freq"][row["Observed Passband"]+" "+str(band_itr)] = row["Frequency"]
                sed["sed_unc_upper"][row["Observed Passband"]+" "+str(band_itr)] = row["Upper limit of uncertainty"]
                sed["sed_unc_lower"][row["Observed Passband"]+" "+str(band_itr)] = row["Lower limit of uncertainty"]
                sed["sed_telescopenames"][row["Observed Passband"]+" "+str(band_itr)] = "Parkes"
                sed["sed_flags"][row["Observed Passband"]+" "+str(band_itr)] = "n"
                band_itr += 1
            # ATCA
            if row["Refcode"] == "2010MNRAS.402.2403M":
                print(colored("yee","yellow"))
                sed["sed_flux"][row["Observed Passband"]+" "+str(band_itr)] = row["Flux Density"]
                sed["sed_freq"][row["Observed Passband"]+" "+str(band_itr)] = row["Frequency"]
                sed["sed_unc_upper"][row["Observed Passband"]+" "+str(band_itr)] = row["Upper limit of uncertainty"]
                sed["sed_unc_lower"][row["Observed Passband"]+" "+str(band_itr)] = row["Lower limit of uncertainty"]
                sed["sed_telescopenames"][row["Observed Passband"]+" "+str(band_itr)] = "ATCA"
                sed["sed_flags"][row["Observed Passband"]+" "+str(band_itr)] = "n"
                band_itr += 1

            # Deal with potentially empty boxes
            for key in sed:
                for fltr in sed[key]:
                    if sed[key][fltr] is np.ma.masked:
                        sed[key][fltr] = 0.0

    return sed
