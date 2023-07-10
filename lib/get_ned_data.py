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
        newrow = False
        # IRAS
        if "quality flag" in row["Qualifiers"] and\
                (("100" in row["Observed Passband"]) or ("60" in row["Observed Passband"]) or ("25" in row["Observed Passband"])) and\
                "IRAS" in row["Observed Passband"]:
            sed["sed_flux"][row["Observed Passband"]+" "+str(band_itr)] = row["Flux Density"]
            sed["sed_freq"][row["Observed Passband"]+" "+str(band_itr)] = row["Frequency"]
            sed["sed_unc_upper"][row["Observed Passband"]+" "+str(band_itr)] = row["Upper limit of uncertainty"]
            sed["sed_unc_lower"][row["Observed Passband"]+" "+str(band_itr)] = row["Lower limit of uncertainty"]
            sed["sed_telescopenames"][row["Observed Passband"]+" "+str(band_itr)] = "IRAS"
            sed["sed_flags"][row["Observed Passband"]+" "+str(band_itr)] = "n"
            band_itr += 1
            newrow = True
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
            newrow = True
        # Spitzer MIPS
        if "MIPS" in row["Observed Passband"] and\
                "Total" in row["Spatial Mode"] and\
                row["Frequency"] < 1.3e13:
            sed["sed_flux"][row["Observed Passband"]+" "+str(band_itr)] = row["Flux Density"]
            sed["sed_freq"][row["Observed Passband"]+" "+str(band_itr)] = row["Frequency"]
            sed["sed_unc_upper"][row["Observed Passband"]+" "+str(band_itr)] = row["Upper limit of uncertainty"]
            sed["sed_unc_lower"][row["Observed Passband"]+" "+str(band_itr)] = row["Lower limit of uncertainty"]
            sed["sed_telescopenames"][row["Observed Passband"]+" "+str(band_itr)] = "Spitzer"
            sed["sed_flags"][row["Observed Passband"]+" "+str(band_itr)] = "n"
            band_itr += 1
            newrow = True
        # Parkes
        if row["Refcode"] == "1990PKS90.C...0000W"\
                or row["Refcode"] == "1977MNRAS.179..235D":
            sed["sed_flux"][row["Observed Passband"]+" "+str(band_itr)] = row["Flux Density"]
            sed["sed_freq"][row["Observed Passband"]+" "+str(band_itr)] = row["Frequency"]
            sed["sed_unc_upper"][row["Observed Passband"]+" "+str(band_itr)] = row["Upper limit of uncertainty"]
            sed["sed_unc_lower"][row["Observed Passband"]+" "+str(band_itr)] = row["Lower limit of uncertainty"]
            sed["sed_telescopenames"][row["Observed Passband"]+" "+str(band_itr)] = "Parkes"
            sed["sed_flags"][row["Observed Passband"]+" "+str(band_itr)] = "n"
            band_itr += 1
            newrow = True
        # ATCA
        if row["Refcode"] == "2010MNRAS.402.2403M"\
                or row["Refcode"] == "2005ApJ...623..815G":
            sed["sed_flux"][row["Observed Passband"]+" "+str(band_itr)] = row["Flux Density"]
            sed["sed_freq"][row["Observed Passband"]+" "+str(band_itr)] = row["Frequency"]
            sed["sed_unc_upper"][row["Observed Passband"]+" "+str(band_itr)] = row["Upper limit of uncertainty"]
            sed["sed_unc_lower"][row["Observed Passband"]+" "+str(band_itr)] = row["Lower limit of uncertainty"]
            sed["sed_telescopenames"][row["Observed Passband"]+" "+str(band_itr)] = "ATCA"
            sed["sed_flags"][row["Observed Passband"]+" "+str(band_itr)] = "n"
            band_itr += 1
            newrow = True
        # VLA
        if (row["Refcode"] == "2002AJ....124..675C"\
                or row["Refcode"] == "2006A&A...449..559B"\
                or row["Refcode"] == "2004AJ....128.2013C"\
                or row["Refcode"] == "2005A&A...435..521N"\
                or row["Refcode"] == "2009ApJ...703..802M"\
                or row["Refcode"] == "2017ApJ...847..136B"\
                or row["Refcode"] == "2006A&A...455..161L"\
                or row["Refcode"] == "1998AJ....115.1693C"\
                or row["Refcode"] == "2011ApJ...731L..41B")\
                and row["Frequency"] < 1e12:
            sed["sed_flux"][row["Observed Passband"]+" "+str(band_itr)] = row["Flux Density"]
            sed["sed_freq"][row["Observed Passband"]+" "+str(band_itr)] = row["Frequency"]
            sed["sed_unc_upper"][row["Observed Passband"]+" "+str(band_itr)] = row["Upper limit of uncertainty"]
            sed["sed_unc_lower"][row["Observed Passband"]+" "+str(band_itr)] = row["Lower limit of uncertainty"]
            sed["sed_telescopenames"][row["Observed Passband"]+" "+str(band_itr)] = "VLA"
            sed["sed_flags"][row["Observed Passband"]+" "+str(band_itr)] = "n"
            band_itr += 1
            newrow = True
        # Effelsberg
        if row["Refcode"] == "2004A&A...418....1V":
            sed["sed_flux"][row["Observed Passband"]+" "+str(band_itr)] = row["Flux Density"]
            sed["sed_freq"][row["Observed Passband"]+" "+str(band_itr)] = row["Frequency"]
            sed["sed_unc_upper"][row["Observed Passband"]+" "+str(band_itr)] = row["Upper limit of uncertainty"]
            sed["sed_unc_lower"][row["Observed Passband"]+" "+str(band_itr)] = row["Lower limit of uncertainty"]
            sed["sed_telescopenames"][row["Observed Passband"]+" "+str(band_itr)] = "Effelsberg"
            sed["sed_flags"][row["Observed Passband"]+" "+str(band_itr)] = "n"
            band_itr += 1
            newrow = True
        # MOST
        if row["Refcode"] == "2003MNRAS.342.1117M"\
                or row["Refcode"] == "2008SUMSS.2.1.....:":
            sed["sed_flux"][row["Observed Passband"]+" "+str(band_itr)] = row["Flux Density"]
            sed["sed_freq"][row["Observed Passband"]+" "+str(band_itr)] = row["Frequency"]
            sed["sed_unc_upper"][row["Observed Passband"]+" "+str(band_itr)] = row["Upper limit of uncertainty"]
            sed["sed_unc_lower"][row["Observed Passband"]+" "+str(band_itr)] = row["Lower limit of uncertainty"]
            sed["sed_telescopenames"][row["Observed Passband"]+" "+str(band_itr)] = "MOST"
            sed["sed_flags"][row["Observed Passband"]+" "+str(band_itr)] = "n"
            band_itr += 1
            newrow = True
        # Arecibo
        if row["Refcode"] == "1978ApJS...36...53D":
            sed["sed_flux"][row["Observed Passband"]+" "+str(band_itr)] = row["Flux Density"]
            sed["sed_freq"][row["Observed Passband"]+" "+str(band_itr)] = row["Frequency"]
            sed["sed_unc_upper"][row["Observed Passband"]+" "+str(band_itr)] = row["Upper limit of uncertainty"]
            sed["sed_unc_lower"][row["Observed Passband"]+" "+str(band_itr)] = row["Lower limit of uncertainty"]
            sed["sed_telescopenames"][row["Observed Passband"]+" "+str(band_itr)] = "Arecibo"
            sed["sed_flags"][row["Observed Passband"]+" "+str(band_itr)] = "n"
            band_itr += 1
            newrow = True

        # If upper limit
        if row["Flux Density"] is np.ma.masked and row["Upper limit of Flux Density"] > 0.0 and newrow:
            sed["sed_unc_upper"][row["Observed Passband"]+" "+str(band_itr-1)] = -1
            sed["sed_flux"][row["Observed Passband"]+" "+str(band_itr-1)] = row["Upper limit of Flux Density"]

        # Deal with potentially empty boxes
        for key in sed:
            for fltr in sed[key]:
                if sed[key][fltr] is np.ma.masked:
                    sed[key][fltr] = 0.0

        # Deal with absence of upper limits
        for key in sed["sed_unc_upper"]:
            if sed["sed_unc_upper"][key] == 0.0:
                sed["sed_unc_upper"][key] = 0.1 * sed["sed_flux"][key]
                sed["sed_unc_lower"][key] = 0.1 * sed["sed_flux"][key]

    return sed
