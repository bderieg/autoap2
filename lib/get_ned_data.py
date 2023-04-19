from astroquery.ipac.ned import Ned as ned
from statistics import mean

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
            "sed_data" : {},
            "sed_unc_upper" : {},
            "sed_unc_lower" : {},
            "sed_telescopenames" : {},
            "sed_filternames" : {},
            "sed_flags" : {}
        }

    iras_itr = 0
    iso_itr = 0
    mips_itr = 0
    for row in sed_table:
        if row["Flux Density"] > 0.0:
            # Increment frequency if already there
            tempfreq = row["Frequency"]
            while str(tempfreq) in sed["sed_data"]:
                tempfreq += 1.0
            # IRAS
            if "quality flag" in row["Qualifiers"] and\
                    (("100" in row["Observed Passband"]) or ("60" in row["Observed Passband"])) and\
                    "IRAS" in row["Observed Passband"]:
                sed["sed_data"][str(tempfreq)] = row["Flux Density"]
                sed["sed_unc_upper"][str(tempfreq)] = row["Upper limit of uncertainty"]
                sed["sed_unc_lower"][str(tempfreq)] = row["Lower limit of uncertainty"]
                sed["sed_telescopenames"][str(tempfreq)] = "IRAS"
                sed["sed_filternames"][str(tempfreq)] = row["Observed Passband"] + " " + str(iras_itr)
                sed["sed_flags"][str(row["Observed Passband"]) + " " + str(iras_itr)] = "n"
                iras_itr += 1
            # ISO
            if "ISO" in row["Observed Passband"] and\
                    "[" not in row["Observed Passband"] and\
                    row["Frequency"] < 1e13:
                sed["sed_data"][str(tempfreq)] = row["Flux Density"]
                sed["sed_unc_upper"][str(tempfreq)] = row["Upper limit of uncertainty"]
                sed["sed_unc_lower"][str(tempfreq)] = row["Lower limit of uncertainty"]
                sed["sed_telescopenames"][str(tempfreq)] = "ISO"
                sed["sed_filternames"][str(tempfreq)] = row["Observed Passband"] + " " + str(iso_itr)
                sed["sed_flags"][str(row["Observed Passband"]) + " " + str(iso_itr)] = "n"
                iso_itr += 1
            # Spitzer MIPS
            if "MIPS" in row["Observed Passband"] and\
                    "Total" in row["Spatial Mode"] and\
                    row["Frequency"] < 1e13:
                sed["sed_data"][str(tempfreq)] = row["Flux Density"]
                sed["sed_unc_upper"][str(tempfreq)] = row["Upper limit of uncertainty"]
                sed["sed_unc_lower"][str(tempfreq)] = row["Lower limit of uncertainty"]
                sed["sed_telescopenames"][str(tempfreq)] = "Spitzer"
                sed["sed_filternames"][str(tempfreq)] = row["Observed Passband"] + " " + str(mips_itr)
                sed["sed_flags"][str(row["Observed Passband"]) + " " + str(mips_itr)] = "n"
                mips_itr += 1

    return sed
