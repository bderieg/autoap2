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
