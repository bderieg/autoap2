from astroquery.ipac.ned import Ned as ned
from statistics import mean

def get_ellipse_parameters(target_name):
    diameters_table = ned.get_table(target_name, table="diameters")
    a_to_b_ratios = []
    pa = []
    for itr in range(len(diameters_table)):
        if diameters_table[itr]["Axis Ratio Flag"] == "(a/b)":
            a_to_b_ratios.append(diameters_table[itr]["Axis Ratio"])
            pa.append(diameters_table[itr]["Position Angle"])
        elif diameters_table[itr]["Axis Ratio Flag"] == "(b/a)":
            a_to_b_ratios.append(1/diameters_table[itr]["Axis Ratio"])
            pa.append(90-diameters_table[itr]["Position Angle"])
    avg_ratio = mean(a_to_b_ratios)
    avg_pa = mean(pa)
    return {"axis_ratio" : avg_ratio, "position_angle" : avg_pa}


def get_coords(target_name):
    main_table = ned.query_object(target_name)
    return float(main_table["RA"]), float(main_table["DEC"])
