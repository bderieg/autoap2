import pandas as pd
import re

# properties = pd.read_csv('/home/ben/Desktop/research/research_boizelle_working/ap_phot_data/galaxy_properties.csv', index_col=0)
properties = pd.read_excel('/home/ben/Desktop/research/research_boizelle_working/kinemetry_working/kinemetry_progress.ods', engine='odf', sheet_name='Target Parameters', index_col=0, usecols='A:V')

outputfile_loc = "/home/ben/Desktop/sample_params.tex"
with open(outputfile_loc,"w") as f:
    for itr, row in properties.iterrows():
        numrows = 1
        try:
            while str(properties.index[properties.index.get_loc(itr)+numrows]) == "nan":
                numrows += 1
        except IndexError:
            pass
        if str(itr) == "nan":
            target_ms = " & " +\
                        "& " +\
                        "& " +\
                        "& " +\
                        "& " +\
                        "& " +\
                        "\\\\"
        elif numrows == 1:
            target_ms = str(itr) + " & " +\
                        str(row['RC3 type']) + " & " +\
                        str(row['luminosity distance (Mpc)']) + " (" + str(row['luminosity distance unc. (Mpc)']) + ") & " +\
                        str(row['dust minor axis (arcsec)']) + "/" + str(row['dust major axis (arcsec)']) + " & " +\
                        str(row['gas minor axis (arcsec)']) + "/" + str(row['gas major axis (arcsec)']) + " & " +\
                        str(round(row['gas sma / beam sma'], 1)) + " & " +\
                        str(int(row['distance reference'])) + " \\\\"
        else:
            target_ms = "\multirow{" + str(numrows) + "}{*}{" + str(itr) + "} & " +\
                        "\multirow{" + str(numrows) + "}{*}{" + str(row['RC3 type']) + "} & " +\
                        "\multirow{" + str(numrows) + "}{*}{" + str(row['luminosity distance (Mpc)']) + " (" + str(row['luminosity distance unc. (Mpc)']) + ")} & " +\
                        "\multirow{" + str(numrows) + "}{*}{" + str(row['dust minor axis (arcsec)']) + "/" + str(row['dust major axis (arcsec)']) + "} & " +\
                        "\multirow{" + str(numrows) + "}{*}{" + str(row['gas minor axis (arcsec)']) + "/" + str(row['gas major axis (arcsec)']) + "} & " +\
                        "\multirow{" + str(numrows) + "}{*}{" + str(round(row['gas sma / beam sma'],1)) + "} & " +\
                        "\multirow{" + str(numrows) + "}{*}{" + str(int(row['distance reference'])) + "} \\\\"

        target_ms = target_ms.replace("nan", "\\nodata")
        target_ms = target_ms.replace(r"\nodata/\nodata", r"\nodata")
        target_ms = target_ms.replace(r"\nodata (\nodata)", r"\nodata")
        target_ms = target_ms.replace(r"$\nodata$", r"\nodata")

        print(target_ms, file=f)
