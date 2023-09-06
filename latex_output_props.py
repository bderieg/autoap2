import pandas as pd

properties = pd.read_excel('/home/ben/Desktop/research/research_boizelle_working/kinemetry_working/kinemetry_progress.ods', engine='odf', sheet_name='Target Parameters', usecols='A:AB')

outputfile_loc = "/home/ben/Desktop/sample_params.tex"
with open(outputfile_loc,"w") as f:
    for itr, row in properties.iterrows():
        numrows = 1
        try:
            while str(properties.iloc[itr+numrows,0]) == "nan":
                numrows += 1
        except IndexError:
            pass
        if str(row['target']) == "nan":
            target_ms = " & " +\
                        "& " +\
                        "& " +\
                        "& " +\
                        "& " +\
                        str(round(row['ROI / beam radius'],2)) + " & " +\
                        str(row['transition']) + " & " +\
                        str(row['project'])
        elif numrows == 1:
            target_ms = str(row['target']) + " & " +\
                        str(row['RC3 type']) + " & " +\
                        str(row['luminosity distance (Mpc)']) + " (" + str(row['luminosity distance unc. (Mpc)']) + ")\\tablenotemark{" + str(row['distance reference']) + "} & " +\
                        str(round(row['stellar velocity dispersion (km/s)'],1)) + " & " +\
                        str(round(row['SMBH mass (solar masses)']/1e7,1)) + " & " +\
                        str(round(row['ROI / beam radius'],2)) + " & " +\
                        str(row['transition']) + " & " +\
                        str(row['project'])
        else:
            target_ms = "\multirow{" + str(numrows) + "}{*}{" + str(row['target']) + "} & " +\
                        "\multirow{" + str(numrows) + "}{*}{" + str(row['RC3 type']) + "} & " +\
                        "\multirow{" + str(numrows) + "}{*}{" + str(row['luminosity distance (Mpc)']) + " (" + str(row['luminosity distance unc. (Mpc)']) + ")\\tablenotemark{" + str(row['distance reference']) + "}} & " +\
                        "\multirow{" + str(numrows) + "}{*}{" + str(round(row['stellar velocity dispersion (km/s)'],1)) + "} & " +\
                        "\multirow{" + str(numrows) + "}{*}{" + str(round(row['SMBH mass (solar masses)']/1e7,1)) + "} & " +\
                        str(round(row['ROI / beam radius'],2)) + " & " +\
                        str(row['transition']) + " & " +\
                        str(row['project'])

        target_ms = target_ms.replace("nan", "\\nodata")
        target_ms = target_ms.replace(r"\nodata/\nodata", r"\nodata")
        target_ms = target_ms.replace(r"\nodata (\nodata)", r"\nodata")
        target_ms = target_ms.replace(r"$\nodata$", r"\nodata")

        if itr < len(properties)-1:
            target_ms += " \\\\"
            # if str(properties.iloc[itr+1,0]) != "nan":
            #     target_ms += " \\hline"

        print(target_ms, file=f)
