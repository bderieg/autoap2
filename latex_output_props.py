import pandas as pd

properties = pd.read_excel('/home/ben/Desktop/research/research_boizelle_working/kinemetry_working/kinemetry_progress.ods', engine='odf', sheet_name='Target Parameters', usecols='A:AB')

outputfile_loc = "/home/ben/Desktop/sample_params.tex"
with open(outputfile_loc,"w") as f:
    targets = []
    for itr, row in properties.iterrows():

        if row['target'] not in targets:
            targets.append(row['target'])

            target_ms = str(row['target']) + " & " +\
                        str(row['RC3 type']) + " & " +\
                        str(row['luminosity distance (Mpc)']) + " (" + str(row['luminosity distance unc. (Mpc)']) + ")\\tablenotemark{" + str(row['distance reference']) + "} & " +\
                        str(round(row['stellar velocity dispersion (km/s)'],1)) + " & " +\
                        str(round(row['SMBH mass (solar masses)']/1e7,1))

            target_ms = target_ms.replace("nan", "\\nodata")
            target_ms = target_ms.replace(r"\nodata/\nodata", r"\nodata")
            target_ms = target_ms.replace(r"\nodata (\nodata)", r"\nodata")
            target_ms = target_ms.replace(r"$\nodata$", r"\nodata")

            if itr < len(properties)-1:
                target_ms += " \\\\"
                # if str(properties.iloc[itr+1,0]) != "nan":
                #     target_ms += " \\hline"

            print(target_ms, file=f)
