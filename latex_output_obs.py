import pandas as pd

properties = pd.read_excel('/home/ben/Desktop/research/research_boizelle_working/kinemetry_working/kinemetry_progress.ods', engine='odf', sheet_name='Target Parameters', usecols='A:AB')

outputfile_loc = "/home/ben/Desktop/obs_params.tex"
with open(outputfile_loc,"w") as f:
    cycle = 0
    for itr, row in properties.iterrows():
        prevcycle = cycle
        cycle = row['cycle']
        try:
            if (prevcycle != cycle) and cycle == 2:
                target_ms = r"\\[-3.6ex] \multicolumn{7}{c}{\textbf{\normalsize Cycle " + str(int(cycle)) + r"}} \\ \hline"
                print(target_ms, file=f)
            elif prevcycle != cycle:
                target_ms = r"\hline \multicolumn{7}{c}{\textbf{\normalsize Cycle " + str(int(cycle)) + r"}} \\ \hline"
                print(target_ms, file=f)
        except ValueError:
            pass
        numrows = 1
        try:
            while str(properties.iloc[itr+numrows,16]) == row['project']:
                numrows += 1
        except IndexError:
            pass
        if str(properties.iloc[itr-1,16]) == row['project']:
            target_ms = "--- & " +\
                        str(row['target']) + " & " +\
                        str(row['transition']) + " & " +\
                        str(round(row['integration time (s)']/60,1)) + " & " +\
                        str(row['velocity bin (km/s)']) + " & " +\
                        str(row['beam major axis (arcsec)']) + "/" + str(row['beam minor axis (arcsec)']) + " & " +\
                        '\\nodata'
        elif numrows == 1:
            target_ms = str(row['project']) + " & " +\
                        str(row['target']) + " & " +\
                        str(row['transition']) + " & " +\
                        str(round(row['integration time (s)']/60,1)) + " & " +\
                        str(row['velocity bin (km/s)']) + " & " +\
                        str(row['beam major axis (arcsec)']) + "/" + str(row['beam minor axis (arcsec)']) + " & " +\
                        '\\nodata'
        else:
            target_ms = str(row['project']) + " & " +\
                        str(row['target']) + " & " +\
                        str(row['transition']) + " & " +\
                        str(round(row['integration time (s)']/60,1)) + " & " +\
                        str(row['velocity bin (km/s)']) + " & " +\
                        str(row['beam major axis (arcsec)']) + "/" + str(row['beam minor axis (arcsec)']) + " & " +\
                        '\\nodata'

        target_ms = target_ms.replace("nan", "\\nodata")
        target_ms = target_ms.replace(r"\nodata/\nodata", r"\nodata")
        target_ms = target_ms.replace(r"\nodata (\nodata)", r"\nodata")
        target_ms = target_ms.replace(r"$\nodata$", r"\nodata")

        if itr < len(properties)-1:
            target_ms += " \\\\"
            # if str(properties.iloc[itr+1,0]) != "nan":
            #     target_ms += " \\hline"

        print(target_ms, file=f)
