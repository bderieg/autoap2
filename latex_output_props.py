import pandas as pd

targets = {
            "GAMA 272990",
            "GAMA 622429",
            "GAMA 064646",
            "IC 1024",
            "NGC 0524",
            "NGC 0708",
            "NGC 0997",
            "NGC 1052",
            "NGC 1332",
            "NGC 1380",
            "NGC 1684",
            "NGC 2872",
            "NGC 3078",
            "NGC 3169",
            "NGC 3258",
            "NGC 3268",
            "NGC 3557",
            "NGC 3599",
            "NGC 3607",
            "NGC 3862",
            "NGC 4061",
            "NGC 4435",
            "NGC 4438",
            "NGC 4477",
            "NGC 4697",
            "NGC 4786",
            "NGC 4876",
            "NGC 5084",
            "NGC 5208",
            "NGC 5838",
            "NGC 6861",
            "NGC 7465",
            "PGC 43387"
        }

properties = pd.read_csv('/home/ben/Desktop/research/research_boizelle_working/ap_phot_data/galaxy_properties.csv', index_col=0)

outputfile_loc = "/home/ben/Desktop/sample_params.tex"
with open(outputfile_loc,"w") as f:
    for target in sorted(targets):
        try:
            target_ms = str(target) + " & " +\
                        str(properties.loc[target]['Type']) + " & " +\
                        str(properties.loc[target]['AGN Type']) + " & " +\
                        "$" + str(properties.loc[target]['Redshift (via NED)']) + "$ & " +\
                        "$" + str(properties.loc[target]['Stellar Velocity Dispersion (km/s)']) + "$ & " +\
                        "$" + str(properties.loc[target]['Angular Scale (pc arcsec^-1)']) + "$ & " +\
                        "$" + str(properties.loc[target]['b/a dust (arcsec)']) + "$ \\\\"

            target_ms = target_ms.replace("nan", "\\nodata")
            target_ms = target_ms.replace("No HST data", "\\nodata")
            target_ms = target_ms.replace("Irreg. dust", "\\nodata")
            target_ms = target_ms.replace("$\\nodata$", "\\nodata")

            print(target_ms, file=f)
        except KeyError:
            print("warning . . . " + target + " not found")
