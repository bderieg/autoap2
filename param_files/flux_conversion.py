import numpy as np

abs_unc = {
    'ALMA' : 0.1,
    'IRAC1' : 0.1,
    'IRAC2' : 0.1,
    'W1' : 0.024,
    'W2' : 0.028,
    'W3' : 0.045,
    'W4' : 0.057,
    'PACS70' : 0.07,
    'PACS100' : 0.07,
    'PACS160' : 0.07,
    'SPIRE250' : 0.1,
    'SPIRE350' : 0.1,
    'SPIRE500' : 0.1,
    'SDSS_u' : 0.055,
    'SDSS_g' : 0.048,
    'SDSS_r' : 0.048,
    'SDSS_i' : 0.048,
    'SDSS_z' : 0.050,
    'GALEX_NUV' : 0.103,
    'GALEX_FUV' : 0.110,
    'DES_u' : 0.011,
    'DES_Y' : 0.011,
    'DES_g' : 0.010,
    'DES_r' : 0.010,
    'DES_i' : 0.010,
    'DES_z' : 0.011,
    'DECam_u' : 0.011,
    'DECam_Y' : 0.011,
    'DECam_g' : 0.010,
    'DECam_r' : 0.010,
    'DECam_i' : 0.010,
    'DECam_z' : 0.011,
    '2MASS_J' : 0.017,
    '2MASS_H' : 0.020,
    '2MASS_K' : 0.019,
    '2MASS_J_tile' : 0.017,
    '2MASS_H_tile' : 0.020,
    '2MASS_K_tile' : 0.019
}

color_correction = {
    'ALMA' : lambda header : 1,
    'IRAC1' : lambda header : 1.0111,
    'IRAC2' : lambda header : 1.0121,
    'W1' : lambda header : 1.0038,
    'W2' : lambda header : 1.0512,
    'W3' : lambda header : 1.0030,
    'W4' : lambda header : 1.0013,
    'PACS70' : lambda header : 1,
    'PACS100' : lambda header : 1,
    'PACS160' : lambda header : 1,
    'SPIRE250' : lambda header : 0.9881,
    'SPIRE350' : lambda header : 0.9757,
    'SPIRE500' : lambda header : 0.9853,
    'SDSS_u' : lambda header : 1,
    'SDSS_g' : lambda header : 1,
    'SDSS_r' : lambda header : 1,
    'SDSS_i' : lambda header : 1,
    'SDSS_z' : lambda header : 1,
    'GALEX_NUV' : lambda header : 1,
    'GALEX_FUV' : lambda header : 1,
    'DES_u' : lambda header : 1,
    'DES_Y' : lambda header : 1,
    'DES_g' : lambda header : 1,
    'DES_r' : lambda header : 1,
    'DES_i' : lambda header : 1,
    'DES_z' : lambda header : 1,
    'DECam_u' : lambda header : 1,
    'DECam_Y' : lambda header : 1,
    'DECam_g' : lambda header : 1,
    'DECam_r' : lambda header : 1,
    'DECam_i' : lambda header : 1,
    'DECam_z' : lambda header : 1,
    '2MASS_J' : lambda header : 1,
    '2MASS_H' : lambda header : 1,
    '2MASS_K' : lambda header : 1,
    '2MASS_J_tile' : lambda header : 1,
    '2MASS_H_tile' : lambda header : 1,
    '2MASS_K_tile' : lambda header : 1
}

brightness_conversion = {
    'ALMA' : lambda header : 1,
    'IRAC1' : lambda header : 1e6 / 4.25452e10,  # MJy->Jy * sr->arcsec^2
    'IRAC2' : lambda header : 1e6 / 4.25452e10,  # MJy->Jy * sr->arcsec^2
    'W1' : lambda header : np.exp(-0.917*header['MAGZP']) * 306.682 /\
            color_correction['W1'](header),  # DN -> Jy
    'W2' : lambda header : np.exp(-0.917*header['MAGZP']) * 170.663 /\
            color_correction['W2'](header),  # DN -> Jy
    'W3' : lambda header : np.exp(-0.917*header['MAGZP']) * 29.045 /\
            color_correction['W3'](header),  # DN -> Jy
    'W4' : lambda header : np.exp(-0.917*header['MAGZP']) * 8.284 /\
            color_correction['W4'](header),  # DN -> Jy
    'PACS70' : lambda header : 1,
    'PACS100' : lambda header : 1,
    'PACS160' : lambda header : 1,
    'SPIRE250' : lambda header : 1,
    'SPIRE350' : lambda header : 1,
    'SPIRE500' : lambda header : 1,
    'SDSS_u' : lambda header : 3.631e-6,  # nanomaggy -> Jy
    'SDSS_g' : lambda header : 3.631e-6,  # nanomaggy -> Jy
    'SDSS_r' : lambda header : 3.631e-6,  # nanomaggy -> Jy
    'SDSS_i' : lambda header : 3.631e-6,  # nanomaggy -> Jy
    'SDSS_z' : lambda header : 3.631e-6,  # nanomaggy -> Jy
    'GALEX_NUV' : lambda header : 36e-6,  # count rate -> Jy
    'GALEX_FUV' : lambda header : 108e-6, # count rate -> Jy
    'DES_u' : lambda header : 1e-9 * 3.63 * 10**((30-header['MAGZERO'])/2.5),  # ADU -> Jy
    'DES_Y' : lambda header : 1e-9 * 3.63 * 10**((30-header['MAGZERO'])/2.5),  # ADU -> Jy
    'DES_g' : lambda header : 1e-9 * 3.63 * 10**((30-header['MAGZERO'])/2.5),  # ADU -> Jy
    'DES_r' : lambda header : 1e-9 * 3.63 * 10**((30-header['MAGZERO'])/2.5),  # ADU -> Jy
    'DES_i' : lambda header : 1e-9 * 3.63 * 10**((30-header['MAGZERO'])/2.5),  # ADU -> Jy
    'DES_z' : lambda header : 1e-9 * 3.63 * 10**((30-header['MAGZERO'])/2.5),  # ADU -> Jy
    'DECam_u' : lambda header : 1e-9 * 3.63 * 10**((30-header['MAGZERO'])/2.5),  # ADU -> Jy
    'DECam_Y' : lambda header : 1e-9 * 3.63 * 10**((30-header['MAGZERO'])/2.5),  # ADU -> Jy
    'DECam_g' : lambda header : 1e-9 * 3.63 * 10**((30-header['MAGZERO'])/2.5),  # ADU -> Jy
    'DECam_r' : lambda header : 1e-9 * 3.63 * 10**((30-header['MAGZERO'])/2.5),  # ADU -> Jy
    'DECam_i' : lambda header : 1e-9 * 3.63 * 10**((30-header['MAGZERO'])/2.5),  # ADU -> Jy
    'DECam_z' : lambda header : 1e-9 * 3.63 * 10**((30-header['MAGZERO'])/2.5),  # ADU -> Jy
    '2MASS_J' : lambda header : 1594 * np.exp(0.917*(-header['JMAGZP'])),  # DN -> Jy
    '2MASS_H' : lambda header : 1024 * np.exp(0.917*(-header['HMAGZP'])),  # DN -> Jy
    '2MASS_K' : lambda header : 666.7 * np.exp(0.917*(-header['KMAGZP'])),  # DN -> Jy
    '2MASS_J_tile' : lambda header : 1594 * np.exp(0.917*(-header['MAGZP'])),  # DN -> Jy
    '2MASS_H_tile' : lambda header : 1024 * np.exp(0.917*(-header['MAGZP'])),  # DN -> Jy
    '2MASS_K_tile' : lambda header : 666.7 * np.exp(0.917*(-header['MAGZP']))  # DN -> Jy
}

beam_size = {
    'ALMA' : lambda header : (4*np.log(2)/np.pi) * (header['BMIN']*header['BMAJ']),
    'IRAC1' : lambda header : 1,
    'IRAC2' : lambda header : 1,
    'W1' : lambda header : 1,
    'W2' : lambda header : 1,
    'W3' : lambda header : 1,
    'W4' : lambda header : 1,
    'PACS70' : lambda header : (5.74*6.26)/(header['PIXSIZE']**2) if header['SCANVELO']==20.0 else (8.80*9.60)/(header['PIXSIZE']**2),  # pixels
    'PACS100' : lambda header : (6.98*7.42)/(header['PIXSIZE']**2) if header['SCANVELO']==20.0 else (9.73*10.69)/(header['PIXSIZE']**2),  # pixels
    'PACS160' : lambda header : (10.546*12.27)/(header['PIXSIZE']**2) if header['SCANVELO']==20.0 else (11.51*13.65)/(header['PIXSIZE']**2),  # pixels
    'SPIRE250' : lambda header : 446.75,  # arcsec^2
    'SPIRE350' : lambda header : 786.38,  # arcsec^2
    'SPIRE500' : lambda header : 1630.96,  # arcsec^
    'SDSS_u' : lambda header : 1,
    'SDSS_g' : lambda header : 1,
    'SDSS_r' : lambda header : 1,
    'SDSS_i' : lambda header : 1,
    'SDSS_z' : lambda header : 1,
    'GALEX_NUV' : lambda header : 1,
    'GALEX_FUV' : lambda header : 1,
    'DES_u' : lambda header : 1,
    'DES_Y' : lambda header : 1,
    'DES_g' : lambda header : 1,
    'DES_r' : lambda header : 1,
    'DES_i' : lambda header : 1,
    'DES_z' : lambda header : 1,
    'DECam_u' : lambda header : 1,
    'DECam_Y' : lambda header : 1,
    'DECam_g' : lambda header : 1,
    'DECam_r' : lambda header : 1,
    'DECam_i' : lambda header : 1,
    'DECam_z' : lambda header : 1,
    '2MASS_J' : lambda header : 1,
    '2MASS_H' : lambda header : 1,
    '2MASS_K' : lambda header : 1,
    '2MASS_J_tile' : lambda header : 1,
    '2MASS_H_tile' : lambda header : 1,
    '2MASS_K_tile' : lambda header : 1
}

pix_size = {
    'ALMA' : lambda header : header['CDELT2']**2,
    'IRAC1' : lambda header : header['PXSCAL2']**2,  # arcsec^2
    'IRAC2' : lambda header : header['PXSCAL2']**2,  # arcsec^2
    'W1' : lambda header : 1,
    'W2' : lambda header : 1,
    'W3' : lambda header : 1,
    'W4' : lambda header : 1,
    'PACS70' : lambda header : (5.74*6.26)/(header['PIXSIZE']**2) if header['SCANVELO']==20.0 else (8.80*9.60)/(header['PIXSIZE']**2),  # beam^-1
    'PACS100' : lambda header : (6.98*7.42)/(header['PIXSIZE']**2) if header['SCANVELO']==20.0 else (9.73*10.69)/(header['PIXSIZE']**2),  # beam^-1
    'PACS160' : lambda header : (10.546*12.27)/(header['PIXSIZE']**2) if header['SCANVELO']==20.0 else (11.51*13.65)/(header['PIXSIZE']**2),  # beam^-1
    'SPIRE250' : lambda header : 6**2,  # arcsec^2
    'SPIRE350' : lambda header : 10**2,  # arcsec^2
    'SPIRE500' : lambda header : 14**2,  # arcsec^2
    'SDSS_u' : lambda header : 1,
    'SDSS_g' : lambda header : 1,
    'SDSS_r' : lambda header : 1,
    'SDSS_i' : lambda header : 1,
    'SDSS_z' : lambda header : 1,
    'GALEX_NUV' : lambda header : 1,
    'GALEX_FUV' : lambda header : 1,
    'DES_u' : lambda header : 1,
    'DES_Y' : lambda header : 1,
    'DES_g' : lambda header : 1,
    'DES_r' : lambda header : 1,
    'DES_i' : lambda header : 1,
    'DES_z' : lambda header : 1,
    'DECam_u' : lambda header : 1,
    'DECam_Y' : lambda header : 1,
    'DECam_g' : lambda header : 1,
    'DECam_r' : lambda header : 1,
    'DECam_i' : lambda header : 1,
    'DECam_z' : lambda header : 1,
    '2MASS_J' : lambda header : 1,
    '2MASS_H' : lambda header : 1,
    '2MASS_K' : lambda header : 1,
    '2MASS_J_tile' : lambda header : 1,
    '2MASS_H_tile' : lambda header : 1,
    '2MASS_K_tile' : lambda header : 1
}

other_correction = {
    'ALMA' : lambda header,eff_radius : 1,
    'IRAC1' : lambda header,eff_radius : 0.82 * np.exp(-(eff_radius*3600)**0.370) + 0.910,  # Extended aperture correction
    'IRAC2' : lambda header,eff_radius : 1.16 * np.exp(-(eff_radius*3600)**0.433) + 0.94,  # Extended aperture correction
    'W1' : lambda header,eff_radius : 1,
    'W2' : lambda header,eff_radius : 1,
    'W3' : lambda header,eff_radius : 1,
    'W4' : lambda header,eff_radius : 1,
    'PACS70' : lambda header,eff_radius : 1,
    'PACS100' : lambda header,eff_radius : 1,
    'PACS160' : lambda header,eff_radius : 1,
    'SPIRE250' : lambda header,eff_radius : 0.997,  # Pixelization correction
    'SPIRE350' : lambda header,eff_radius : 0.993,  # Pixelization correction
    'SPIRE500' : lambda header,eff_radius : 0.902,  # Pixelization correction
    'SDSS_u' : lambda header,eff_radius : 1,
    'SDSS_g' : lambda header,eff_radius : 1,
    'SDSS_r' : lambda header,eff_radius : 1,
    'SDSS_i' : lambda header,eff_radius : 1,
    'SDSS_z' : lambda header,eff_radius : 1,
    'GALEX_NUV' : lambda header,eff_radius : 1,
    'GALEX_FUV' : lambda header,eff_radius : 1,
    'DES_u' : lambda header,eff_radius: 1,
    'DES_Y' : lambda header,eff_radius: 1,
    'DES_g' : lambda header,eff_radius: 1,
    'DES_r' : lambda header,eff_radius: 1,
    'DES_i' : lambda header,eff_radius: 1,
    'DES_z' : lambda header,eff_radius: 1,
    'DECam_u' : lambda header,eff_radius: 1,
    'DECam_Y' : lambda header,eff_radius: 1,
    'DECam_g' : lambda header,eff_radius: 1,
    'DECam_r' : lambda header,eff_radius: 1,
    'DECam_i' : lambda header,eff_radius: 1,
    'DECam_z' : lambda header,eff_radius: 1,
    '2MASS_J' : lambda header,eff_radius : 1,
    '2MASS_H' : lambda header,eff_radius : 1,
    '2MASS_K' : lambda header,eff_radius : 1,
    '2MASS_J_tile' : lambda header,eff_radius : 1,
    '2MASS_H_tile' : lambda header,eff_radius : 1,
    '2MASS_K_tile' : lambda header,eff_radius : 1
}
