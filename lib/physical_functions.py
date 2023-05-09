import numpy as np
import astropy.units as u


# Modified blackbody function

def mb(params, nu, data):
    # Define constants in SI units
    h = 6.62607015e-34
    c = 299792458
    k = 1.3806503e-23
    kappa500 = 0.051

    # Final conversion factor
    si_to_Jy = 1e26
   
    beta = params['beta'].value
    kappa = [ (kappa500/(6e11**beta))*n**beta for n in nu ]
    mass = params['mass'].value
    temperature = params['temperature'].value
    distance = params['distance'].value
    blackbody = [ (2*h*n**3)/(c**2) / ( np.exp(h*n/k/temperature)-1 ) for n in nu ]

    model = [ si_to_Jy * kap * mass * b / distance**2 for b,kap in zip(blackbody,kappa) ]

    return [ (m-d)**2 for m,d in zip(model,data) ]

def mb_basic(params, nu):
    h = 6.62607015e-34
    c = 299792458
    k = 1.3806503e-23
    kappa500 = 0.051

    si_to_Jy = 1e26
   
    beta = params['beta']
    kappa = [ (kappa500/(6e11**beta))*n**beta for n in nu ]
    mass = params['mass']
    temperature = params['temperature']
    distance = params['distance']
    blackbody = np.zeros_like(nu)
    for ind, n in zip( range(len(blackbody)), nu ):
        try:
            blackbody[ind] = (2*h*n**3)/(c**2) / ( np.exp(h*n/k/temperature)-1 )
        except RuntimeWarning:
            blackbody[ind] = 0.0

    model = [ si_to_Jy * kap * mass * b / distance**2 for b,kap in zip(blackbody,kappa) ]

    return model

# Power law function

def powerlaw(params, nu, data):
    alpha = params['alpha'].value
    b = params['b'].value

    model = [ b*n**alpha for n in nu ]

    return [ (m-d)**2 for m,d in zip(model,data) ]

def powerlaw_basic(params, nu):
    alpha = params['alpha']
    b = params['b']

    model = [ b*n**alpha for n in nu ]

    return model
