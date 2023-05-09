import numpy as np
import astropy.units as u


# Modified blackbody function

def mb(params, nu, data, unc_upper, unc_lower):
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

    lower_weights = [ 1/(np.sqrt(s)) for s in unc_lower ]
    upper_weights = [ 1/(np.sqrt(s)) for s in unc_upper ]
    lower_weights = [ -1 if np.isnan(w) else w for w in lower_weights ]
    upper_weights = [ -1 if np.isnan(w) else w for w in upper_weights ]
    lower_weights = [ 0 if np.isinf(w) else w for w in lower_weights ]
    upper_weights = [ 0 if np.isinf(w) else w for w in upper_weights ]
    lower_weights = [ max(lower_weights) if w==0 else w for w in lower_weights ]
    upper_weights = [ max(upper_weights) if w==0 else w for w in upper_weights ]
    lower_weights = [ lw if uw>0 else 0 for lw,uw in zip(lower_weights,upper_weights) ]
    upper_weights = [ uw if uw>0 else max(upper_weights) for uw in upper_weights ]

    return [ lw*(m-d)**2 if m<d else uw*(m-d)**2 for m,d,lw,uw in zip(model,data,lower_weights,upper_weights) ]

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

def powerlaw(params, nu, data, unc_upper, unc_lower):
    alpha = params['alpha'].value
    b = params['b'].value

    model = [ b*n**alpha for n in nu ]

    lower_weights = [ 1/(np.sqrt(s)) for s in unc_lower ]
    upper_weights = [ 1/(np.sqrt(s)) for s in unc_upper ]
    lower_weights = [ -1 if np.isnan(w) else w for w in lower_weights ]
    upper_weights = [ -1 if np.isnan(w) else w for w in upper_weights ]
    lower_weights = [ 0 if np.isinf(w) else w for w in lower_weights ]
    upper_weights = [ 0 if np.isinf(w) else w for w in upper_weights ]
    lower_weights = [ max(lower_weights) if w==0 else w for w in lower_weights ]
    upper_weights = [ max(upper_weights) if w==0 else w for w in upper_weights ]
    lower_weights = [ lw if uw>0 else 0 for lw,uw in zip(lower_weights,upper_weights) ]
    upper_weights = [ uw if uw>0 else max(upper_weights) for uw in upper_weights ]

    return [ lw*(m-d)**2 if m<d else uw*(m-d)**2 for m,d,lw,uw in zip(model,data,lower_weights,upper_weights) ]

def powerlaw_basic(params, nu):
    alpha = params['alpha']
    b = params['b']

    model = [ b*n**alpha for n in nu ]

    return model
