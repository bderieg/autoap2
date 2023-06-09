import numpy as np
import astropy.units as u

def mb_model(p, x=None):

    nu = x  #in Ghz
    Md = p[0]  # in Solar Masses
    Td = p[1]  # in Kelvin
    beta = p[2]  # dust emissivity value
    Da = p[3]  # in Mpc

    model = 5.913E-13*nu**3*(nu/856.6)**beta*Md/Da**2*(np.exp(0.048*(nu/Td))-1)**(-1)

    return model

def mb_fit(p, fjac=None, x=None, y=None, err=None):

    model = mb_model(p, x)
    status = 0

    return [status, (y-model)/err]

def pl_model(nu, a, s0):
    return s0 * nu**a
