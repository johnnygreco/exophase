import numpy as np

def ang_sep(orbdist, TA, argperi, incl, D, input_deg=False):
    """
    Angular separation in arcsec
    orbitdist must be in AU, and D must be in pc
    """
    f = TA + argperi
    if input_deg:
        f *= np.pi/180.0
    x = np.cos(f)
    y = np.sin(f)*np.cos(incl)
    sep = orbdist*np.sqrt(x**2 + y**2)/D
    return sep

def orb_dist(TA, a, ecc, input_deg=False):
    """
    Calculate orbital distance given true anamoly (TA)
    semi-major axis (a) and eccentricity (ecc).
    """
    if input_deg:
        TA *= np.pi/180.0
    return a*(1.0 - ecc*ecc)/(1.0 + ecc*np.cos(TA))

def flux_ratio(phi, dist, Ag, Rp, scale=1.0):
    """
    Planet/star flux ratio. 

    Parameters
    ----------
    phi : float or ndarray, phase function
    dist : float or ndarray, orbital distance in AU
    Ag : float, geometric albedo
    Rp : float, planet radius in R_Jup
    scale : float, scale factor (optional)
     
    Returns
    -------
    fr : float or ndarray, flux ratio
    """
    dist *= 2092.51204 # convert orbital distance to R_Jup
    return Ag*np.power(Rp/dist, 2)*phi*scale

def lambert(alpha, input_deg=False):
    """
    Analytic Lambert phase function
    by default, alpha is in radians
    """
    if input_deg:
        alpha = alpha*(np.pi/180.0)
    phi = (np.sin(alpha) + (np.pi - alpha)*np.cos(alpha))/np.pi
    return phi
