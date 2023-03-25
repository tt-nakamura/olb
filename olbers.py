import numpy as np
from scipy.integrate import quad,odeint,cumtrapz
from scipy.constants import G,h,c,k,parsec

def CosmicTime(a, Omega0=1):
    """ age of the universe
    a = scale factor (scalar or 1d-array)
    Omega0 = density parameter of (dark) matter
    return H0*t = (hubble constant)*(age at a)
           H0*t has the same shape as a
    assume flat geometiry of the universe
    so that cosmological const = 1 - Omega0
    """
    t = odeint(lambda _,x:
               np.sqrt(x/(Omega0 + (1-Omega0)*x**3)),
               0, np.r_[0,a])[1:]
    return np.squeeze(t)

def K_param(Omega, H0=70, Mstar=2e30, Rstar=7e8):
    """ n0*sigma*c/H0
    Omega = density parameter of stars
    H0 = Hubble constant / km/s/Mpc
    Mstar = mass of stars / kg
    Rstar = radius of stars / m
    """
    H0 /= 1000*parsec # Hubble constant / 1/s
    return 3/8*Omega*H0*c*Rstar**2/G/Mstar

def q_param(wavelen, Tstar=6e3):
    """ h*c/lambda/k/T
    wavelen = wavelength / m
    Tstar = surface temperature of stars / K
    """
    return h*c/wavelen/k/Tstar

def olbers(a, Omega0, K, q1=0, q2=np.inf, N=128):
    """ Olbers' paradox on darkness of night sky
    a = scale factor (scalar or 1d-array)
    Omega0 = density parameter of (dark matter)
    K = K-parameter = n0*sigma*c/H0
    q1,q2 = q-parameters = hc/lambda/k/T
            of observable wavelength range
    N = number of trapezoids for integration
    return i = I/I* where I,I* = brightness of
               night sky and daytime sky, resp.
           i has the same shape as a
    assume flat geometry of the universe
    """
    def planck(x): return -x**3*np.exp(-x)/np.expm1(-x)
    def hubble(x): return np.sqrt(x/(Omega0 + (1-Omega0)*x**3))

    j0 = quad(planck, q1, q2)[0]

    i = [] # vectorize for a
    for a in np.atleast_1d(a).flat:
        x = np.linspace(0, a, N)[1:]
        tau = -K*cumtrapz(hubble(x)/x**3, x, initial=0)
        tau -= tau[-1]
        y = a/x
        dj = q1*planck(q1*y) if q1 else 0
        if np.isfinite(q2): dj -= q2*planck(q2*y)
        j = cumtrapz(dj*y/x, x, initial=0)
        j += j0 - j[-1]
        di = x*hubble(x)*np.exp(-tau)*j
        i1 = np.trapz(np.r_[0,di], np.r_[0,x])
        i.append(K*i1/a**4/j0)

    return np.squeeze(i)
