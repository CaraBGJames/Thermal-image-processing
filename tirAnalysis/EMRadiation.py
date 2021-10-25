#!/usr/bin/python

""" EMRadiation.py
    A set of routines that are useful for calculations involving electro-
    magnetic (EM) radiation.

    This package provides the following routines:
    wavelength2wavenumber
    wavenumber2wavelength
    planckFunction
    rayleighJeans
    wienApprox
"""

from numpy import pi, exp, divide

import scipy.constants as const
import numpy as np

# A list of physical constants.  All units in MKS-SI
kB    = const.k     # Boltzmann constant/[J/K]
sigSB = const.sigma # Stefan-Boltzmann constant/[W/m**2/K**4]
c     = const.c     # Speed of light in vacuuo/[m/s]
h     = const.h     # Planck constant/[J.s]
b     = 2.8977729e-3  # constant from Wein's displacement law/[m.K]


def wavelength2wavenumber(lam, inputunits='micron', outputunits='percm',
                          rev=False):
    """ Converts wavelength to wavenumber according to 
        nu/[cm^-1] = 10000 / lam/[micron] """
    if inputunits == 'micron':
        # Use of numpy.divide allows lists to be dealt with
        nu = divide(10000, lam)
    else:
        raise TypeError('Unit of measure not currently supported.  ' + 
              'Convert wavelength to microns and try again.  ')
    # Express units in cm**-1 if required
    # if requnits == 'percm':
    #     nu *= 100

    # Reverse the list if required
    if rev:
        nu = nu[::-1]  
    return nu


def wavenumber2wavelength(nu, inputunits='percm'):
    """ Converts wavenumber to wavelength according to lam = 1. / nu """
    num = 1.
    if inputunits == 'percm':
        num /= 100.
    if inputunits == 'permm':
        num /= 1000.
    
    return np.divide(num, nu)


def planckFunction(lam, T):
    """ Returns the Planck function at wavelength lam and temperature T.  
    Currently accepts only SI units """
    denom = exp(h * c / (lam * kB * T)) - 1.
    return 2. * h * c**2 / (lam**5 * denom)


def planckFunctionCGS(lam, T):
    """ Returns the Planck function at wavelength lam and temperature T in CGS 
    units.  
    
    Parameters:
    -----------
    lam         wavelength/[cm]
    T           temperature/[K]

    
    """
    # Constant definitions (in CGS base)
    h  = 6.62606957e-27                # Planck's constant/[ergs]
    c  = 2.99792458e+10                # Speed of light in vaccuo/[cm/s]
    k  = 1.3806488e-16                 # Boltzmann's constant/[erg/K]
    denom = exp(h * c / (lam * k * T)) - 1.
    return divide(2. * h * c**2, lam**5 * denom)


def planckFuncWavenumber(nu, T):
    """ Returns the Planck function at wavenumber nu and temperature T.  Note
    that the wavenumber is expressed in cm**-1 so that the output is not in 
    true MKS.
    """
    return 2. * h * c**2 * nu**3 / (exp(h * c * nu / (kB * T)) - 1.)


def rayleighJeans(lam, T):
    return 2. * c * kB * T / lam**4


def weinApprox(lam, T):
    return 2. * h * c**2 / lam**5 * exp(- h * c / (lam * kB * T))


def weinDisplacementLaw(T, units='metres'):
    """ Returns the wavelength corresponding to the peak of the Planck function
    given by Wein's Displacement Law, \lambda_{max} = b / T """
    lamMax = b / T
    if units == 'microns':
        lamMax *= 1e6
    return lamMax


def stefanBoltzmann(T, emissivity=1.):
    """ Returns the radiative flux for a body at absolute temperature, T, 
    according to the Stefan-Boltzmann law.  By default, the wavelength-averaged 
    emissivity is taken to be unity.
    """
    return emissivity * sigSB * T**4

