
import numpy as np

def meador_smart():
    """Laminar and turbulent methods"""

    pass

def eckert_method(mache, Tw, Te):
    """Eckert Reference temperature method, calculates the temperature at
    which you should reference properties for calculations of a
    comrpessible boundary layer

    Parameters
    ----------
    mache : mach number at the edge of the boundary layer Tw :
    temperature of the wall Te : temperature of the flow at the edge of
    the boundary layer
    
    Returns
    -------
    Tstar
    """

    Tstar = Te*(1+0.032*mache**2+0.58*(Tw/Te-1))
    return Tstar

def suth_viscosity(temperature, gas='air'):
    """Sutherland's law of viscosity. Is strictly valid only for
    single-component gases at low pressure. Works decent for air"""
    match gas:
        case 'air':
            T0 = 273
            mu0 = 1.716e-5
            Sm = 111
        case 'Ar':
            T0 = 273
            mu0 = 2.125e-5
            Sm = 114
        case 'CO2':
            T0 = 273
            mu0 = 1.370e-5
            Sm = 222
        case 'CO':
            T0 = 273
            mu0 = 1.657-5
            Sm = 136
        case 'N2':
            T0 = 273
            mu0 = 1.663e-5
            Sm = 107
        case 'O2':
            T0 = 273
            mu0 = 1.919e-5
            Sm = 139
        case 'H2':
            T0 = 273
            mu0 = 8.411-5
            Sm = 97
        case 'steam':
            T0 = 350
            mu0 = 1.12e-5
            Sm = 1064    

    mu = mu0*(temperature/T0)**(3/2)*(T0+Sm)/(temperature+Sm)
    return mu

def suth_thermcond(temperature, gas='air'):
    """Sutherland's law of thermal conductivity. Is strictly valid only for
    single-component gases at low pressure. Works decent for air"""
    match gas:
        case 'air':
            T0 = 273
            k0 = .0241
            Sk = 194
        case 'Ar':
            T0 = 273
            k0 = .0163
            Sk = 170
        case 'CO2':
            T0 = 273
            k0 = .0146
            Sk = 1800
        case 'CO':
            T0 = 273
            k0 = .0232
            Sk = 180
        case 'N2':
            T0 = 273
            k0 = .0242
            Sk = 150
        case 'O2':
            T0 = 273
            k0 = .0244
            Sk = 240
        case 'H2':
            T0 = 273
            k0 = .168
            Sk = 120
        case 'steam':
            T0 = 300
            k0 = .0181
            Sk = 2200    

    k = k0*(temperature/T0)**(3/2)*(T0+Sk)/(temperature+Sk)
    return k

def viscosity_water(temp):
    """White Viscous Fluid Flow 2e. Returns the viscosity of water in
    milliPascal-seconds"""

    mu = (20-temp)/(temp+96)*(1.2378-(1.303e-3)*(20-temp)+(3.06e-6)*(20-temp)**2+(2.55e-8)*(20-temp)**3)
    mu = (10**mu)*1.002
    return mu

def bl_compressible(Minf, Tinf, Pinf, Rgas=287.05):
    """Pressure is constant throughout the boundary layer and is equal
    to the free stream"""

    Tstar = eckert_method(Minf, Twall, Tinf)
    Dstar = Pinf/Rgas/Tinf
    Vstar = suth_viscosity(Tstar)
    Restar = Dsar*ue*x/Vstar

    # Laminar compressible solution
    coeff_frc_lam = 0.664/np.sqrt(Restar)

    # Turbulent compressible solution
    coeff_frc_turb = 0.592/Restar**0.2

    shearstress = 1/2*Dstar*ue**2*coeff_frc


import matplotlib.pyplot as plt

t = np.linspace(0, 50, 11)
v = viscosity_water(t)

print(t, v)