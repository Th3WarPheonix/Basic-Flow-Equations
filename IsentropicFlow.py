
"""Module containing basic isentropic flow relation equations e.g. stagnation temrperature from Mach number and pressure and vice versa, for all properties"""

# Constants
Ru = 8.314 # J/mol-K Universal Gas Constant
Avo = 6.02214076e23 # mol^−1 Avogadro's Number
Boltz = 1.380649e-23 # J/K Boltzmann Constant
Euler = 2.71828182845904523536 # Euler's number

def p0(M, p, gamma=1.4):
    """Returns stagnation pressure from static pressure"""
    P0P1 = (1 + (gamma-1)/2.0*M**2.0)**(gamma/(gamma-1))
    return P0P1*p

def T0(M, T, gamma=1.4):
    """Returns stagnation temperature from static temperature"""
    T0T1 = (1 + (gamma-1)/2.0*M**2.0)
    return T0T1*T

def density0(M, r, gamma=1.4):
    """Returns stagnation density from static density"""
    rho0rho1 = (1 + (gamma-1)/2.0*M**2.0)**(1.0/(gamma-1))
    return rho0rho1*r

def p(M, p0, gamma=1.4):
    """Returns static pressure from stagnation pressure"""
    P0P1 = (1 + (gamma-1)/2.0*M**2.0)**(gamma/(gamma-1))
    return P0P1**-1*p0

def T(M, T0, gamma=1.4):
    """Returns static temperature from stagnation temperature"""
    T0T1 = (1 + (gamma-1)/2.0*M**2.0)
    return T0T1**-1*T0

def density(M, r0, gamma=1.4):
    """Returns static density from stagnation density"""
    rho0rho1 = (1 + (gamma-1)/2.0*M**2.0)**(1.0/(gamma-1))
    return rho0rho1**-1*r0