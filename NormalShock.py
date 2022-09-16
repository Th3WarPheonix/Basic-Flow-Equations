
"""Module containing normal shock equations e.g. property ratios across a shock"""

import math

# Constants
Ru = 8.314 # J/mol-K Universal Gas Constant
Avo = 6.02214076e23 # mol^−1 Avogadro's Number
Boltz = 1.380649e-23 # J/K Boltzmann Constant
Euler = 2.71828182845904523536 # Euler's number

def Mach2(M1, gamma=1.4):
    """Returns Mach number after a normal shock from incident Mach number"""
    top = (gamma-1)*M1**2+2
    bottom = 2*gamma*M1**2-(gamma-1)
    M2 = math.sqrt(top/bottom)
    return M2

def stag_pressure_ratio(M1, gamma=1.4):
    """Returns P02/P01 through a normal shock from incident Mach number"""
    top1 = (gamma+1)*M1**2
    bottom1 = (gamma-1)*M1**2+2
    top2 = gamma+1
    bottom2 = 2*gamma*M1**2-(gamma-1)
    P02P01 = (top1/bottom1)**(gamma/(gamma-1))*(top2/bottom2)**(1/(gamma-1))
    return P02P01

def pressure_ratio(M, gamma=1.4):
    """Returns static pressure ratio (P2/P1) across a normal shock"""
    top = 2*gamma*M**2-(gamma-1)
    bottom = gamma+1
    P2P1 = top/bottom
    return P2P1

def temperature_ratio(M, gamma=1.4):
    """Returns static temperature ratio (T2/T1) across a normal shock"""
    top1 = 2*gamma*M**2-(gamma-1)
    top2 = (gamma-1)*M**2+2
    bottom = (gamma+1)**2*M**2
    T2T1 = top1*top2/bottom
    return T2T1

def density_ratio(M, gamma=1.4):
    """Returns static desnity ratio (rho2/rho1) across a normal shock"""
    top = (gamma+1)*M**2
    bottom = (gamma+1)*M**2
    rho2rho1 = top/bottom
    return rho2rho1

def solve_Nshock(M1, gamma=1.4):
    """Returns dictionary of ratios of perintent values across the normal shock from incident Mach number

    Dictionary Keys
    ---------------
    'M2'\n
    'P1P01'\n
    'T1T01'\n
    'rho1rho01'\n
    'P2P02'\n
    'T2T02'\n
    'rho2rho02'\n
    'P02P01'\n
    'P2P1'\n
    'T2T1'\n
    'rho2rho1'\n
    """
    M2 = Mach2(M1, gamma)

    P02P01 = stag_pressure_ratio(M1, gamma)
    P2P1 = pressure_ratio(M1, gamma)
    T2T1 = temperature_ratio(M1, gamma)
    rho2rho1 = density_ratio(M1, gamma)

    result = {}
    result['M2'] = M2
    result['P02P01'] = P02P01
    result['P2P1'] = P2P1
    result['T2T1'] = T2T1
    result['rho2rho1'] = rho2rho1

    return result

def print_solve_Nshock(result, decimals=4):
    """Prints out ouput of solve_Nshock in a structured manner using SI units"""

    print('Mach number after shock', round(result['M2'], decimals))
    print('Stagnation pressure ratio (P02/P01)', round(result['P02P01'], decimals))
    print('Pressure ratio (P2/P1)', round(result['P2P1'], decimals))
    print('Temperature ratio (T2/T1)', round(result['T2T1'], decimals))
    print('Desnity ratio (rho2/rho1)', round(result['rho2rho1'], decimals))