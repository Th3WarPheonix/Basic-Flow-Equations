
"""Module containing converging-diverging nozzle equations e.g. property ratios to sonic velocity values, Mach number from area ratio"""

import numpy as np

PRECISION = 1e-7

def SonicTemperatureRatio(M, gamma=1.4):
    '''Returns T*/T from from Mach number'''
    top = 1+(gamma-1)/2*M**2
    bottom = (gamma+1)/2
    TstarT = (top/bottom)
    return TstarT

def SonicPressureRatio(M, gamma=1.4):
    '''Returns P*/P from Mach number'''
    top = 1+(gamma-1)/2*M**2
    bottom = (gamma+1)/2
    PstarP = (top/bottom)**(gamma/(gamma-1))
    return PstarP

def SonicAreaRatio(M, gamma=1.4):
    '''Returns A*/A from Mach number'''
    top = (gamma+1)/2
    bottom = 1+(gamma-1)/2*M**2
    AstarA = M*(top/bottom)**((gamma+1)/(2*(gamma-1)))
    return AstarA

def MachfromAreaRatio(A, gamma=1.4, case=0):
    """Returns Mach number from area ration (A/A*)
    Supersonic solution: case = 0
    Subsonic solution: case = 1
    Source: https://www.grc.nasa.gov/www/winddocs/utilities/b4wind_guide/mach.html"""

    P = 2/(gamma+1)
    Q = 1 - P
    
    if case == 0: # Supersonic case
        R = A**(2*Q/P)
        a = Q**(1/P)

    elif case == 1: # Subsonic case
        R = A**2
        a = P**(1/Q)

    r = (R-1)/(2*a)
    X_new = 1/(1+r+np.sqrt(r*r+2*r))
    X = 0
    
    while abs(X - X_new) > PRECISION: #Newton-raphson
        X = X_new

        if case == 0: # Supersonic case
            F = (P*X+Q)**(1/P) - R*X
            f = (P*X+Q)**((1/P)-1) - R

        elif case == 1: # Subsonic case
            F = (P+Q*X)**(1/Q) - R*X
            f = (P+Q*X)**((1/Q)-1) - R
    
        X_new = X - F/f

    if case == 0:
        M = 1/np.sqrt(X)
    elif case == 1: # Subsonic case
        M = np.sqrt(X)
    
    return M

def find_Astar(Ps, Pt, M, area, gamma=1.4, R_gas=287.05):
    """LOOKIN INTO DEPRECATION IF EQUATION IS NOT VERIFIED
    Find throat area from static and total pressure"""
    Tt = (1+(gamma-1)/2*M**2)**.5
    b = M/Pt
    c = (gamma+1)/(gamma-1)
    d = np.sqrt(gamma/R_gas*(2/(gamma+1))**c)**-1

    Astar = Ps/R_gas*area*np.sqrt(gamma*R_gas)*Tt*b*d

    return Astar

def massflow(Tt, Pt, area, gamma=1.4, R_gas=287.05):
    """Returns the most amount of mass that can pass through a converging diverging nozzle"""
    massflow = area*Pt/np.sqrt(gamma/R_gas/Tt)*((gamma+1)/2)**(-(gamma+1)/2/(gamma-1))
    return massflow