
"""Module containing normal shock equations e.g. property ratios across
a shock"""

import numpy as np

def mach2(M1, gamma=1.4):
    """Returns Mach number after a normal shock from incident Mach
    number"""
    top = (gamma-1)*M1**2+2
    bottom = 2*gamma*M1**2-(gamma-1)
    M2 = np.sqrt(top/bottom)
    return M2

def total_pressure_ratio(M1, gamma=1.4):
    """Returns P02/P01 through a normal shock from incident Mach
    number"""
    top1 = (gamma+1)*M1**2
    bottom1 = (gamma-1)*M1**2+2
    top2 = gamma+1
    bottom2 = 2*gamma*M1**2-(gamma-1)
    Pt2Pt1 = (top1/bottom1)**(gamma/(gamma-1))*(top2/bottom2)**(1/(gamma-1))
    return Pt2Pt1

def pressure_ratio(M, gamma=1.4):
    """Returns static pressure ratio (Ps2/Ps1) across a normal shock"""
    Ps2Ps1 = (2*gamma*M**2-(gamma-1))/(gamma+1)
    return Ps2Ps1

def temperature_ratio(M, gamma=1.4):
    """Returns static temperature ratio (T2/T1) across a normal shock"""
    top1 = 2*gamma*M**2-(gamma-1)
    top2 = (gamma-1)*M**2+2
    bottom = (gamma+1)**2*M**2
    Ts2Ts1 = top1*top2/bottom
    return Ts2Ts1

def density_ratio(M, gamma=1.4):
    """Returns static desnity ratio (rho2/rho1) across a normal shock"""
    top = (gamma+1)*M**2
    bottom = (gamma+1)*M**2
    rho2rho1 = top/bottom
    return rho2rho1

def solve(M1, gamma=1.4):
    """Returns post-shock Mach number and property ratios across the
    shock
    
    Rerturns
    --------
    0: post shock mach
    1 : static pressure ratio Ps2/Ps1
    2 : static temperature ratio Ts2/Ts1
    3 : total pressure ratio Pt2/Pt1
    4 : static density ratio Ds2/Ds1"""
    mach2 = mach2(M1, gamma)
    Pt2Pt1 = total_pressure_ratio(M1, gamma)
    Ps2Ps1 = pressure_ratio(M1, gamma)
    Ts2Ts1 = temperature_ratio(M1, gamma)
    Ds2Ds1 = density_ratio(M1, gamma)

    return (mach2, Ps2Ps1, Ts2Ts1, Pt2Pt1, Ds2Ds1)