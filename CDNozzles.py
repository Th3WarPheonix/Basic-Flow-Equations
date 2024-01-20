
"""Module containing converging-diverging nozzle equations e.g. property
ratios to sonic velocity values, Mach number from area ratio"""

import numpy as np
import NormalShock as nsh
import IsentropicFlow as isenf

def SonicTemperatureRatio(M, gamma=1.4):
    """Returns T*/Tt from from incident Mach number"""
    top = 1+(gamma-1)/2*M**2
    bottom = (gamma+1)/2
    TstarT = (top/bottom)
    return TstarT

def SonicPressureRatio(M, gamma=1.4):
    """Returns P*/Pt from incident Mach number"""
    top = 1+(gamma-1)/2*M**2
    bottom = (gamma+1)/2
    PstarP = (top/bottom)**(gamma/(gamma-1))
    return PstarP

def SonicAreaRatio(M, gamma=1.4):
    """eturns A*/A from incident Mach number.
    If the contraction is less than the returned value flow is not
    choked"""
    top = (gamma+1)/2
    bottom = 1+(gamma-1)/2*M**2
    AstarA = M*(top/bottom)**((gamma+1)/(2*(gamma-1)))
    return AstarA

def MachfromAreaRatio(AR, gamma=1.4, case=0, precision=1e-7):
    """Returns Mach number from area ratio (A/A*)
    Supersonic solution: case = 0
    Subsonic solution: case = 1
    Source:
    https://www.grc.nasa.gov/www/winddocs/utilities/b4wind_guide/mach.html"""

    P = 2/(gamma+1)
    Q = 1 - P
    if case == 0: # Supersonic case
        R = AR**(2*Q/P)
        a = Q**(1/P)

        r = (R-1)/(2*a)
        X_new = 1/(1+r+np.sqrt(r*r+2*r))
        X = 0
        while abs(X - X_new) > precision: #Newton-raphson
            X = X_new

            F = (P*X+Q)**(1/P) - R*X
            Fprime = (P*X+Q)**((1/P)-1) - R
        
            X_new = X - F/Fprime

        M = 1/np.sqrt(X)

    elif case == 1: # Supersonic case
        R = AR**2
        a = P**(1/Q)

        r = (R-1)/(2*a)
        X_new = 1/(1+r+np.sqrt(r*r+2*r))
        X = 0
        while abs(X - X_new) > precision: #Newton-raphson
            X = X_new

            F = (P+Q*X)**(1/Q) - R*X
            Fprime = (P+Q*X)**((1/Q)-1) - R
        
            X_new = X - F/Fprime

        M = np.sqrt(X)
    
    return M

def find_Astar(Ps, Pt, M, area, gamma=1.4, R_gas=287.05):
    """LOOK INTO DEPRECATION IF EQUATION IS NOT VERIFIED
    Find throat area from static and total pressure"""
    Tt = (1+(gamma-1)/2*M**2)**.5
    b = M/Pt
    c = (gamma+1)/(gamma-1)
    d = np.sqrt(gamma/R_gas*(2/(gamma+1))**c)**-1

    Astar = Ps/R_gas*area*np.sqrt(gamma*R_gas)*Tt*b*d

    return Astar

def massflow(Tt, Pt, area, gamma=1.4, R_gas=287.05):
    """Returns the most amount of mass that can pass through a
    converging diverging nozzle. This assumes the flow is choked"""
    massflow = area*Pt/np.sqrt(gamma/R_gas/Tt)*((gamma+1)/2)**(-(gamma+1)/2/(gamma-1))
    return massflow

def backpressure(Tt, massflow, area, mach=1, gamma=1.4, R_gas=287.05):
    """Returns the total pressure created from a certain amount of mass
    flow through an orifice."""
    Pt = massflow/area*np.sqrt(gamma/R_gas/Tt)*mach*((gamma*1)/2*mach**2)**(-(gamma+1)/2/(gamma-1))
    return Pt

def WTAP(mach, gamma=1.4, R_gas=287.05):
    """Returns the mass flow parameter as a function of mach number.
    Mass flow parameter WTAP = W*sqrt(Tt)/(Pt*area). To get the mass
    flow parameter that chokes the flow set mach number  equal to 1"""
    wtap = np.sqrt(gamma/R_gas)*mach*((gamma+1)/2*mach**2)**(-(gamma+1)/2/(gamma-1))
    return wtap

def shock_location(AR, PbPt1, gamma=1.4, precision=1e-6):
    """Assuming choked at throat"""
    
    ARns1 = (AR-1)/2 + 1
    delta = precision/10
    ARns2 = ARns1 + precision*10
    while abs(ARns2-ARns1)>precision:
        ARns1 = ARns2
        y1 = (_shock_loc(AR, ARns1+delta, gamma)-_shock_loc(AR, ARns1, gamma))/delta
        dx = (PbPt1-_shock_loc(AR, ARns1, gamma))/y1
        ARns2 += dx
    
    return ARns2

def _shock_loc(AR, ARguess, gamma=1.4):
    mach1 = MachfromAreaRatio(ARguess, gamma)
    Ps1Pt1 = isenf.pressure_ratio(mach1, gamma)
    mach2 = nsh.mach2(mach1, gamma)
    Ps2Ps1 = nsh.pressure_ratio(mach1, gamma)
    Ps2Pt2 = isenf.pressure_ratio(mach2, gamma)
    AR2ns = 1/SonicAreaRatio(mach2, gamma)
    AR2e = AR/ARguess*AR2ns
    mache = MachfromAreaRatio(AR2e, gamma, 1)
    PsePt2 = isenf.pressure_ratio(mache)
    PsePt1 = Ps1Pt1*Ps2Ps1/Ps2Pt2*PsePt2
    return PsePt1

def exit_shock(AR, gamma=1.4):
    """Calculates the pressure ratio Pe/Pt7 if there was a normal shock
    right at the exit of the nozzle. Returns the pressure ratio Pe/Pt7
    and the exit mach number"""
    mache = MachfromAreaRatio(AR, gamma)
    p2p1 = nsh.pressure_ratio(mache)
    p1pt1 = isenf.pressure_ratio(mache)
    PePt1 = p2p1*p1pt1
    return (PePt1, mache)

def exit_pressure_ratios(AR, gamma=1.4):
    """Calculates the nozzle exit pressure ratio for the two conditions
    if the throat is choked but remains subsonic to the exit, and if
    the flow does go supersonic with no shocks. Returns Ps/Ptsub, Ps/Ptsup"""
    machsub = MachfromAreaRatio(AR, gamma, 1)
    machsup = MachfromAreaRatio(AR, gamma)
    PsPtsub = isenf.pressure_ratio(machsub, gamma)
    PsPtsup = isenf.pressure_ratio(machsup, gamma)
    return (PsPtsub, PsPtsup)

def exit_pressure_delineation(AR, gamma=1.4):
    """Calculates the nozzle exit pressure ratio (PR) for the three
    conditions if the throat is choked but reamins subsonic to the exit,
    if the flow does go supersonic with no shocks, and if there is a
    shock in the nozzle. The flow is subsonic if the PR is equal to or
    below Ps/Ptsub, there is a shock in the nozzle if the PR is between
    Ps/Ptsub and PsPtshock, the flow is over expanded with oblique
    shocks outside the nozzle if the PR is between PsPtshock and
    Ps/Ptsup, the flow is perfectly expanded if the PR is equal to
    Ps/Ptsup and the flow is under expanded if the PR is below Ps/Ptsup.
    Returns Ps/Ptsub, PsPtshock, Ps/Ptsup"""
    
    PsPtshock, _ = exit_shock(AR, gamma)
    PsPtsub, PsPtsup = exit_pressure_ratios(AR, gamma)
    return (PsPtsub, PsPtshock, PsPtsup)