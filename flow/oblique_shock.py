
"""
Module containing oblique shock equations e.g. property ratios across
a shock, wave angle from Mach number and deflection angle
Beta is wave angle
Theta is deflection angle
All angles input as degrees and returned as degrees
"""

import numpy as np
import flow.isentropic_flow as isenflow

def normal_mach1(M, beta):
    """Returns the mach number BEFORE an oblique shock normal to the
    oblique shock"""
    Mn1 = M*np.sin(beta*np.pi/180)
    return Mn1

def normal_mach2(Mn1, gamma=1.4):
    """Returns the mach number AFTER an oblique shock normal to the
    oblique shock"""
    top = Mn1**2+(2/(gamma-1))
    bottom = 2*gamma/(gamma-1)*Mn1**2-1
    Mn2 = np.sqrt(top/bottom)
    return Mn2

def mach2(Mn2, beta, theta):
    """Returns Mach number after oblique shock from post-shock normal
    Mach number, and wave, and deflection angle. Give all angles in
    degrees"""
    M2 = Mn2/(np.sin((beta-theta)*np.pi/180))
    return M2

def static_pressure_ratio(Mn1, gamma=1.4):
    """Returns static pressure ratio Ps2/Ps1 from incident normal Mach
    number"""
    Ps2Ps1 = 1 + 2*gamma/(gamma+1)*(Mn1**2-1)
    return Ps2Ps1

def static_density_ratio(Mn1, gamma=1.4):
    """Returns static density ratio Ds2/Ds1 from incident normal Mach
    number"""
    Ds2Ds1 = (gamma+1)*Mn1**2/((gamma-1)*Mn1**2+2)
    return Ds2Ds1

def static_temperature_ratio(Mn1, gamma=1.4):
    """Returns static temperature ratio Ts2/Ts1 from incident normal Mach
    number"""
    Ts2Ts1 = (2*gamma*Mn1**2-(gamma-1))*((gamma-1)*Mn1**2+2)/((gamma+1)**2*Mn1**2)
    return Ts2Ts1

def total_pressure_ratio(Mn1, gamma=1.4):
    """Returns P02/P01, T2/T1, D2/D1, P2/P1"""
    Pt2Pt1 = (((gamma+1)*Mn1**2/((gamma-1)*Mn1**2+2))**(gamma/(gamma-1))*
    ((gamma+1)/(2*gamma*Mn1**2-gamma+1))**(1/(gamma-1)))
    return Pt2Pt1

def zero_OBM(theta, M1, gamma=1.4, n=0):
    """Returns the wave angle from the deflection angle and incident
    mach number.
    - Weak: n = 0, defualt because more common
    - Strong: n = 1
    """
    # derivation of equation found in link below
    # https://www.npworks.com/matlabcentral/fileexchange/32777-theta-beta-mach-analytic-relation
    theta = theta*np.pi/180
    mu = np.arcsin(1/M1)
    c = np.tan(mu)**2
    a = ((gamma-1)/2+(gamma+1)*c/2)*np.tan(theta)
    b = ((gamma+1)/2+(gamma+3)*c/2)*np.tan(theta)
    d = np.sqrt(4*(1-3*a*b)**3/((27*a**2*c+9*a*b-2)**2)-1)
    beta = np.arctan((b+9*a*c)/(2*(1-3*a*b))-(d*(27*a**2*c+9*a*b-2))/
                   (6*a*(1-3*a*b))*
                   np.tan(n*np.pi/3+1/3*np.arctan(1/d)))*180/np.pi
    return beta

def deflection_angle(beta, mach, gamma=1.4):
    """Returns the deflection angle"""
    beta *= np.pi/180
    top = 2*(mach**2*np.sin(beta)*np.sin(beta)-1)
    bottom = np.tan(beta)*(2+mach**2*(gamma+np.cos(2*beta)))
    theta = np.arctan(top/bottom)
    theta *= 180/np.pi
    return theta

def beta_max(mach, gamma=1.4):
    """Maximum wave angle for attached oblique shock wave. Plug in
    result into deflection_angle() for maximum
    deflection angle that will have an attached shock"""
    first = mach**2*(gamma-1)+4
    second = 16*(gamma+1)
    third = 8*mach**2*(gamma**2-1)
    fourth = mach**4*(gamma+1)**2
    fifth = 2*gamma*mach**2
    beta_max = 0.5*np.arccos((first-np.sqrt(second+third+fourth))/fifth)
    return beta_max*180/np.pi

def sonic_line(mach, gamma=1.4):
    """Returns the wave angle at which the flow behind a detached shock
    wave will be sonic. Plug this into deflection_angle() to obtain the
    angle on the surface at which the flow will be sonic"""
    first = (gamma-3)*mach**2 + (gamma+1)*mach**4
    second = 16*gamma*mach**4
    fourth = 4*gamma*mach**4
    beta_sonic = np.arcsin(np.sqrt((first+np.sqrt(second+first**2))/fourth))
    return beta_sonic*180/np.pi

def standoff_2d(mach):
    """Returns shock standoff distance for a cylinder in air with
    gamma = 1.4
    Source: Solomon, G. E., Shock Wave Detachment
    Distances for Plane and Axially Symmetric
    Flow, NACA Tech Note 3213 (1952)"""
    distance = 0.193*np.exp(4.67/mach**2)
    return distance

def curvature_2d(mach):
    """Returns normalized shock curvature (curvature/diameter of
    cylinder) for a cylinder in air with gamma = 1.4
    Source: Solomon, G. E., Shock Wave Detachment
    Distances for Plane and Axially Symmetric
    Flow, NACA Tech Note 3213 (1952)"""
    distance = 0.693*np.exp(1.8/(mach-1)**.75)
    return distance

def shape_2d(theta, mach, diameter, gamma=1.4):
    """Returns shock shape for a cylinder in air with gamma = 1.4
    Source: Solomon, G. E., Shock Wave Detachment
    Distances for Plane and Axially Symmetric
    Flow, NACA Tech Note 3213 (1952)"""
    distance = standoff_2d(mach)
    curv = curvature_2d(mach)
    beta = zero_OBM(theta, mach, gamma)
    y = np.linspace(0, 1.5*diameter)
    shape = 1/2*diameter+distance*diameter-curv*diameter*(np.sqrt(1+y**2*np.tan(beta)*np.tan(beta)/curv**2)-1)/(np.tan(beta)*np.tan(beta))

    return shape

def calc_mach2(mach1, theta, gamma=1.4):
    """Returns mach number after oblique shock wave and the wave angle"""
    beta = zero_OBM(theta, mach1, gamma)
    Mn1 = normal_mach1(mach1, beta)
    Mn2 = normal_mach2(Mn1, gamma)
    return mach2(Mn2, beta, theta)

def solve(M1, theta, gamma=1.4):
    """Returns post-shock Mach number and property ratios across the
    shock
    
    Rerturns
    --------
    0: post shock mach
    1 : wave angle
    2 : static pressure ratio Ps2/Ps1
    3 : static temperature ratio Ts2/Ts1
    4 : total pressure ratio Pt2/Pt1
    5 : static density ratio Ds2/Ds1"""
    beta = zero_OBM(theta, M1)
    Mn1 = normal_mach1(M1, beta)
    mach2 = calc_mach2(M1, theta, gamma)
    Ts2Ts1 = static_temperature_ratio(Mn1, gamma)
    Ds2Ds1 = static_density_ratio(Mn1, gamma)
    Ps2Ps1 = static_pressure_ratio(Mn1, gamma)
    Pt2Pt1 = total_pressure_ratio(Mn1, gamma)
    

    return (mach2, beta, Ps2Ps1, Ts2Ts1, Pt2Pt1, Ds2Ds1)