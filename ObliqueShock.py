
"""Module containing oblique shock equations e.g. property ratios across a shock, wave angle from Mach number and deflection angle"""

import numpy as np
import IsentropicFlow as IF

def normal_mach1(M, beta):
    """Returns the mach number BEFORE an oblique shock normal to the oblique shock"""
    Mn1 = M*np.sin(beta*np.pi/180)
    return Mn1

def normal_mach2(Mn1, gamma=1.4):
    """Returns the mach number AFTER an oblique shock normal to the oblique shock"""
    top = Mn1**2+(2/(gamma-1))
    bottom = 2*gamma/(gamma-1)*Mn1**2-1
    Mn2 = np.sqrt(top/bottom)
    return Mn2

def Mach2(Mn2, beta, theta):
    """Returns mach number after oblique shock. Give all angles in degrees"""
    M2 = Mn2/(np.sin((beta-theta)*np.pi/180))
    return M2

def press_ratio_Oshock(Mn1, gamma=1.4):
    """Returns P2/P1"""
    p2p1 = 1 + 2*gamma/(gamma+1)*(Mn1**2-1)
    return p2p1

def dens_ratio_Oshock(Mn1, gamma=1.4):
    """Returns rho2/rho1"""
    rho2rho1 = (gamma+1)*Mn1**2/((gamma-1)*Mn1**2+2)
    return rho2rho1

def temp_ratio_Oshock(Mn1, gamma=1.4):
    """Returns T2/T1, rho2/rho1, P2/P1"""
    p2p1 = press_ratio_Oshock(Mn1, gamma)
    rho2rho1 = dens_ratio_Oshock(Mn1, gamma)
    T2T1 = p2p1*rho2rho1**-1
    return T2T1, rho2rho1, p2p1

def stagpress_ratio_Oshock(Mn1, gamma=1.4):
    """Returns P02/P01, T2/T1, rho2/rho1, P2/P1"""
    T2T1, rho2rho1, p2p1 = temp_ratio_Oshock(Mn1, gamma)
    P02P01 = p2p1*(T2T1**-1)**(gamma/(gamma-1))
    return P02P01, T2T1, rho2rho1, p2p1

def zero_OBM(theta, M1, gamma=1.4, n = 0):
    """Returns the wave angle from the deflection angle and incident mach number. Weak: n = 0 Strong: n = 1"""
    # derivation of equation found in link below
    # https://www.npworks.com/matlabcentral/fileexchange/32777-theta-beta-mach-analytic-relation
    theta = theta*np.pi/180
    mu = np.asin(1/M1)
    c = np.tan(mu)**2
    a = ((gamma-1)/2+(gamma+1)*c/2)*np.tan(theta)
    b = ((gamma+1)/2+(gamma+3)*c/2)*np.tan(theta)
    d = np.sqrt(4*(1-3*a*b)**3/((27*a**2*c+9*a*b-2)**2)-1)
    beta = np.atan((b+9*a*c)/(2*(1-3*a*b))-(d*(27*a**2*c+9*a*b-2))/(6*a*(1-3*a*b))*np.tan(n*np.pi/3+1/3*np.atan(1/d)))*180/np.pi

    return beta

def solve_Oshock(M1, p1, T1, density1, theta, gamma=1.4):
    """Returns dictionary of all calculated values of the oblique shock from incident Mach number, pressure and temperature

    Dictionary Keys
    ---------------
    'beta'\n
    'p2'\n
    'T2'\n
    'density1'\n
    'density2'\n
    'Mn1'\n
    'Mn2'\n
    'M2'\n
    'p01'\n
    'p02'\n
    'T01'\n
    'T02'"""
    beta = zero_OBM(theta, M1, gamma)
    Mn1 = normal_mach1(M1, beta)
    p02p01, T2T1, rho2rho1, p2p1 = stagpress_ratio_Oshock(Mn1, gamma)
    p01 = IF.p0(p1, M1, gamma)
    p02 = p02p01*p01
    p2 = p2p1*p1
    T2 = T2T1*T1
    density2 = rho2rho1*density1
    Mn2 = normal_mach2(Mn1, gamma)
    M2 = Mach2(Mn2, beta, theta)
    T01 = IF.T0(T1, M1, gamma)
    T02 = IF.T0(T2, M2, gamma)
    
    result = {}
    result['beta'] = beta
    result['p2'] = p2
    result['T2'] = T2
    result['density1'] = density1
    result['density2'] = density2
    result['Mn1'] = Mn1
    result['Mn2'] = Mn2
    result['M2'] = M2
    result['p01'] = p01
    result['p02'] = p02
    result['T01'] = T01
    result['T02'] = T02
    
    return result

def print_solve_Oshock(result, decimals=4, state1='1', state2='2'):
    """Prints out ouput of solve_Oshock in a structured manner using SI units """

    print('Properties for states {} and {}'.format(state1, state2))
    print('Beta', round(result['beta'], decimals))
    print('Normal mach number before shock', round(result['Mn1'], decimals))
    print('Stagnation pressure before shock', round(result['p01'], decimals))
    print('Stagnation temperature before shock', round(result['T01'], decimals))

    print('\nMach number after shock', round(result['M2'], decimals))
    print('Normal mach number after shock', round(result['Mn2'], decimals))
    print('Pressure after shock', round(result['p2'], decimals))
    print('Temperature after shock', round(result['T2'], decimals))
    print('Density after shock', round(result['density2'], decimals))
    print('Stagnation pressure after shock', round(result['p02'], decimals))
    print('Stagnation temperature after shock', round(result['T02'], decimals))

    print('\nStagnation pressure ratio', round(result['p02']/result['p01'], decimals))
    print('Stagnation temperature ratio', round(result['T02']/result['T01'], decimals))