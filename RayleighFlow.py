
"""Set of equations for frictionless heat addition in a constant area
duct"""

def pressure_ratio(M1, gamma=1.4):
    """Calculates the inlet pressure to sonic pressure ratio"""
    Pt1Ptstar = (gamma+1)/(1+gamma*M1**2)
    return Pt1Ptstar

def temperature_ratio(M1, gamma=1.4):
    """Calculates the inlet temperature to sonic temperature ratio"""
    Ts1Tsstar = M1**2*((gamma+1)/(1+gamma*M1**2))**2
    return Ts1Tsstar

def density_ratio(M1, gamma=1.4):
    """Calculates the inlet desnity to sonic density ratio"""
    Ds1Dsstar = 1/M1**2*(1+gamma*M1**2)/(gamma+1)
    return Ds1Dsstar

def total_temperature_ratio(M1, gamma=1.4):
    """Calculates the inlet total temperature to sonic total temperature
    ratio"""
    Tt1Ttstar = M1**2*((gamma+1)/(1+gamma*M1**2))**2*(
        1+(gamma-1)/2*M1**2)/((gamma+1)/2)
    return Tt1Ttstar

def total_pressure_ratio(M1, gamma=1.4):
    """Calculates the inlet total pressure to sonic total pressure ratio"""
    Pt1Ptstar = M1**2*((gamma+1)/(1+gamma*M1**2))**2*(
        (1+(gamma-1)/2*M1**2)/((gamma+1)/2))**(gamma/(gamma+1))
    return Pt1Ptstar

def heat_addition(M1, Tt1, cp=1000, gamma=1.4):
    """Calculates the heat addition to drive the flow to sonic
    conditions. In a real flow if the heat addition is greater than the
    flow at the inlet will spill outside of the duct. Heat is returned
    as per unit mass"""
    qstar = (total_temperature_ratio(M1, gamma)-1)*cp*Tt1
    return qstar
