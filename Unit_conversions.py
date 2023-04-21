
import numpy as np

def convert_temperature(temps, to:str='K'):
    """Convert temperature between Kelvin and Rankine
    to = 'K' for converting to Kelvin
    to = 'R' for converting to Rankine"""
    factor = 1.8
    if to == 'R':
        try:
            temps *= factor
            return temps
        except:
            return np.array(temps)*factor
    elif to == 'K':
        try:
            temps = temps/factor
            return temps
        except:
            return np.array(temps)/factor
    else:
        print('Did not convert temperature')

def convert_mass(mass, to:str='kg'):
    """Convert mass between lbm and kg
    to = 'lbm' for converting to lbm
    to = 'kg' for converting to kg"""
    factor = 2.20462
    if to == 'lbm':
        try:
            mass *= factor
            return mass
        except:
            return np.array(mass)*factor
    elif to == 'kg':
        try:
            mass = mass/factor
            return mass
        except:
            return np.array(mass)/factor
    else:
        print('Did not convert mass')

def convert_pressure(pressures, to:str='Pa'):
    """Convert pressure between Pascals and PSI
    to = 'Pa' for converting to Pascals
    to = 'psi' for converting to PSI"""
    factor = 6895.0
    if to == 'psi':
        try:
            pressures = pressures/factor
            return pressures
        except:
            return np.array(pressures)/factor
    elif to == 'Pa':
        try:
            pressures = pressures*factor
            return pressures
        except:
            return np.array(pressures)*factor
    else:
        print('Did not convert pressure')

def convert_energy(energy, to:str='J'):
    """Convert mass specific energy/work between Btu/lbm and J/kg
    to = 'J' for converting to J/kg
    to = 'BTU' for converting to BTU/lbm"""
    factor = 1055 * 2.205
    if to == 'BTU':
        try:
            energy = energy/factor
            return energy
        except:
            return np.array(energy)/factor
    elif to == 'J':
        try:
            energy = energy*factor
            return energy
        except:
            return np.array(energy)*factor
    else:
        print('Did not convert energy')

def convert_force(force, to:str='N'):
    """Convert force between Newtons and lbf
    to = 'N' for converting to Newtons
    to = 'lbf' for converting to lbf"""
    factor = 4.44822
    if to == 'lbf':
        try:
            force = force/factor
            return force
        except:
            return np.array(force)/factor
    elif to == 'N':
        try:
            force = force*factor
            return force
        except:
            return np.array(force)*factor
    else:
        print('Did not convert force')

def convert_length(length, to:str='m', from:str='in'):
    """Convert length from meters to inches or meters to feet NOT between feet and inches
    to = 'm' for converting to meters
    to = 'in' for converting to inches
    to = 'ft' for converting to feet"""
    factor1 = .0254
    factor2 = .0254*12
    if to == 'm':
        try:
            length = length/factor1
            return length
        except:
            return np.array(length)/factor1
    elif to == 'in':
        try:
            length = length*factor1
            return length
        except:
            return np.array(length)*factor1
    elif to == 'ft':
        try:
            length = length*factor2
            return length
        except:
            return np.array(length)*factor2
    else:
        print('Did not convert length')

def convert_area(area, to:str='SI'):
    """Convert length from sq meters to sq inches or sq meters to sq feet NOT between sq feet and sq inches
    to = 'm' for converting to sq meters
    to = 'in' for converting to sq inches
    to = 'ft' for converting to sq feet"""
    factor1 = .0254**2
    factor2 = (.0254*12)**2
    if to == 'm':
        try:
            area = area/factor1
            return area
        except:
            return np.array(area)/factor1
    elif to == 'in':
        try:
            area = area*factor1
            return area
        except:
            return np.array(area)*factor1
    elif to == 'ft':
        try:
            area = area*factor2
            return area
        except:
            return np.array(area)*factor2
    else:
        print('Did not convert area')

def convert_speed(speed, to:str='SI'):
    """Convert length from meters/sec to inches/sec or meters/sec to feet/sec NOT between feet and inches
    to = 'm' for converting to meters/sec
    to = 'in' for converting to inches/sec
    to = 'ft' for converting to feet/sec"""
    factor1 = .0254**2
    factor2 = (.0254*12)**2
    if to == 'm':
        try:
            speed = speed/factor1
            return speed
        except:
            return np.array(speed)/factor1
    elif to == 'in':
        try:
            speed = speed*factor1
            return speed
        except:
            return np.array(speed)*factor1
    elif to == 'ft':
        try:
            speed = speed*factor2
            return speed
        except:
            return np.array(speed)*factor2
    else:
        print('Did not convert speed')


