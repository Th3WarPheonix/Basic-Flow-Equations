import numpy as np
import matplotlib.pyplot as plt

def four_digit(x:np.ndarray, thickness:float, maxcamber:float, camber_loc:float, opt='base'):
    """
    NACA 4-digit series airfoil. Thickness is always maximum at 30% of
    chord. Use option 'cfd' for a sharp trailing edge with minimal
    deviation from base airfoil.

    Parameters
    ----------
    thickness : maximum thickness as percentage of chord
    maxcamber : maximum camber as percentage of chord
    camber_loc : location of maximum camber as percentage of chord
    
    Returns
    -------
    upperx : x points for upper surface
    uppery : y points for upper surface
    lowerx : x points for lower surface
    lowery : y points for lower surface
    """

    firsthalf = (x<=camber_loc).astype(int)
    sechalf = (x>=camber_loc).astype(int)

    match opt:
        case 'base':
            a0 = 0.2969
            a1 = -0.1260
            a2 = -0.3516
            a3 = 0.2843
            a4 = -0.1015

        case 'cfd':
            a0 = 0.2969
            a1 = -0.1260
            a2 = -0.3516
            a3 = 0.2843
            a4 = -0.1036
    
    # half thickness at a given x
    yt = thickness/0.2*(a0*np.sqrt(x) + a1*x + a2*x**2 + a3*x**3 + a4*x**4)

    # mean camber line
    yc = (maxcamber/camber_loc**2*(2*camber_loc*x-x**2))*firsthalf + maxcamber/(1-camber_loc)**2*((1-2*camber_loc)+2*camber_loc*x-x**2)*sechalf
    # derivative of the mean camber line
    dycdx = (2*maxcamber/camber_loc**2*(camber_loc-x))*firsthalf + 2*maxcamber/(1-camber_loc)**2*(camber_loc-x)*sechalf
    theta = np.arctan(dycdx)

    upperx = x - yt*np.sin(theta)
    lowerx = x + yt*np.sin(theta)

    uppery = yc + yt*np.cos(theta)
    lowery = yc - yt*np.cos(theta)

    return (upperx, uppery, lowerx, lowery)

def four_digit_modified(x:np.ndarray, thickness:float, maxcamber:float, camber_loc:float):
    """
    NACA 4-digit modified series airfoil. 

    Parameters
    ----------
    thickness : maximum thickness as percentage of chord
    maxcamber : maximum camber as percentage of chord
    camber_loc : location of maximum camber as percentage of chord
    
    Returns
    -------
    upperx : x points for upper surface
    uppery : y points for upper surface
    lowerx : x points for lower surface
    lowery : y points for lower surface
    """

    firsthalf = (x<=camber_loc).astype(int)
    sechalf = (x>=camber_loc).astype(int)


    a0 = 0.2969
    a1 = -0.1260
    a2 = -0.3516
    a3 = 0.2843
    a4 = -0.1015

    d0 = 0.2
    d1 = 0.234
    d2 = 0.315
    d3 = 0.465
    d4 = 0.7

    
    # half thickness at a given x
    yt = thickness/0.2*(a0*np.sqrt(x) + a1*x + a2*x**2 + a3*x**3 + a4*x**4)*firsthalf
    yt += thickness/0.20*(d0 + d1*(1-x) + d2*(1-x)**2 + d3*(1-x)**3)*sechalf

    # the camber line equation can be subsituted for another function
    # but it will no longer be a four digit airfoil
    # mean camber line
    yc = (maxcamber/camber_loc**2*(2*camber_loc*x-x**2))*firsthalf 
    yc += maxcamber/(1-camber_loc)**2*((1-2*camber_loc)+2*camber_loc*x-x**2)*sechalf
    # derivative of the mean camber line
    dycdx = (2*maxcamber/camber_loc**2*(camber_loc-x))*firsthalf + 2*maxcamber/(1-camber_loc)**2*(camber_loc-x)*sechalf
    theta = np.arctan(dycdx)

    upperx = x - yt*np.sin(theta)
    lowerx = x + yt*np.sin(theta)

    uppery = yc + yt*np.cos(theta)
    lowery = yc - yt*np.cos(theta)

    return (upperx, uppery, lowerx, lowery)

def five_digit(x:np.ndarray, thickness:float, maxcamber:float, camber_loc:float,  opt='base'):
    """
    NACA 5-digit series airfoil. Thickness is always maximum at 30% of
    chord. Use option 'cfd' for a sharp trailing edge with minimal
    deviation from base airfoil.

    Parameters
    ----------
    thickness : maximum thickness as percentage of chord
    maxcamber : maximum camber as percentage of chord
    camber_loc : location of maximum camber as percentage of chord
    
    Returns
    -------
    upperx : x points for upper surface
    uppery : y points for upper surface
    lowerx : x points for lower surface
    lowery : y points for lower surface
    """

    firsthalf = (x<=camber_loc).astype(int)
    sechalf = (x>=camber_loc).astype(int)

    match opt:
        case 'base':
            a = 0.2969
            b = -0.1260
            c = -0.3516
            d = 0.2843
            e = -0.1015

        case 'cfd':
            a = 0.2969
            b = -0.1260
            c = -0.3516
            d = 0.2843
            e = -0.1036
    
    # half thickness at a given x
    yt = thickness*(a*np.sqrt(x) + b*x + c*x**2 + d*x**3 + e*x**4)/0.2

    # mean camber line
    yc = (maxcamber/camber_loc**2*(2*camber_loc*x-x**2))*firsthalf + maxcamber/(1-camber_loc)**2*((1-2*camber_loc)+2*camber_loc*x-x**2)*sechalf
    # derivative of the mean camber line
    dycdx = (2*maxcamber/camber_loc**2*(camber_loc-x))*firsthalf + 2*maxcamber/(1-camber_loc)**2*(camber_loc-x)*sechalf
    theta = np.arctan(dycdx)

    upperx = x - yt*np.sin(theta)
    lowerx = x + yt*np.sin(theta)

    uppery = yc + yt*np.cos(theta)
    lowery = yc - yt*np.cos(theta)

    return (upperx, uppery, lowerx, lowery)

def generate_coefficients():
    """This function generates the coefficients for the equations for
    the NACA airfoils as described in "Computer Program To Obtain
    Ordinates for NACA Airfoils", and hopefully gives enough
    information to generate the coefficients for other airfoils.
    (https://ntrs.nasa.gov/api/citations/19970008124/downloads/19970008124.pdf)
    """

    # NACA 4-digit airfoil constraints
    # maximum ordinate
    # x/c = 0.30 y/c = 0.10
    # x/c = 0.30 dy/dx = 0
    # ordinate at trailing edge
    # x/c = 1.00 y/c = 0.002
    # Magnitude of trailing edge angle
    # x/c = 1 |dy/dx| = 0.234
    # Nose shape
    # x/c = 0.10 y/c = 0.078
    x1 = 0.3
    y1 = 0.1
    dydx1 = 0

    x2 = 1
    y2 = 0.002
    dydx2 = -0.234

    x3 = 0.1
    y3 = 0.078
    A = np.array(
        [[x1**(1/2), x1, x1**2, x1**3, x1**4],
        [0.5*x1**(-1/2), 1, 2*x1, 3*x1**2, 4*x1**3],
        [x2**(1/2), x2, x2**2, x2**3, x2**4],
        [0.5*x2**(-1/2), 1, 2*x2, 3*x2**2, 4*x2**3],
        [x3**(1/2), x3, x3**2, x3**3, x3**4]])
    b = np.array(
        [y1, dydx1, y2, dydx2, y3])
    
    # NACA 4-digit modified airfoil
    radius_of_curvature = (1+(0.5*x1**(-1/2)+1+2*x1+3*x1**2+4*x1**3)**2)**(3/2)/((-0.25*x1**(-3/2)+2+6*x1+12*x1))

    m = 0.2 # the point at which the two portions of the airfoil meet
    d1 = np.interp(m, [0.2, 0.3, 0.4, 0.5, 0.6], [0.2, 0.234, 0.315, 0.465, .070])
    
    x1 = m
    y1 = 0.1
    dydx1 = 0
    R = (1-m)**2/(2*d1*(1-m)-0.588)

    x2 = 1
    y2 = 0.002
    dydx2 = d1

    x3 = 0.1
    y3 = 0.078
    A = np.array(
        [[x1**(1/2), x1, x1**2, x1**3, x1**4],
        [0.5*x1**(-1/2), 1, 2*x1, 3*x1**2, 4*x1**3],
        [(1+(0.5*x1**(-1/2)+1+2*x1+3*x1**2+4*x1**3)**2)**(3/2)/((-0.25*x1**(-3/2)+2+6*x1+12*x1))]
        [x2**(1/2), x2, x2**2, x2**3, x2**4],
        [0.5*x2**(-1/2), 1, 2*x2, 3*x2**2, 4*x2**3],
        [x3**(1/2), x3, x3**2, x3**3, x3**4]])
    b = np.array(
        [y1, dydx1, y2, dydx2, y3])
    
    coefficients = np.linalg.solve(A, b)
    return coefficients

def main():
    # thickness = .12 # decimal
    # camber = .09 # decimal
    # camberloc = .4 # decimal
    # x = np.linspace(0, 1, 1000)
    # upperx, uppery, lowerx, lowery = four_digit(x, thickness, camber, camberloc)
    # plt.plot(upperx, uppery)
    # plt.plot(lowerx, lowery)
    # plt.axis('equal')
    # plt.title(f'NACA {int(camber*100):01d}{int(camberloc*10):01d}{int(100*thickness):02d}')
    # plt.show()

    print(generate_coefficients())


if __name__ == '__main__':
    main()