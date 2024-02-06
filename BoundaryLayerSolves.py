import numpy as np

def Y1(y2):
    """Add in definitions"""
    return y2

def Y2(y3):
    """Add in definitions"""
    return y3

def Y3(y1,y3,y4,y5,Sm,Tinf):
    """Add in definitions"""
    return -y3*((y5/(2*(y4)))-(y5/(y4+Sm/Tinf)))-y1*y3*((y4+Sm/Tinf)/(np.sqrt(y4)*(1+Sm/Tinf)))

def Y4(y5):
    """Add in definitions"""
    return y5

def Y5(y1,y3,y4,y5,Sm,Tinf,Minf,Pr,gamma):
    """Add in definitions"""
    return -y5**2*((0.5/y4)-(1/(y4+Sm/Tinf)))-Pr*y1*y5/np.sqrt(y4)*(y4+Sm/Tinf)/(1+Sm/Tinf)-(gamma-1)*Pr*Minf**2*y3**2

def rungekutta(N,h,y1,y2,y3,y4,y5,Sm,Tinf,Pr,gamma,Minf):
    """Fourth order Runge-Kutta method specifically for compressible
    boundary layer flow functions"""
    for i in range(N):
        k11 = Y1(y2[i])
        k21 = Y2(y3[i])
        k31 = Y3(y1[i], y3[i], y4[i], y5[i],Sm,Tinf)
        k41 = Y4(y5[i])
        k51 = Y5(y1[i], y3[i], y4[i], y5[i],Sm,Tinf,Minf,Pr,gamma)
        
        
        k12 = Y1(y2[i]+0.5*h*k21)
        k22 = Y2(y3[i]+0.5*h*k31)
        k32 = Y3(y1[i]+0.5*h*k11, y3[i]+0.5*h*k31, y4[i]+0.5*h*k41, y5[i]+0.5*h*k51,Sm,Tinf)
        k42 = Y4(y5[i]+0.5*h*k51)
        k52 = Y5(y1[i]+0.5*h*k11, y3[i]+0.5*h*k31, y4[i]+0.5*h*k41, y5[i]+0.5*h*k51,Sm,Tinf,Minf,Pr,gamma)
        
        k13 = Y1(y2[i]+0.5*h*k22)
        k23 = Y2(y3[i]+0.5*h*k32)
        k33 = Y3(y1[i]+0.5*h*k12, y3[i]+0.5*h*k32, y4[i]+0.5*h*k42, y5[i]+0.5*h*k52,Sm,Tinf)
        k43 = Y4(y5[i]+0.5*h*k52)
        k53 = Y5(y1[i]+0.5*h*k12, y3[i]+0.5*h*k32, y4[i]+0.5*h*k42, y5[i]+0.5*h*k52,Sm,Tinf,Minf,Pr,gamma)
        
        k14 = Y1(y2[i]+h*k23)
        k24 = Y2(y3[i]+h*k33)
        k34 = Y3(y1[i]+h*k13, y3[i]+h*k33, y4[i]+h*k43, y5[i]+h*k53,Sm,Tinf)
        k44 = Y4(y5[i]+h*k53)
        k54 = Y5(y1[i]+h*k13, y3[i]+h*k33, y4[i]+h*k43, y5[i]+h*k53,Sm,Tinf,Minf,Pr,gamma)
        
        y5[i+1] = y5[i] + (1/6)*(k51 + 2*k52 + 2*k53 + k54)*h
        y4[i+1] = y4[i] + (1/6)*(k41 + 2*k42 + 2*k43 + k44)*h
        y3[i+1] = y3[i] + (1/6)*(k31 + 2*k32 + 2*k33 + k34)*h 
        y2[i+1] = y2[i] + (1/6)*(k21 + 2*k22 + 2*k23 + k24)*h
        y1[i+1] = y1[i] + (1/6)*(k11 + 2*k12 + 2*k13 + k14)*h

    return y1,y2,y3,y4,y5

def compressible_solve(Minf=1.0, Tinf=300, TwTe=2.0, etamax=10, adi=1, Pr=0.72, gamma=1.4, N=100, precision=1e-9):
    """
    Source: https://github.com/frkanz/A-CFD-Tutorial-in-Julia-Compressible-Blasius

    Parameters
    ----------
    Minf : mach number at the edge of the boundary layer
    Tinf : temperature at the edge of the boundary layer
    TwTe : ratio of wall temperature to boundary layer edge temperature,
    only necessary if wall is not adiabatic
    adi : specify 1 for an adiabatic wall, 0 if wall is not adiabatic
    gamma : ratio of specific heats
    PR : Prandtl's number
    etamax : dimensionless distance away from the wall
    N : number of vertical divisions at which to split the boundary layer
    precision : for Newton iteration convergence

    Returns
    -------
    eta : vertical similarity parameter defined below
    y : actual vertical distance
    U : u/ue dimensionless velocity
    T : T/Te dimensionless temperature
    N : number of vertical divisions

    Definitions
    -----------
    eta : ue*De/sqrt(2s)*int([0 y](Te/T)dy)
    s : Ve*De*ue*x
    Ve : viscosity at BL edge
    De : density at BL edge
    ue : velocity at BL edge

    """
    deltaeta = etamax/N
    Sm = 110.4 # Sutherland law coefficient for [Kelvin]
    delta = precision/10 # Small number for shooting method

    # Initializing the solution vectors
    y1 = np.zeros(N+1) # f
    y2 = np.zeros(N+1) # f'
    y3 = np.zeros(N+1) # f''
    y4 = np.zeros(N+1) # ρ(eta)
    y5 = np.zeros(N+1) # ρ(eta)'
    eta  = [i*deltaeta for i in range(N+1)]

    if adi == 1:
        y1[0] = 0
        y2[0] = 0
        y5[0] = 0
        alpha0 = 0.1 # Initial Guess
        beta0 = 2.3118 # Initial Guess

    elif adi == 0:
        y1[0] = 0
        y2[0] = 0
        y4[0] = TwTe
        alpha0 = 0.1 # Initial Guess
        beta0 = 3.0 # Initial Guess

    maxiternum = 500
    y2n = 1
    y4n = 1
    for i in range(maxiternum):
        # Solve one
        if adi == 1:
            y1[0] = 0
            y2[0] = 0
            y5[0] = 0
            y3[0] = alpha0 # Initial Guess
            y4[0] = beta0 # Initial Guess
        elif adi == 0:
            y1[0] = 0
            y2[0] = 0
            y4[0] = TwTe

            y3[0] = alpha0 # Initial Guess
            y5[0] = beta0 # Initial Guess
    
		# First solution for Newton's iteration
        y1,y2,y3,y4,y5 = rungekutta(N,deltaeta,y1,y2,y3,y4,y5,Sm,Tinf,Pr,gamma,Minf)

        # Storing the freestream values for Newton's iteration method
        y20 = y2[N]
        y40 = y4[N]
        
        # Small number addition for Newton's iteration method
        if adi == 1:
            y1[0] = 0
            y2[0] = 0
            y5[0] = 0

            y3[0] = alpha0+delta # Initial Guess + Small number
            y4[0] = beta0 # Initial Guess
        elif adi == 0:
            y1[0] = 0
            y2[0] = 0
            y4[0] = TwTe

            y3[0] = alpha0+delta # Initial Guess + Small number
            y5[0] = beta0 # Initial Guess
    
		# Second solution for Newton's iteration
        y1,y2,y3,y4,y5 = rungekutta(N,deltaeta,y1,y2,y3,y4,y5,Sm,Tinf,Pr,gamma,Minf)

        # Storing the freestream values for Newton's iteration method
        y201 = y2[N]
        y401 = y4[N]

        # Small number addition for Newton's iteration method
        # Small number addition for Newton's iteration method
        if adi == 1:
            y1[0] = 0
            y2[0] = 0
            y5[0] = 0
            y3[0] = alpha0 # Initial Guess
            y4[0] = beta0+delta # Initial Guess + Small number
        elif adi == 0:
            y1[0] = 0
            y2[0] = 0
            y4[0] = TwTe
            y3[0] = alpha0 # Initial Guess
            y5[0] = beta0+delta # Initial Guess + Small number
    
		# Third solution for Newton's iteration
        y1,y2,y3,y4,y5 = rungekutta(N,deltaeta,y1,y2,y3,y4,y5,Sm,Tinf,Pr,gamma,Minf)

        # Storing the freestream values for Newton's iteration method
        y202 = y2[N]
        y402 = y4[N]

        # Calculation of the next initial guess with Newton's iteration method
        p11 = (y201-y20)/delta
        p21 = (y401-y40)/delta
        p12 = (y202-y20)/delta
        p22 = (y402-y40)/delta
        r1 = 1-y20
        r2 = 1-y40
        deltaalpha = (p22*r1-p12*r2)/(p11*p22-p12*p21)
        deltabeta = (p11*r2-p21*r1)/(p11*p22-p12*p21)
        alpha0 = alpha0 + deltaalpha
        beta0 = beta0 + deltabeta

        if ((np.abs(y2[N]-1.0)<precision) and 
            (np.abs(y4[N]-1.0)<precision) and 
            (np.abs(y2n-np.linalg.norm(y2))<precision) and 
            (np.abs(y4n-np.linalg.norm(y4))<precision)):
            break
    
        y2n = np.linalg.norm(y2)
        y4n = np.linalg.norm(y4)

	# Copying values for logical names
    U = y2
    T = y4

    # Integration for eta --> y transformation
    y = np.zeros(N+1)
    for i in range(N):
        y[i] = y[i-1] + y4[i]*(eta[i]-eta[i-1])

    y = y*np.sqrt(2)

    return eta,y,U,T,N

def dfpn(n, f, w1, w2):
    """Add in definitions"""
    return w1
def dw1n(n, f, w1, w2):
    """Add in definitions"""
    return w2
def dw2n(n, f, w1, w2):
    """Add in definitions"""
    return -f*w2/2

def rungekutta(N, n, f, w1, w2, step):
    """Fourth order Runge-Kutta method specifically for incompressible
    boundary layer flow functions"""
    for i in range(N):
        K1f  = dfpn(n[i], f[i], w1[i], w2[i])
        K1w1 = dw1n(n[i], f[i], w1[i], w2[i])
        K1w2 = dw2n(n[i], f[i], w1[i], w2[i])
    
        K2f  = dfpn(n[i]+1/2, f[i]+K1f*step/2, w1[i]+K1w1*step/2, w2[i]+K1w2*step/2)
        K2w1 = dw1n(n[i]+1/2, f[i]+K1f*step/2, w1[i]+K1w1*step/2, w2[i]+K1w2*step/2)
        K2w2 = dw2n(n[i]+1/2, f[i]+K1f*step/2, w1[i]+K1w1*step/2, w2[i]+K1w2*step/2)
    
        K3f  = dfpn(n[i]+1/2, f[i]+K2f*step/2, w1[i]+K2w1*step/2, w2[i]+K2w2*step/2)
        K3w1 = dw1n(n[i]+1/2, f[i]+K2f*step/2, w1[i]+K2w1*step/2, w2[i]+K2w2*step/2)
        K3w2 = dw2n(n[i]+1/2, f[i]+K2f*step/2, w1[i]+K2w1*step/2, w2[i]+K2w2*step/2)
    
        K4f  = dfpn(n[i], f[i]+K3f*step, w1[i]+K3w1*step, w2[i]+K3w2*step)
        K4w1 = dw1n(n[i], f[i]+K3f*step, w1[i]+K3w1*step, w2[i]+K3w2*step)
        K4w2 = dw2n(n[i], f[i]+K3f*step, w1[i]+K3w1*step, w2[i]+K3w2*step)
        
        n[i+1]  = n[i] + step
        f[i+1]  = f[i]  + (step/6)*(K1f  + 2*K2f  + 2*K3f  + K4f)
        w1[i+1] = w1[i] + (step/6)*(K1w1 + 2*K2w1 + 2*K3w1 + K4w1)
        w2[i+1] = w2[i] + (step/6)*(K1w2 + 2*K2w2 + 2*K3w2 + K4w2)
    
    return n, f, w1, w2


def incompressible_solve(precision=1e-6):
    step = .4
    stop = 10
    N = int(stop/step)

    n = np.zeros(N+1)
    f = np.zeros(N+1)
    w1 = np.zeros(N+1)
    w2 = np.zeros(N+1)

    f[0] = 0
    w1[0] = 0

    w21 = .2 # Just some random initial guess
    precision = 1e-6
    delta = precision/10

    maxiternum = 500
    for i in range(maxiternum):
        w2[0] = w21

        v = rungekutta(N, n, f, w1, w2, step)[2][-1]
        w2[0] += delta
        t = rungekutta(N, n, f, w1, w2, step)[2][-1]

        y1 = (t-v)/delta
        dx = (1-v)/y1
        w21 += dx
        
        if (np.abs(w1[-1]-1.0)<precision):
            break

    return