import numpy as np
import matplotlib.pyplot as plt

def fdJ(X, syseqs, xidx=None, fxidx=None, h:float=1e-10, **inputs):
    """
    Calculates the Jacobian for a system of equations via a forward
    finite difference. Allows separate specification of which inputs are
    to be perturbed and which outputs are to be considered. This
    function can be used with a shooting method and if an n x m input
    array X is given then X[xidx, 0] will be perturbed by h and X[fxidx,
    -1] will be considered as the ouput of interest
    
    Parameters
    ----------
    X : a column vector of that will based into the system of equations
    syseqs : a used-defined function that takes in a vector the same 
        size as X and outputs a vector the same size as X
    xidx : indices of X that are to be perturbed, 
        not a boolean matrix use actual indices
    fxidx : indices of the output of the system of equations that are to
        be considered when calculating the result of the perturbation,
        not a boolean matrix use actual indices
    h : the amount of perturbation inputs : (optional) any inputs that
    must be passed into the 
        system of equations

    Example
    -------
    fdJ(X, FX = F(X), h = 10**-4, n = size(X,1))

    (Forward) Finite Difference approximation of the Jacobian. Requires
    starting vector X. Optional h and n. FX is optional in case F(X) was
    already computed for efficiency. # Example ```julia-repl julia>
    fdJ([1.5;2.0]) 2Xmid Array{Float64,2}: 3.0001  -1.0 1.0     -4.0001
    ``` Source:
    https://github.com/osveliz/numerical-veliz/blob/master/src/systems/QuasiNewton.jl
    
    """

    xidx = np.array(xidx)
    fxidx = np.array(fxidx)
    if xidx is not None:
        if xidx.ndim == 1:
            xidx = [[xidxi] for i, xidxi in enumerate(xidx)]
            fxidx = [[fxidxi] for i, fxidxi in enumerate(fxidx)]
        nx = len(xidx)
        nf = len(fxidx)
    else:
        print('Specify input and output indices')
        return
        
    # Establish arrays in memory
    H = np.zeros_like(X[:, [0]])
    J = np.zeros((nf, nx))
    # Get initial function value
    FX = syseqs(X, **inputs)
    # Loop through the indices that are to be perturbed and 
    # calculate the difference
    for i, xidxi in enumerate(xidx):
        H[[xidxi], 0] = h
        J[:, [i]] = syseqs(X + H, **inputs)[fxidx, [-1]] - FX[fxidx, [-1]]
        H = np.zeros_like(X[:, [0]])

    J = J/h # removes placeholder and divide by h
    return J

def broyden2(X, syseqs, goal=0, xidx=None, fxidx=None, precision:float=1e-10, 
             maxiters:int=500, **inputs):
    """
    The function solves a system of nonlinear equations using Broyden's
    method. This can be used to solve shooting method problems if a
    Runge-Kutta-like system of equations is passed in

    Parameters
    ----------
    X : starting guess for iterative procedure, must be a column vector
    syseqs : a function that contains the system of equations for
        solving, output must be a column vector
    goal : a vector of equal size to xidx of the values the system of
        equations should converge to
    xidx : indices of X that are to be perturbed, 
        not a boolean matrix use actual indices
    midx : indices of X that need to be matched in the middle shoots, 
        should be the equations that are one below the highest derivative
    fxidx : indices of the output of the system of equations that are to
        be considered when calculating the result of the perturbation, 
        not a boolean matrix use actual indices
    **inputs : other inputs that will be passed into the system of equations

    Notes
    -----
    This function implements the "Good" Broyden method for a system of
    nonlinear equations. The initial J approximated by differencing, the
    updated ones are determined by the Sherman-Morrison formula

    Example
    -------
    X = [[0], [0], [1.5], [.75], [.5]]
    syseqs = [[X[1]], [X[1]**2+X[3]**2], [X[2]+5-X[4], X[0]+X[2], X[3]]]
    goal = [[1], [2]]
    xidx = [[2], [4]]
    fxidx = [[3], [4]]
    The exmaple inputs will be interpreted by this function as follows:
    1. X will be input into syseqs for the base value
    2. X[2] and X[4] will be perturbed by a really small an donput into syseqs
    3. syseqs[3] and syseqs[4] will be used to calculate the Jacobian
    4. X[2] and X[4] will be updated via Broyden's method
    5. Repeat steps uintil convergence


    # Example
    ```julia-repl
    julia> Broyden([1.5;2.0])
    2-element Array{Float64,1}:
    1.6180340721138766
    1.6180340507162783
    ```
    Source: 
    https://github.com/osveliz/numerical-veliz/blob/master/src/systems/QuasiNewton.jl
    
    """
    if xidx:
        if len(xidx) != len(fxidx):
            raise('List of indices of input and output must be same length')
        xidx = [[xidxi] for i, xidxi in enumerate(xidx)]
        fxidx = [[fxidxi] for i, fxidxi in enumerate(fxidx)]
    else:
        print('Specify input and output indices')
        return
    
    # y2n = 1 # Can be used for judging convergence
    # y4n = 1 # Can be used for judging convergence
    # get the base value and numerical get the jacobian
    FX = syseqs(X, **inputs)
    invJ = np.linalg.inv(fdJ(X, syseqs, xidx, fxidx,**inputs))
    for iter in range(maxiters):
        Xold_jac = X[xidx, [0]]
        X_jac = Xold_jac - np.matmul(invJ, FX[fxidx, [-1]]-goal)
        deltaX_jac = X_jac - Xold_jac
        X[xidx, [0]] = X_jac

        FXold_jac = FX[fxidx, [-1]]
        FX = syseqs(X, **inputs)
        deltaF_jac = FX[fxidx, [-1]] - FXold_jac
        trans = deltaX_jac[:,0] # tranpose
        # Next line is Broydens method, it is just a way to update the
        # jacobian without doing all of the perturbations
        invJ = invJ + (deltaX_jac - np.matmul(invJ,deltaF_jac))/np.matmul(np.matmul(trans,invJ), deltaF_jac)*np.matmul(trans,invJ)
        
        if np.all(np.abs(FX[fxidx, -1]-goal)<precision): #and
            # (np.abs(y2n-np.linalg.norm(FX[1, -1]))<precision) and 
            # (np.abs(y4n-np.linalg.norm(FX[3, -1]))<precision)):
            break
        # y2n = np.linalg.norm(FX[1, -1])
        # y4n = np.linalg.norm(FX[3, -1])
        
    return X, iter

def S1(fp):
    """Add in definitions"""
    return fp

def S2(fpp):
    """Add in definitions"""
    return fpp

def S3(f,fpp,g,gp,Sm,Tinf):
    """Add in definitions"""
    return -fpp*(gp/(2*g)-gp/(g+Sm/Tinf))-f*fpp*((g+Sm/Tinf)/(np.sqrt(g)*(1+Sm/Tinf)))

def S4(gp):
    """Add in definitions"""
    return gp

def S5(f,fpp,g,gp,Sm,Tinf,Minf,Pr,gamma):
    """Add in definitions"""
    return -gp**2*(0.5/g-1/(g+Sm/Tinf))-Pr*f*gp/np.sqrt(g)*(g+Sm/Tinf)/(1+Sm/Tinf)-(gamma-1)*Pr*Minf**2*fpp**2

def rungekutta_sphere(N,h,f,fp,fpp,g,gp,Sm,Tinf,Pr,gamma,Minf):
    """Fourth order Runge-Kutta method specifically for compressible
    boundary layer flow functions"""
    for i in range(N):
        k11 = S1(fp[i])
        k21 = S2(fpp[i])
        k31 = S3(f[i], fpp[i], g[i], gp[i],Sm,Tinf)
        k41 = S4(gp[i])
        k51 = S5(f[i], fpp[i], g[i], gp[i],Sm,Tinf,Minf,Pr,gamma)
        
        
        k12 = S1(fp[i]+0.5*h*k21)
        k22 = S2(fpp[i]+0.5*h*k31)
        k32 = S3(f[i]+0.5*h*k11, fpp[i]+0.5*h*k31, g[i]+0.5*h*k41, gp[i]+0.5*h*k51,Sm,Tinf)
        k42 = S4(gp[i]+0.5*h*k51)
        k52 = S5(f[i]+0.5*h*k11, fpp[i]+0.5*h*k31, g[i]+0.5*h*k41, gp[i]+0.5*h*k51,Sm,Tinf,Minf,Pr,gamma)
        
        k13 = S1(fp[i]+0.5*h*k22)
        k23 = S2(fpp[i]+0.5*h*k32)
        k33 = S3(f[i]+0.5*h*k12, fpp[i]+0.5*h*k32, g[i]+0.5*h*k42, gp[i]+0.5*h*k52,Sm,Tinf)
        k43 = S4(gp[i]+0.5*h*k52)
        k53 = S5(f[i]+0.5*h*k12, fpp[i]+0.5*h*k32, g[i]+0.5*h*k42, gp[i]+0.5*h*k52,Sm,Tinf,Minf,Pr,gamma)
        
        k14 = S1(fp[i]+h*k23)
        k24 = S2(fpp[i]+h*k33)
        k34 = S3(f[i]+h*k13, fpp[i]+h*k33, g[i]+h*k43, gp[i]+h*k53,Sm,Tinf)
        k44 = S4(gp[i]+h*k53)
        k54 = S5(f[i]+h*k13, fpp[i]+h*k33, g[i]+h*k43, gp[i]+h*k53,Sm,Tinf,Minf,Pr,gamma)
        
        gp[i+1] = gp[i] + (1/6)*(k51 + 2*k52 + 2*k53 + k54)*h
        g[i+1] = g[i] + (1/6)*(k41 + 2*k42 + 2*k43 + k44)*h
        fpp[i+1] = fpp[i] + (1/6)*(k31 + 2*k32 + 2*k33 + k34)*h 
        fp[i+1] = fp[i] + (1/6)*(k21 + 2*k22 + 2*k23 + k24)*h
        f[i+1] = f[i] + (1/6)*(k11 + 2*k12 + 2*k13 + k14)*h

    return f,fp,fpp,g,gp

def sphere_solve(Minf=1.0, Tinf=300, TwTe=2.0, etamax=10, adi=1, Pr=0.72, gamma=1.4, N=100, precision=1e-9):
    """
    Self-similar compressible laminar solution at the stagnation point
    for a cylinder for either adiabatic or specifed wall temperature.
    Does not consider entropy layer or viscous interactions, and uses a
    constant pressure normal to the wall. Uses the Lees-Dorodnitsyn
    transformation

    Parameters
    ----------
    Minf : mach number at the edge of the boundary layer
    Tinf : temperature at the edge of the boundary layer
    TwTe : ratio of wall temperature to boundary layer edge temperature,
    only necessary if wall is not adiabatic
    adi : specify 1 for an adiabatic wall, 0 if wall is not adiabatic
    gamma : ratio of specific heats
    PR : Prandtl number
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
    f = np.zeros(N+1) # f
    fp = np.zeros(N+1) # f'
    fpp = np.zeros(N+1) # f''
    g = np.zeros(N+1) # ρ(eta)
    gp = np.zeros(N+1) # ρ(eta)'
    eta  = [i*deltaeta for i in range(N+1)]

    if adi == 1:
        f[0] = 0
        fp[0] = 0
        gp[0] = 0
        alpha0 = 0.1 # Initial Guess
        beta0 = 2.3118 # Initial Guess

    elif adi == 0:
        f[0] = 0
        fp[0] = 0
        g[0] = TwTe
        alpha0 = 0.1 # Initial Guess
        beta0 = 3.0 # Initial Guess

    maxiternum = 500
    y2n = 1
    y4n = 1
    for i in range(maxiternum):
        # Solve one
        if adi == 1:
            f[0] = 0
            fp[0] = 0
            gp[0] = 0
            fpp[0] = alpha0 # Initial Guess
            g[0] = beta0 # Initial Guess
        elif adi == 0:
            f[0] = 0
            fp[0] = 0
            g[0] = TwTe

            fpp[0] = alpha0 # Initial Guess
            gp[0] = beta0 # Initial Guess
    
		# First solution for Newton's iteration
        f,fp,fpp,g,gp = rungekutta_cylinder(N,deltaeta,f,fp,fpp,g,gp,Sm,Tinf,Pr,gamma,Minf)

        # Storing the freestream values for Newton's iteration method
        y20 = fp[N]
        y40 = g[N]
        
        # Small number addition for Newton's iteration method
        if adi == 1:
            f[0] = 0
            fp[0] = 0
            gp[0] = 0

            fpp[0] = alpha0+delta # Initial Guess + Small number
            g[0] = beta0 # Initial Guess
        elif adi == 0:
            f[0] = 0
            fp[0] = 0
            g[0] = TwTe

            fpp[0] = alpha0+delta # Initial Guess + Small number
            gp[0] = beta0 # Initial Guess
    
		# Second solution for Newton's iteration
        f,fp,fpp,g,gp = rungekutta_comp(N,deltaeta,f,fp,fpp,g,gp,Sm,Tinf,Pr,gamma,Minf)

        # Storing the freestream values for Newton's iteration method
        y201 = fp[N]
        y401 = g[N]

        # Small number addition for Newton's iteration method
        if adi == 1:
            f[0] = 0
            fp[0] = 0
            gp[0] = 0
            fpp[0] = alpha0 # Initial Guess
            g[0] = beta0+delta # Initial Guess + Small number
        elif adi == 0:
            f[0] = 0
            fp[0] = 0
            g[0] = TwTe
            fpp[0] = alpha0 # Initial Guess
            gp[0] = beta0+delta # Initial Guess + Small number
    
		# Third solution for Newton's iteration
        f,fp,fpp,g,gp = rungekutta_comp(N,deltaeta,f,fp,fpp,g,gp,Sm,Tinf,Pr,gamma,Minf)

        # Storing the freestream values for Newton's iteration method
        y202 = fp[N]
        y402 = g[N]

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

        if ((np.abs(fp[N]-1.0)<precision) and 
            (np.abs(g[N]-1.0)<precision) and 
            (np.abs(y2n-np.linalg.norm(fp))<precision) and 
            (np.abs(y4n-np.linalg.norm(g))<precision)):
            break
    
        y2n = np.linalg.norm(fp)
        y4n = np.linalg.norm(g)

	# Copying values for logical names
    U = fp
    T = g

    # Integration for eta --> y transformation
    y = np.zeros(N+1)
    for i in range(N):
        y[i] = y[i-1] + g[i]*(eta[i]-eta[i-1])

    y = y*np.sqrt(2)

    return eta,y,U,T,N

def fcyl(fp):
    """Add in definitions"""
    return fp

def fpcyl(fpp):
    """Add in definitions"""
    return fpp

def fppcyl(f,fp,fpp,g,gp,Sm,Tinf):
    """Add in definitions"""
    
    return -fpp*(gp/(2*g)-gp/(g+Sm/Tinf))-f*fpp*((g+Sm/Tinf)/(np.sqrt(g)*(1+Sm/Tinf)))+(fp**2-g)*((g+Sm/Tinf)/(np.sqrt(g)*(1+Sm/Tinf)))
    # return -fpp*(gp/(2*g)-gp/(g+Sm/Tinf))-f*fpp*((g+Sm/Tinf)/(np.sqrt(g)*(1+Sm/Tinf)))

def gcyl(gp):
    """Add in definitions"""
    return gp

def gpcyl(f,fpp,g,gp,Sm,Tinf,Minf,Pr,gamma):
    """Add in definitions"""
    # return -gp**2*(0.5/g-1/(g+Sm/Tinf))-Pr*f*gp/np.sqrt(g)*(g+Sm/Tinf)/(1+Sm/Tinf)
    return -gp**2*(0.5/g-1/(g+Sm/Tinf))-Pr*f*gp/np.sqrt(g)*(g+Sm/Tinf)/(1+Sm/Tinf)-(gamma-1)*Pr*Minf**2*fpp**2

def rungekutta_cylinder(N,h,f,fp,fpp,g,gp,Sm,Tinf,Pr,gamma,Minf):
    """Fourth order Runge-Kutta method specifically for compressible
    boundary layer flow functions"""
    # print(f,fp,fpp,g,gp,Sm,Tinf,Pr)
    for i in range(N):
        k11 = fcyl(fp[i])
        k21 = fpcyl(fpp[i])
        k31 = fppcyl(f[i], fp[i], fpp[i], g[i], gp[i],Sm,Tinf)
        k41 = gcyl(gp[i])
        k51 = gpcyl(f[i], fpp[i], g[i], gp[i],Sm,Tinf,Minf,Pr,gamma)  
        # print(k31, k51)
        k12 = fcyl(fp[i]+0.5*h*k21)
        k22 = fpcyl(fpp[i]+0.5*h*k31)
        k32 = fppcyl(f[i]+0.5*h*k11, fp[i]+0.5*h*k31, fpp[i]+0.5*h*k31, g[i]+0.5*h*k41, gp[i]+0.5*h*k51,Sm,Tinf)
        k42 = gcyl(gp[i]+0.5*h*k51)
        k52 = gpcyl(f[i]+0.5*h*k11, fpp[i]+0.5*h*k31, g[i]+0.5*h*k41, gp[i]+0.5*h*k51,Sm,Tinf,Minf,Pr,gamma)
        # print(k32, k52)
        k13 = fcyl(fp[i]+0.5*h*k22)
        k23 = fpcyl(fpp[i]+0.5*h*k32)
        k33 = fppcyl(f[i]+0.5*h*k12, fp[i]+0.5*h*k32, fpp[i]+0.5*h*k32, g[i]+0.5*h*k42, gp[i]+0.5*h*k52,Sm,Tinf)
        k43 = gcyl(gp[i]+0.5*h*k52)
        k53 = gpcyl(f[i]+0.5*h*k12, fpp[i]+0.5*h*k32, g[i]+0.5*h*k42, gp[i]+0.5*h*k52,Sm,Tinf,Minf,Pr,gamma)
        # print(k33, k53)
        k14 = fcyl(fp[i]+h*k23)
        k24 = fpcyl(fpp[i]+h*k33)
        k34 = fppcyl(f[i]+h*k13, fp[i]+0.5*h*k33, fpp[i]+h*k33, g[i]+h*k43, gp[i]+h*k53,Sm,Tinf)
        k44 = gcyl(gp[i]+h*k53)
        k54 = gpcyl(f[i]+h*k13, fpp[i]+h*k33, g[i]+h*k43, gp[i]+h*k53,Sm,Tinf,Minf,Pr,gamma)
        # print(k34, k54)
        gp[i+1] = gp[i] + (1/6)*(k51 + 2*k52 + 2*k53 + k54)*h
        g[i+1] = g[i] + (1/6)*(k41 + 2*k42 + 2*k43 + k44)*h
        fpp[i+1] = fpp[i] + (1/6)*(k31 + 2*k32 + 2*k33 + k34)*h 
        fp[i+1] = fp[i] + (1/6)*(k21 + 2*k22 + 2*k23 + k24)*h
        f[i+1] = f[i] + (1/6)*(k11 + 2*k12 + 2*k13 + k14)*h
    
    return f,fp,fpp,g,gp

def cylinder_solve_ig(Minf=2, Tinf=300, TwTe=2.0, etamax=10, adi=1, Pr=0.72, gamma=1.4, N=100, precision=1e-9):
    """
    Self-similar compressible laminar solution at the stagnation point
    of a sphere for either adiabatic or specifed wall temperature. Does
    not consider entropy layer or viscous interactions, and uses a
    constant pressure normal to the wall. Uses the Lees-Dorodnitsyn
    transformation

    Parameters
    ----------
    Minf : mach number at the edge of the boundary layer
    Tinf : temperature at the edge of the boundary layer
    TwTe : ratio of wall temperature to boundary layer edge temperature,
    only necessary if wall is not adiabatic
    adi : specify 1 for an adiabatic wall, 0 if wall is not adiabatic
    gamma : ratio of specific heats
    PR : Prandtl number
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
    f = np.zeros(N+1) # f
    fp = np.zeros(N+1) # f'
    fpp = np.zeros(N+1) # f''
    g = np.zeros(N+1) # ρ(eta)
    gp = np.zeros(N+1) # ρ(eta)'
    eta  = [i*deltaeta for i in range(N+1)]
    alpha = 1
    beta  = 1
    if adi == 1:
        f[0] = 0
        fp[0] = 0
        gp[0] = 0
        alpha0 = alpha # Initial Guess
        beta0 = beta # Initial Guess

    elif adi == 0:
        f[0] = 0
        fp[0] = 0
        g[0] = TwTe
        alpha0 = alpha # Initial Guess
        beta0 = beta # Initial Guess

    maxiternum = 500
    y2n = 1
    y4n = 1
    for i in range(maxiternum):
        # Solve one
        if adi == 1:
            f[0] = 0
            fp[0] = 0
            gp[0] = 0
            fpp[0] = alpha0 # Initial Guess
            g[0] = beta0 # Initial Guess
        elif adi == 0:
            f[0] = 0
            fp[0] = 0
            g[0] = TwTe
            fpp[0] = alpha0 # Initial Guess
            gp[0] = beta0 # Initial Guess
    
		# First solution for Newton's iteration
        f,fp,fpp,g,gp = rungekutta_cylinder(N,deltaeta,f,fp,fpp,g,gp,Sm,Tinf,Pr,gamma,Minf)
        # print(np.count_nonzero(np.isnan(f)),np.count_nonzero(np.isnan(fp)),np.count_nonzero(np.isnan(fpp)),np.count_nonzero(np.isnan(g)),np.count_nonzero(np.isnan(gp)))
        # plt.plot(fp, eta, label='1fp', color='blue')
        # plt.plot(f, eta, label='1g', color='blue', linestyle='--')
        # Storing the freestream values for Newton's iteration method
        y20 = np.count_nonzero(np.isnan(fp))
        y40 = np.count_nonzero(np.isnan(g))
        
        # Small number addition for Newton's iteration method
        if adi == 1:
            f[0] = 0
            fp[0] = 0
            gp[0] = 0

            fpp[0] = alpha0+.05 # Initial Guess + Small number
            g[0] = beta0 # Initial Guess
        elif adi == 0:
            f[0] = 0
            fp[0] = 0
            g[0] = TwTe

            fpp[0] = alpha0+.05 # Initial Guess + Small number
            gp[0] = beta0 # Initial Guess
    
		# Second solution for Newton's iteration
        f,fp,fpp,g,gp = rungekutta_cylinder(N,deltaeta,f,fp,fpp,g,gp,Sm,Tinf,Pr,gamma,Minf)
        # print(np.count_nonzero(np.isnan(f)),np.count_nonzero(np.isnan(fp)),np.count_nonzero(np.isnan(fpp)),np.count_nonzero(np.isnan(g)),np.count_nonzero(np.isnan(gp)))
        # plt.plot(fp, eta, label='2fp', color='red')
        # plt.plot(g, eta, label='2g', color='red', linestyle='--')

        # Storing the freestream values for Newton's iteration method
        y201 = np.count_nonzero(np.isnan(fp))
        y401 = np.count_nonzero(np.isnan(g))

        # Small number addition for Newton's iteration method
        if adi == 1:
            f[0] = 0
            fp[0] = 0
            gp[0] = 0
            fpp[0] = alpha0 # Initial Guess
            g[0] = beta0+.05 # Initial Guess + Small number
        elif adi == 0:
            f[0] = 0
            fp[0] = 0
            g[0] = TwTe
            fpp[0] = alpha0 # Initial Guess
            gp[0] = beta0+.05 # Initial Guess + Small number
    
		# Third solution for Newton's iteration
        f,fp,fpp,g,gp = rungekutta_cylinder(N,deltaeta,f,fp,fpp,g,gp,Sm,Tinf,Pr,gamma,Minf)
        # print(np.count_nonzero(np.isnan(f)),np.count_nonzero(np.isnan(fp)),np.count_nonzero(np.isnan(fpp)),np.count_nonzero(np.isnan(g)),np.count_nonzero(np.isnan(gp)))

        # plt.plot(fp, eta, label='3fp', color='green')
        # plt.plot(g, eta, label='3g', color='green', linestyle='--')
        # plt.legend()
        # plt.title(i)
        # plt.show()
        # Storing the freestream values for Newton's iteration method
        y202 = np.count_nonzero(np.isnan(fp))
        y402 = np.count_nonzero(np.isnan(g))
        
        if len(fp) == len(fp)-y20:
            break
        elif len(fp) == len(fp)-y201:
            alpha0 = alpha0+.05
            break
        elif len(fp) == len(fp)-y202:
            beta0 = beta0+.05
            break
        else:
            print(i, 'old', alpha0, beta0)
            if y201 > y20: # if number of NaNs increases when increasing alpha
                alpha0 -= (0.5 * abs((y20-len(fp)))/len(fp))
                print('decrease alpha fp', alpha0)
            else: # if number of NaNs decreases when increasing alpha
                alpha0 += (0.5 * abs((y20-len(fp)))/len(fp))
                print('increase alpha fp', alpha0)

            if y401 > y40: # if number of NaNs increases when increasing alpha
                alpha0 -= (0.5 * abs((y40-len(fp)))/len(fp))
                print('decrease alpha g', alpha0)

            else:  # if number of NaNs decreases when increasing alpha
                alpha0 += (0.5 * abs((y40-len(fp)))/len(fp))
                print('increase alpha g', alpha0)


            # if y202 > y20: # if number of NaNs increases when increasing beta
            #     beta0 -= (0.5 * abs((y20-len(fp)))/len(fp))
            #     print('decrease beta fp', alpha0)

            # else: # if number of NaNs decreases when increasing beta
            #     beta0 += (0.5 * abs((y20-len(fp)))/len(fp))
            #     print('increase beta fp', alpha0)


            # if y402 > y40: # if number of NaNs increases when increasing beta
                
            #     beta0 -= (0.5 * abs((y40-len(fp)))/len(fp))
            #     print('decrease beta g', alpha0)
            # else:  # if number of NaNs decreases when increasing beta
                
            #     beta0 += (0.5 * abs((y40-len(fp)))/len(fp))
            #     print('increase beta g', alpha0)

            print(i, 'new', alpha0, beta0)
        
        if i == 10:
            break

    return alpha0, beta0

def cylinder_solve(Minf=2, Tinf=300, TwTe=2.0, etamax=10, adi=1, Pr=0.72, gamma=1.4, N=200, precision=1e-6, alpha=1, beta=1):
    """
    Self-similar compressible laminar solution at the stagnation point
    of a sphere for either adiabatic or specifed wall temperature. Does
    not consider entropy layer or viscous interactions, and uses a
    constant pressure normal to the wall. Uses the Lees-Dorodnitsyn
    transformation

    Parameters
    ----------
    Minf : mach number at the edge of the boundary layer
    Tinf : temperature at the edge of the boundary layer
    TwTe : ratio of wall temperature to boundary layer edge temperature,
    only necessary if wall is not adiabatic
    adi : specify 1 for an adiabatic wall, 0 if wall is not adiabatic
    gamma : ratio of specific heats
    PR : Prandtl number
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
    delta = precision/100 # Small number for shooting method

    # Initializing the solution vectors
    f = np.zeros(N+1) # f
    fp = np.zeros(N+1) # f'
    fpp = np.zeros(N+1) # f''
    g = np.zeros(N+1) # ρ(eta)
    gp = np.zeros(N+1) # ρ(eta)'
    eta  = [i*deltaeta for i in range(N+1)]
    # alpha = alpha
    # beta = alpha
    if adi == 1:
        f[0] = 0
        fp[0] = 0
        gp[0] = 0
        alpha0 = alpha # Initial Guess
        beta0 = beta # Initial Guess

    elif adi == 0:
        f[0] = 0
        fp[0] = 0
        g[0] = TwTe
        alpha0 = alpha # Initial Guess
        beta0 = beta # Initial Guess

    maxiternum = 500
    y2n = 1
    y4n = 1
    for iter in range(maxiternum):
        # Solve one
        if adi == 1:
            f[0] = 0
            fp[0] = 0
            gp[0] = 0
            fpp[0] = alpha0 # Initial Guess
            g[0] = beta0 # Initial Guess
        elif adi == 0:
            f[0] = 0
            fp[0] = 0
            g[0] = TwTe

            fpp[0] = alpha0 # Initial Guess
            gp[0] = beta0 # Initial Guess
    
		# First solution for Newton's iteration
        f,fp,fpp,g,gp = rungekutta_cylinder(N,deltaeta,f,fp,fpp,g,gp,Sm,Tinf,Pr,gamma,Minf)
        # print(i, np.count_nonzero(np.isnan(f)),np.count_nonzero(np.isnan(fp)),np.count_nonzero(np.isnan(fpp)),np.count_nonzero(np.isnan(g)),np.count_nonzero(np.isnan(gp)))

        # Storing the freestream values for Newton's iteration method
        y20 = fp[N]
        y40 = g[N]
        
        # Small number addition for Newton's iteration method
        if adi == 1:
            f[0] = 0
            fp[0] = 0
            gp[0] = 0

            fpp[0] = alpha0+delta # Initial Guess + Small number
            g[0] = beta0 # Initial Guess
        elif adi == 0:
            f[0] = 0
            fp[0] = 0
            g[0] = TwTe

            fpp[0] = alpha0+delta # Initial Guess + Small number
            gp[0] = beta0 # Initial Guess
    
		# Second solution for Newton's iteration
        f,fp,fpp,g,gp = rungekutta_cylinder(N,deltaeta,f,fp,fpp,g,gp,Sm,Tinf,Pr,gamma,Minf)
        # print(i, np.count_nonzero(np.isnan(f)),np.count_nonzero(np.isnan(fp)),np.count_nonzero(np.isnan(fpp)),np.count_nonzero(np.isnan(g)),np.count_nonzero(np.isnan(gp)))

        # Storing the freestream values for Newton's iteration method
        y201 = fp[N]
        y401 = g[N]

        # Small number addition for Newton's iteration method
        if adi == 1:
            f[0] = 0
            fp[0] = 0
            gp[0] = 0
            fpp[0] = alpha0 # Initial Guess
            g[0] = beta0+delta # Initial Guess + Small number
        elif adi == 0:
            f[0] = 0
            fp[0] = 0
            g[0] = TwTe
            fpp[0] = alpha0 # Initial Guess
            gp[0] = beta0+delta # Initial Guess + Small number
    
		# Third solution for Newton's iteration
        f,fp,fpp,g,gp = rungekutta_cylinder(N,deltaeta,f,fp,fpp,g,gp,Sm,Tinf,Pr,gamma,Minf)
        # print(i, np.count_nonzero(np.isnan(f)),np.count_nonzero(np.isnan(fp)),np.count_nonzero(np.isnan(fpp)),np.count_nonzero(np.isnan(g)),np.count_nonzero(np.isnan(gp)))

        # Storing the freestream values for Newton's iteration method
        y202 = fp[N]
        y402 = g[N]

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

        # print(iter, np.abs(fp[N]-1.0), np.abs(g[N]-1.0), np.abs(y2n-np.linalg.norm(fp)), np.abs(y4n-np.linalg.norm(g)), precision)

        if ((np.abs(fp[N]-1.0)<precision) and 
            (np.abs(g[N]-1.0)<precision) and 
            (np.abs(y2n-np.linalg.norm(fp))<precision) and 
            (np.abs(y4n-np.linalg.norm(g))<precision)):
            break
        
        y2n = np.linalg.norm(fp)
        y4n = np.linalg.norm(g)
        if np.isnan(y2n):
            iter = np.NaN
            break
    
	# Copying values for logical names
    U = fp
    T = g

    # Integration for eta --> y transformation
    y = np.zeros(N+1)
    for i in range(N):
        y[i] = y[i-1] + g[i]*(eta[i]-eta[i-1])

    y = y*np.sqrt(2)
    if iter == 499:
        iter = np.NaN
    print(iter)
    return eta,y,U,T,N,iter

def compfpeqs(X, Sm, Tinf, Pr, gamma, Minf):
    """
    Parameters
    ----------
    f  = X[0]
    f' = X[1]
    f" = X[2]
    g  = X[3]
    g' = X[4]

    Returns
    ----------
    Y[0] = f' 
    Y[1] = f" 
    Y[2] = f"'
    Y[3] = g' 
    Y[4] = g" 

    """
    return np.vstack((
        [X[1]],
        [X[2]],
        [-X[2]*(X[4]/(2*X[3])-X[4]/(X[3]+Sm/Tinf))-X[0]*X[2]*((X[3]+Sm/Tinf)/(np.sqrt(X[3])*(1+Sm/Tinf)))],
        [X[4]],
        [-X[4]**2*(0.5/X[3]-1/(X[3]+Sm/Tinf))-Pr*X[0]*X[4]/np.sqrt(X[3])*(X[3]+Sm/Tinf)/(1+Sm/Tinf)-(gamma-1)*Pr*Minf**2*X[2]**2]
        ))

def rungekutta_comp(X,N,deltaeta,Sm,Tinf,Pr,gamma,Minf):
    """Fourth order Runge-Kutta method specifically for compressible
    boundary layer flow functions"""
    for i in range(N):
        k1s = compfpeqs(X[:, [i]], Sm, Tinf, Pr, gamma, Minf)
        k2s = compfpeqs(X[:, [i]]+deltaeta/2*k1s, Sm, Tinf, Pr, gamma, Minf)
        k3s = compfpeqs(X[:, [i]]+deltaeta/2*k2s, Sm, Tinf, Pr, gamma, Minf)
        k4s = compfpeqs(X[:, [i]]+deltaeta*k3s, Sm, Tinf, Pr, gamma, Minf)
        X[:, [i+1]] = X[:, [i]] + (1/6)*(k1s + 2*k2s + 2*k3s + k4s)*deltaeta
    return X

def compressible_solve(Minf=2, Tinf=300, TwTe=2.0, etamax=10, adi=1, Pr=0.72, gamma=1.4, N=100, precision=1e-9):
    """
    Self-Similar compressible laminar solution over flat plate for
    either adiabatic or specifed wall temperature. Does not consider
    entropy layer or viscous interactions, and uses a constant pressure
    normal to the wall. Uses the Lees-Dorodnitsyn transformation
    Source:
    https://github.com/frkanz/A-CFD-Tutorial-in-Julia-Compressible-Blasius

    Parameters
    ----------
    Minf : mach number at the edge of the boundary layer
    Tinf : temperature at the edge of the boundary layer
    TwTe : ratio of wall temperature to boundary layer edge temperature,
    only necessary if wall is not adiabatic
    adi : specify 1 for an adiabatic wall, 0 if wall is not adiabatic
    gamma : ratio of specific heats
    PR : Prandtl number
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

    # Initializing the solution vectors
    X = np.zeros((5, N+1))
    # The rows of X are a variable at each eta
    # f   = X[0, :]
    # f'  = X[1, :]
    # f"  = X[2, :]
    # g   = X[3, :]
    # g'  = X[4, :]
    
    eta  = [i*deltaeta for i in range(N+1)]

    if adi == 1:
        X[0, 0] = 0
        X[1, 0] = 0
        X[4, 0] = 0
        alpha0 = 0.1 # Initial Guess
        beta0 = 2.3118 # Initial Guess

    elif adi == 0:
        X[0, 0] = 0
        X[1, 0] = 0
        X[4, 0] = TwTe
        alpha0 = 0.1 # Initial Guess
        beta0 = 3.0 # Initial Guess
    
    # Solve one
    if adi == 1:
        X[0, 0] = 0
        X[1, 0] = 0
        X[4, 0] = 0
        X[2, 0] = alpha0 # Initial Guess
        X[3, 0] = beta0 # Initial Guess
    elif adi == 0:
        X[0, 0] = 0
        X[1, 0] = 0
        X[4, 0] = TwTe

        X[2, 0] = alpha0 # Initial Guess
        X[3, 0] = beta0 # Initial Guess

    xidx = [2, 3]
    fxidx = [1, 3]
    broyden2(X, rungekutta_comp, 
             goal=1, xidx=xidx, fxidx=fxidx, N=N, 
             Tinf=Tinf,Pr=Pr,gamma=gamma,Minf=Minf,deltaeta=deltaeta,Sm=Sm)

    # Copying values for logical names
    U = X[1, :]
    T = X[3, :]

    # Integration for eta --> y transformation
    y = np.zeros(N+1)
    for i in range(N):
        y[i] = y[i-1] + X[3, :][i]*(eta[i]-eta[i-1])

    y = y*np.sqrt(2)

    return eta,y,U,T,N

def incompfpeqs(X):
    """
    Parameters
    ----------
    f  = X[0]
    w1 = X[1]
    w2 = X[2]
    
    """
    return np.vstack((
        [X[1]],
        [X[2]],
        [-X[0]*X[2]/2]
        ))

def rungekutta_incomp(X, N, step):
    """Fourth order Runge-Kutta method specifically for incompressible
    boundary layer flow functions"""
    for i in range(N):
        k1s = incompfpeqs(X[:, [i]])
        k2s = incompfpeqs(X[:, [i]]+step/2*k1s)
        k3s = incompfpeqs(X[:, [i]]+step/2*k2s)
        k4s = incompfpeqs(X[:, [i]]+step*k3s)
        
        X[:, [i+1]] = X[:, [i]] + (1/6)*(k1s + 2*k2s + 2*k3s + k4s)*step

    return X

def incompressible_solve(step=0.1, precision=1e-6):
    """Self-similar incompressible laminar boundary layer on a flat 
    plate using the Blasius transformation. No parameters required.

    Parameters
    ----------
    step : how refined the vertical direction is
    precision : accuracy for convergence

    Returns
    -------
    eta : the non-dimenional height of the boundary layer
        : ysqrt(uinf/nu/x)
    X : the 3 rows of X are described as follows
    X[0, :] = f : related to the stream function
    X[1, :] = w1 : the velocity in the boundary normalized to the freestream 
        : u/uinf
    X[2, :] = w2 : v*sqrt(x/nu/uinf)

    (nu : the kinematic viscosity)
    """
    stop = 10
    N = int(stop/step)

    X = np.zeros((3, N+1))
    X[2, 0] = .2 # initial guess
    # The rows of X are a variable at each eta
    # f   = X[0, :]
    # w1  = X[1, :]
    # w2  = X[2, :]

    broyden2(X, rungekutta_incomp, 
             goal=1, xidx=[2], fxidx=[1], N=N, step=step, precision=precision)
    eta = np.linspace(0, stop, N+1)
    return eta, X

def mstmp(X):
    """
    Parameters
    ----------
    Y = X[0]

    Returns
    ----------
    X[end] = Y"
    """
    return np.vstack((
        [X[1]],
        [X[2]],
        [5*np.sinh(5*X[0])]
        ))

def rk_ms(X,N,deltaeta):
    """Fourth order Runge-Kutta method specifically for multiple 
    shotting method development"""
    for i in range(N):
        k1s = mstmp(X[:, [i]])
        k2s = mstmp(X[:, [i]]+deltaeta/2*k1s)
        k3s = mstmp(X[:, [i]]+deltaeta/2*k2s)
        k4s = mstmp(X[:, [i]]+deltaeta*k3s)
        
        X[:, [i+1]] = X[:, [i]] + (1/6)*(k1s + 2*k2s + 2*k3s + k4s)*deltaeta
    return X

def multiple_shooting(syseqs, numsubd:int=10, numsteps:int=50):
    """
    Parameters
    ----------
    syseqs : a function that contains the system of equations for
        solving, output must be a column vector
    numsubd : number of subdivisions spanning between the integration 
        bounds, must be at least 2
    numsteps : number of points that will span one subdivision
    
    
    Returns
    -------
    The final integrated solution
    """
    
    numuk = 2 # how many variables are being guessed per subdivision
    numtotuk = numuk*numsubd-1 # total number of unknowns

    # knowns
    y0 = 0
    yf = 1
    x0 = 0
    xf = .5
    xsubd = np.linspace(x0, xf, numsubd+1)
    # initial guesses
    svar1 = xsubd[0:-1]/(xf-x0)*(yf-y0) # y values
    svar2 = np.zeros_like(svar1) # first derivative
    svar = np.vstack((svar1, svar2))
    # flatten in way the variables are ordered by increasing shoot
    # number and increasing derivative
    svar = svar.flatten('F')

    for idx in range(numsubd):
        # First shoot
        # set the initial guesses for the shoot
        Xfirst = np.zeros((3, numsteps+1))
        Xfirst[0:numuk, 0] = svar[0:numuk]
        # set the step distance for the shoot
        stepdist = (xsubd[1] - xsubd[0])/numsteps
        # calculate the result of the shoot
        syseqs(Xfirst, numsteps, stepdist)
        Ffirsts = Xfirst[0:numuk, -1] - svar[numuk:2*numuk]

        # calculate the fintie difference matrix for the guessed
        # variables
        fdJfirst = fdJ(Xfirst, syseqs, xidx=[1], fxidx=[0, 1], N=numsteps, deltaeta=stepdist)
        FDJfirstrow = np.hstack((fdJfirst, -np.eye(numuk), np.zeros((numuk, numtotuk-numuk-1))))
        
        # Middle shoots
        # same pattern as above
        FDJmidrows = np.zeros((numtotuk-numuk-1, numtotuk))
        Fmids = np.zeros(numtotuk-numuk-1)
        for mdx in range(numsubd-2):
            Xmid = np.zeros((3, numsteps+1))
            Xmid[0:numuk, 0] = svar[numuk*(mdx+1):numuk*(mdx+1)+2]
            stepdist = (xsubd[mdx+2] - xsubd[mdx+1])/numsteps

            syseqs(Xmid, numsteps, stepdist)
            Fmids[mdx*numuk:(1+mdx)*numuk] = Xmid[0:numuk, -1] - svar[(mdx+2)*numuk:(mdx+3)*numuk]
            
            fdJmid = fdJ(Xmid, syseqs, xidx=[0, 1], fxidx=[0, 1], N=numsteps, deltaeta=stepdist)
            FDJmidrows[mdx*numuk:(mdx+1)*numuk, mdx*numuk+1:(2+mdx)*numuk+1] = np.hstack((fdJmid, -np.eye(numuk)))

        # Last shoot
        Xlast = np.zeros((3, numsteps+1))
        Xlast[0:numuk, 0] = svar[-numuk:]
        stepdist = (xsubd[-1] - xsubd[-2])/numsteps
        syseqs(Xlast, numsteps, stepdist)
        Flasts = Xlast[0, -1] - yf
        
        fdJlast = fdJ(Xlast, syseqs, xidx=[0, 1], fxidx=[0], N=numsteps, deltaeta=stepdist)
        FDJlastrow = np.hstack((np.zeros(numtotuk-numuk), fdJlast[0, :]))
        
        # Function result matrix and finite difference matrix
        F = np.hstack(
            [Ffirsts,
             Fmids,
             Flasts])
        FDJ = np.vstack([FDJfirstrow,
                        FDJmidrows, 
                        FDJlastrow])
        invFDJ = np.linalg.inv(FDJ)
        supdate = np.matmul(invFDJ, F)
        svar[1:] = svar[1:] - supdate
    
    # Final calculation to return result to user
    X = np.zeros((3, numsubd*numsteps+1))
    X[0:numuk, 0] = svar[0:numuk]
    N = numsubd*numsteps
    deltaeta = stepdist
    syseqs(X,N,deltaeta)
    
    return np.linspace(x0, xf, numsubd*numsteps+1), X

def multiple_shooting2(syseqs, numsubd:int=10, numsteps:int=50):
    """
    Parameters
    ----------
    numsubd : number of subdivisions spanning between the integration 
        bounds, must be at least 2
    numsteps : number of points that will span one subdivision
    
    
    Returns
    -------
    The final integrated solution
    """
    
    numuk = 2 # how many variables are being guessed per subdivision
    numtotuk = numuk*numsubd-1 # total number of unknowns

    # knowns
    y0 = 0
    yf = 1
    x0 = 0
    xf = .5
    xsubd = np.linspace(x0, xf, numsubd+1)
    # initial guesses
    svar1 = xsubd[0:-1]/(xf-x0)*(yf-y0) # y values
    svar2 = np.zeros_like(svar1) # first derivative
    svar = np.vstack((svar1, svar2))
    # flatten in way the variables are ordered by increasing shoot
    # number and increasing derivative
    svar = svar.flatten('F')

    # First shoot
    # set the initiual guesses for the shoot
    Xfirst = np.zeros((3, numsteps+1))
    Xfirst[0:numuk, 0] = svar[0:numuk]
    # set the step distance for the shoot
    stepdist = (xsubd[1] - xsubd[0])/numsteps
    # calculate the result of the shoot
    syseqs(Xfirst, numsteps, stepdist)
    Ffirsts = Xfirst[0:numuk, -1] - svar[numuk:2*numuk]

    # calculate the fintie difference matrix for the guessed
    # variables
    fdJfirst = fdJ(Xfirst, syseqs, xidx=[1], fxidx=[0, 1], N=numsteps, deltaeta=stepdist)
    FDJfirstrow = np.hstack((fdJfirst, -np.eye(numuk), np.zeros((numuk, numtotuk-numuk-1))))
    
    # Middle shoots
    # same pattern as above
    FDJmidrows = np.zeros((numtotuk-numuk-1, numtotuk))
    Fmids = np.zeros(numtotuk-numuk-1)
    for mdx in range(numsubd-2):
        Xmid = np.zeros((3, numsteps+1))
        Xmid[0:numuk, 0] = svar[numuk*(mdx+1):numuk*(mdx+1)+2]
        stepdist = (xsubd[mdx+2] - xsubd[mdx+1])/numsteps

        syseqs(Xmid, numsteps, stepdist)
        Fmids[mdx*numuk:(1+mdx)*numuk] = Xmid[0:numuk, -1] - svar[(mdx+2)*numuk:(mdx+3)*numuk]
        
        fdJmid = fdJ(Xmid, syseqs, xidx=[0, 1], fxidx=[0, 1], N=numsteps, deltaeta=stepdist)
        FDJmidrows[mdx*numuk:(mdx+1)*numuk, mdx*numuk+1:(2+mdx)*numuk+1] = np.hstack((fdJmid, -np.eye(numuk)))

    # Last shoot
    Xlast = np.zeros((3, numsteps+1))
    Xlast[0:numuk, 0] = svar[-numuk:]
    stepdist = (xsubd[-1] - xsubd[-2])/numsteps
    syseqs(Xlast, numsteps, stepdist)
    Flasts = Xlast[0, -1] - yf
    
    fdJlast = fdJ(Xlast, syseqs, xidx=[0, 1], fxidx=[0], N=numsteps, deltaeta=stepdist)
    FDJlastrow = np.hstack((np.zeros(numtotuk-numuk), fdJlast[0, :]))
    
    # Function result matrix and finite difference matrix
    F = np.hstack(
        [Ffirsts,
            Fmids,
            Flasts])
    FDJ = np.vstack([FDJfirstrow,
                    FDJmidrows, 
                    FDJlastrow])
    # invFDJ = np.linalg.inv(FDJ)
    # supdate = np.matmul(invFDJ, F)
    # svar[1:] = svar[1:] - supdate
    
    # y2n = 1 # Can be used for judging convergence
    # y4n = 1 # Can be used for judging convergence
    # get the base value and numerical get the jacobian
    FX = F
    invJ = np.linalg.inv(FDJ)
    for iter in range(1):
        Xold_jac = svar
        X_jac = Xold_jac - np.matmul(invJ, FX)
        deltaX_jac = X_jac - Xold_jac
        svar = X_jac

        FXold_jac = FX
        FX = syseqs(X)
        deltaF_jac = FX[fxidx, [-1]] - FXold_jac
        trans = deltaX_jac[:,0] # tranpose
        # Next line is Broydens method, it is just a way to update the
        # jacobian without doing all of the perturbations
        invJ = invJ + (deltaX_jac - np.matmul(invJ,deltaF_jac))/np.matmul(np.matmul(trans,invJ), deltaF_jac)*np.matmul(trans,invJ)
        
        if np.all(np.abs(FX[fxidx, -1]-goal)<1e-6): #and
            # (np.abs(y2n-np.linalg.norm(FX[1, -1]))<precision) and 
            # (np.abs(y4n-np.linalg.norm(FX[3, -1]))<precision)):
            break
        # y2n = np.linalg.norm(FX[1, -1])
        # y4n = np.linalg.norm(FX[3, -1])
        
    return X, iter
    
    # Final calculation to return result to user
    X = np.zeros((3, numsubd*numsteps+1))
    X[0:numuk, 0] = svar[0:numuk]
    N = numsubd*numsteps
    deltaeta = stepdist
    syseqs(X,N,deltaeta)

    return np.linspace(x0, xf, numsubd*numsteps+1), X
    
def main():

    
    # plt.plot(h[0], h[3])
    # plt.show()
    x, y = multiple_shooting(rk_ms)
    plt.plot(x, y[0, :])
    plt.show()

if __name__ == '__main__':
    main()
