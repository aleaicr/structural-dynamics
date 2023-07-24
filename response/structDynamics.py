# Module Newmarkbeta method for solving the equation of motion of a SDOF system
# This function is used to solve the equation of motion of a SDOF system using the Newmarkbeta method.
# The method is an implicit method, which means that the solution is obtained by solving a system of equations.

# Comments:
# 

import numpy as np
from structFunctions import computeK


# Newmark method for linear linear SDOF systems
def newmark(upp, t_vect, beta, xi, u0, up0, T_vect):
    # This function returns the displacement, velocity and acceleration time series
    # of a SDOF system using the Newmark's method.

    # Inputs:
    # upp       --> Acceleration time serie (vector)
    # t_vect    --> Time vector of the sample (vector or scalar)
    #           --> If the time vector is a scalar is assumed that this is the time step of the acceleration time serie
    # beta      --> Newmark's method parameter
    # xi        --> Damping ratio (scalar)
    # u0        --> Initial displacement (scalar)
    # up0       --> Initial velocity (scalar)
    # T_vect    --> Periods to be evaluated (vector)

    # Outputs:
    # u         --> Displacement time serie (vector) relative to the ground
    # ud        --> Velocity time serie (vector) relative to the ground
    # udd       --> Acceleration time serie (vector) relative to the ground
    # udd_tot   --> Total (absolute) acceleration time serie (vector)
    
    # Other parameters
    gamma = 0.5             # gamma = 1/2 average acceleration method so that the method is unconditionally stable
    upp_length = len(upp)   # Length of the acceleration time serie
    T_length = len(T_vect)  # Length of the periods vector

    # Check time series
    if len(t_vect) == 1:
        dt = t_vect.copy()  # time step
        t_vect = np.arange(0, upp_length * dt, dt)  # time serie vector
    elif len(t_vect) != upp_length:
        raise ValueError('The length of the time vector must be the same as the acceleration time serie or '
                         'provide the time step instead')
    elif len(t_vect) == upp_length:
        pass

    # Initialize the response vectors
    u = np.zeros(upp_length)
    ud = np.zeros(upp_length)
    udd = np.zeros(upp_length)

    # Calculate the response for each period
    for j in range(T_length):
        # Initialize the response vectors
        u[:] = 0
        ud[:] = 0
        udd[:] = 0
        
        # Initial conditions
        u[0] = u0
        ud[0] = up0
        
        # If the period is zero, the response is same as the 
        if T_vect[j] == 0:
            udd = upp.copy()                                  
        else: 
            # Calculate the parameters of the Newmarkbeta method
            wn = 2 * np.pi / T_vect[j]
            udd[0] = upp[0] - 2 * xi * wn * up0 - wn ** 2 * up0
            dt = t_vect[1] - t_vect[0]
            a1 = 1 / (beta * dt ** 2) + 2 * xi * wn * gamma / (beta * dt)
            a2 = 1 / (beta * dt) + 2 * xi * wn * (gamma / beta - 1)
            a3 = (1 / (2 * beta) - 1) + 2 * xi * wn * dt * (gamma / (2 * beta) - 1)
            k_ton = a1 + wn ** 2

            # Calculate the response for each time step
            for i in range(1, upp_length):
                dt = t_vect[i] - t_vect[i - 1]
                p_ton = -upp[i] + a1 * u[i - 1] + a2 * ud[i - 1] + a3 * udd[i - 1]
                u[i] = p_ton / k_ton
                ud[i] = gamma / (beta * dt) * (u[i] - u[i - 1]) + (1 - gamma / beta) * ud[i - 1] + dt * (1 - gamma / (2 * beta)) * udd[i - 1]
                udd[i] = (u[i] - u[i - 1]) / (beta * dt ** 2) - ud[i - 1] / (beta * dt) - (1 / (2 * beta) - 1) * udd[i - 1]
            
        # Calculate the total acceleration (including the acceleration of the support)
        udd_tot = udd + upp

    return u, ud, udd, udd_tot


# Peak displacement and time using Newmark method for linear SDOF systems
def dispPeak(upp, t_vect, beta, xi, u0, up0, T_vect):
    # This function returns the peak displacement for each period 

    # Other parameters
    gamma = 0.5  # gamma = 1/2 average acceleration method so that the method is unconditionally stable
    upp_length = len(upp)  # Length of the acceleration time serie
    T_length = len(T_vect)  # Length of the periods vector
    t_vect = np.arange(0, upp_length * dt, dt)  # time serie vector

    # Initialize the response vectors
    u = np.zeros(upp_length)
    ud = np.zeros(upp_length)
    udd = np.zeros(upp_length)
    uPeak = np.zeros(T_length)
    tPeak = np.zeros(T_length)
    signPeak = np.zeros(T_length)
    Pos_tPeak = np.zeros(T_length)

    # Calculate the response for each period
    for j in range(T_length):
        # Initialize the response vectors
        u[:] = 0
        ud[:] = 0
        udd[:] = 0
        
        # Initial conditions
        u[0] = u0[j]
        ud[0] = up0[j]
        
        # If the period is zero, the response is same as the 
        if T_vect[j] == 0:
            udd = upp.copy()                                  
        else: 
            # Calculate the parameters of the Newmarkbeta method
            wn = 2 * np.pi / T_vect[j]
            udd[0] = upp[0] - 2 * xi * wn * up0 - wn ** 2 * up0
            dt = t_vect[1] - t_vect[0]
            a1 = 1 / (beta * dt ** 2) + 2 * xi * wn * gamma / (beta * dt)
            a2 = 1 / (beta * dt) + 2 * xi * wn * (gamma / beta - 1)
            a3 = (1 / (2 * beta) - 1) + 2 * xi * wn * dt * (gamma / (2 * beta) - 1)
            k_ton = a1 + wn ** 2

            # Calculate the response for each time step
            for i in range(1, upp_length):
                dt = t_vect[i] - t_vect[i - 1]
                p_ton = -upp[i] + a1 * u[i - 1] + a2 * ud[i - 1] + a3 * udd[i - 1]
                u[i] = p_ton / k_ton
                ud[i] = gamma / (beta * dt) * (u[i] - u[i - 1]) + (1 - gamma / beta) * ud[i - 1] + dt * (1 - gamma / (2 * beta)) * udd[i - 1]
                udd[i] = (u[i] - u[i - 1]) / (beta * dt ** 2) - ud[i - 1] / (beta * dt) - (1 / (2 * beta) - 1) * udd[i - 1]

        # Extra calculations if needed (these are for spectral matching) # 
        uPeak[j], Pos_tPeak[j] = np.max(np.abs(u)), np.argmax(np.abs(u))
        signPeak[j] = np.sign(u[Pos_tPeak[j]])
        tPeak[j] = t_vect[Pos_tPeak[j]]

    return uPeak, tPeak, signPeak, Pos_tPeak


# Newmark's method for non-linear SDOF systems
def newmarkBilineal(upp, t_vect, beta, xi, u0, up0, Tn, alpha, Fy, R_tol):
    # Newmark's method for non-linear SDOF systems returns the response of the system
    # and the resisting force and the ductility demand history
    # upp           --> Acceleration time serie
    # t_vect        --> Time vector
    # beta          --> Beta parameter of the Newmark's method
    # xi            --> Damping ratio
    # u0            --> Initial displacement
    # up0           --> Initial velocity
    # Tn            --> Natural period of the structure
    # alpha         --> Post-yield stiffness ratio
    # Fy            --> Yield force
    # R_tol         --> Tolerance for the Newton-Raphson method

    # Output
    # u            --> History of relative displacements
    # up           --> History of relative velocities
    # upp          --> History of absolute accelerations
    # fs           --> History of resisting forces
    # mu           --> Ductility demand history

    # Parameters
    gamma = 0.5                 # Gamma de Métodos de Newmark
    dt = t_vect[1] - t_vect[0]  # Paso temporal
    t_length = len(t_vect)      # Amount of time steps
    max_j_counter = 10000       # Maximum amount of iterations for Newton-Raphson method

    # Dynamic properties
    wn = 2 * np.pi / Tn         # Frecuencia natural estructura
    k = wn ** 2                 # Equivalent structure's stiffness
    uy = Fy / k                 # Yield displacement
    # m = 1                      # Equivalent structure's mass
    # c = 2 * xi * wn            # Equivalent structure's damping

    # Initialize the response vectors
    u = np.zeros(t_length)
    du = np.zeros(t_length)
    ddu = np.zeros(t_length)
    yield_ = np.zeros(t_length)

    # Initial conditions
    u[0] = u0                   # Displacement
    du[0] = up0                 # Velocity
    if abs(u0) <= uy:
        ddu[0] = upp[0] - 2 * xi * wn * up0 - wn ** 2 * u0
    else:
        ddu[0] = upp[0] - 2 * xi * wn * up0 - Fy * np.sign(u0) 

    # Newmark's Method
    # Initial calculations
    fs = np.zeros(t_length)     # Resisting force initialisation
    fs[0] = k * u[0]            # Resisting force initial value
    kt = k                      # Initial tangent stiffness
    a1 = 1 / (beta * dt ** 2) + 2 * xi * wn * gamma / (beta * dt)
    a2 = 1 / (beta * dt) + 2 * xi * wn * (gamma / beta - 1)
    a3 = (1 / (2 * beta) - 1) + 2 * xi * wn * dt * (gamma / (2 * beta) - 1)

    # Calculations for each time instant 
    for i in range(t_length - 1):
        j_counter = 0           # Initial counter for the iterations NewtonRaphson
        u[i + 1] = u[i]         # Initial guess for the displacement
        fs[i + 1] = fs[i]       # Initial guess for the resisting force
        p_tongo = -upp[i + 1] + a1 * u[i] + a2 * du[i] + a3 * ddu[i]
        R_tongo = p_tongo - fs[i + 1] - a1 * u[i + 1]
        while abs(R_tongo) > R_tol:     # Iteration loop until the convergence criteria is met
            yield_[i + 1] = 0
            kt_tongo = kt + a1
            Du = R_tongo / kt_tongo
            u[i + 1] += Du              # Displacement estimation update
            fs[i + 1] += k * Du         # Resisting force estmation update

            # State determination
            if abs(fs[i + 1] - u[i + 1] * alpha * k) > Fy * (1 - alpha):    # Plastic state
                fs[i + 1] = Fy * np.sign(fs[i + 1]) * (1 - alpha) + u[i + 1] * alpha * k
                yield_[i + 1] = np.sign(fs[i + 1])
                kt = alpha * k
            else:                                                           # Elastic state
                kt = k
            
            # If there is no convergence
            if j_counter >= max_j_counter:
                print(f'No convergence in {i}-th time step (time: {i*dt} [s]) in iteration {j_counter}')
                print(f'Last displacement: {u[i + 1]}')
                print(f'Last resisting force: {fs[i + 1]}')
                break

            # New iteration
            R_tongo = p_tongo - fs[i + 1] - a1 * u[i + 1]
            j_counter += 1
        
        # Calculations for velocity and acceleration
        du[i + 1] = gamma / (beta * dt) * (u[i + 1] - u[i]) + (1 - gamma / beta) * du[i] + dt * (1 - gamma / (2 * beta)) * ddu[i]
        ddu[i + 1] = (u[i + 1] - u[i]) / (beta * dt ** 2) - du[i] / (beta * dt) - (1 / (2 * beta) - 1) * ddu[i]
    
    mu = np.max(np.abs(u)) / uy         # Ductility demand

    return u, du, ddu, fs, yield_, mu

# MDOF Shear Ibarra-Medina-Krawinkler seismic response
def MDOF_Shear_IMK_seismic(h, wi, Pi, k, Xi, ug, dt, do, Vo, Fy, a_s, dcdy, a_c, tol, g, Name, StrengthLimitCheck):
    N = len(wi)     # Number of stories
    MaxIter = 20
    
    # Obtain T and phi
    
    # Note: first index corresponds to 1st floor, and last index to roof.
    # M and K matrices
    M = np.diag(wi)/g         # Mass Matrix
    K = computeK(k)           # Stiffness Matrix
    
    # Eigenvalue analysis
    w2, phi = np.linalg.eig(np.linalg.inv(M).dot(K))
    w = np.sqrt(w2)          # Undamped frequencies
    idx = w.argsort()
    w = w[idx]
    T = 2 * np.pi / w        # Undamped periods
    
    # Sort vectors (modal shapes) and normalize them at roof: phi_roof = 1.0
    sphi = np.copy(phi)
    for i in range(N):
        sphi[:, i] = phi[:, idx[i]] / phi[-1, idx[i]]
    phi = sphi             # Normalized modal shapes
    
    # C matrix
    Mi = np.diag(phi.T.dot(M).dot(phi))
    if Xi.shape == Mi.shape:
        Ci = 2 * Mi * w * Xi
    else:
        Ci = 2 * Mi * w * Xi.T
    C = np.linalg.inv(phi).dot(np.diag(Ci)).dot(np.linalg.inv(phi))
    
    # Check stability of the method
    Np = len(ug)     # Length of the record
    time = np.arange(0, dt * Np, dt)
    
    # If the time step is too large, display a warning.
    if dt > T[0] / 30:
        print('Warning: The time step used (dt) for the ground motion "' + Name + '" is greater than T_1/30. This is not recommended for representing the response correctly.')
        dt_ = dt / np.ceil(30 * dt / T[0])
        time_ = np.arange(0, dt_ * Np, dt_)
        ug = np.interp(time_, time, ug)
        time = time_
        dt = dt_
        Np = len(ug)
    
    # Initial Calculations
    r = np.ones(N)            # Note that this assumes horizontal excitation
    P = -M.dot(r) * ug * g    # Equivalent external load
    
    dy = Fy / k               # Yielding displacement of each story
    dc = dcdy * dy            # Capping displacement of each story
    Fmax = (1 - a_s) * Fy + a_s * k * dc    # Positive Strength Limit
    Fmin = -Fmax             # Negative Strength Limit
    LimMax = np.zeros(N)
    LimMin = np.zeros(N)
    
    # Initialize vectors
    fs_st = np.zeros((N, Np))        # Story restoring force
    fs = np.copy(fs_st)              # Total floor restoring force
    fs[:, 0] = np.concatenate((np.diff(fs_st[:, 0]), [fs_st[-1, 0]]))
    Kt = np.copy(K)                   # Initial tangent stiffness
    kt = np.copy(k)
    
    u = np.zeros((N, Np))            # Relative displacement time history
    v = np.zeros((N, Np))            # Relative velocity time history
    a = np.zeros((N, Np))            # Relative acceleration time history
    F = np.zeros((N, Np))            # Story force time history
    F[:, 0] = np.copy(P)             # Store equivalent seismic load
    c = np.zeros((N, Np))            # Damping force
    c[:, 0] = C.dot(v[:, 0])         # Store damping forces
    
    # Check to avoid zero dividing
    dt_ = np.copy(dt)
    if np.min(T) * 0.02 < dt:
        dt_ = np.min(T) * 0.02
    
    # Initialize the parameters
    i = 0
    nf = np.copy(i)
    AbsD = np.zeros((N, Np))
    AbsDMax = np.zeros((N, Np))
    AbsDMax[:, 0] = np.copy(np.abs(u[:, 0]))
    AbsDMax[:, 1] = np.copy(AbsDMax[:, 0])
    Drift = np.copy(AbsD)
    
    # Iteration of Newmark
    for i in range(1, Np):
        for j in range(N):
            # Determine force limits
            if Drift[j, i - 1] >= dc[j]:
                LimMax[j] = Fmax[j]
                LimMin[j] = Fmin[j]
            else:
                LimMax[j] = np.interp(Drift[j, i - 1], [0, dc[j]], [Fmax[j], a_c * Fmax[j]])
                LimMin[j] = np.interp(Drift[j, i - 1], [0, dc[j]], [Fmin[j], a_c * Fmin[j]])
            
            # Calculate stiffness of each story at each time step
            if (fs_st[j, i - 1] >= LimMax[j] and P[j] > 0) or (fs_st[j, i - 1] <= LimMin[j] and P[j] < 0):
                kt[j] = 0
                Kt[j, j] = 0
            else:
                kt[j] = k
                Kt[j, j] = K[j, j]
            
            # Acceleration
            a[j, i] = ((1 - a_s) / (1 + a_s)) * (P[j] - fs_st[j, i - 1] - c[j, i - 1] - kt[j] * u[j, i - 1]) / M[j, j]
            
            # Velocity
            v[j, i] = v[j, i - 1] + (1 - a_s) * dt_ * a[j, i] + a_s * dt_ * a[j, i - 1]
            
            # Displacement
            u[j, i] = u[j, i - 1] + v[j, i - 1] * dt_ + (1 / (2 * xi[j])) * ((1 - a_c) * a[j, i - 1] + a_c * a[j, i]) * dt_**2
            
            # Restoring force
            fs_st[j, i] = kt[j] * (u[j, i] - u[j, i - 1])
            fs[:, i] = np.concatenate((np.diff(fs_st[:, i]), [fs_st[-1, i]]))
            
            # Damping force
            c[:, i] = C.dot(v[:, i])
            
            # Total force
            F[:, i] = P + fs[:, i] + c[:, i]
            
            # Absolute displacement and drift
            AbsD[j, i] = np.abs(u[j, i])
            AbsDMax[j, i] = max(AbsDMax[j, i - 1], AbsD[j, i])
            Drift[j, i] = AbsDMax[j, i] / h[j]
        
        # Check if strength limit state is reached
        if StrengthLimitCheck == 1:
            if np.max(Drift[:, i]) >= 0.75 * np.max(dc):
                nf = np.copy(i)
                break
        
    # Return results
    return u, v, a, F, fs, AbsD, Drift, nf
