## Extra functions for structural analysis

import numpy as np

# Compute stiffness matrix
def computeK(k):
    # Creata a stiffness matrix from a vector of stiffnesses (for MDOF systems)
    # Inputs:
    # k     --> vector of stiffnesses
    #           k = [k1, k2, ..., kn], where n is the number of floors
    # Outputs:
    # K     --> stiffness matrix 
    # Comments:
    # Note that the first element of the vector is the stiffness of the first floor, the second element is the stiffness of the first floor, and so on.
    if len(k) > 1:
        k_aux = k[1:]   # k_aux = [k2, k3, ..., kn]
        k_aux = np.append(k_aux, 0) # k_aux = [k2, k3, ..., kn, 0]
        K = np.diag(k + k_aux) - np.diag(k[1:], 1) - np.diag(k[1:], -1) # K = [k1+k2, -k2, 0, ..., 0; -k2, k2+k3, -k3, ..., 0; 0, ..., ..., -k(n-1), k(n-1)+kn]
    else:
        K = k # if there is only one floor, the stiffness matrix is just the stiffness of the floor
    return K
